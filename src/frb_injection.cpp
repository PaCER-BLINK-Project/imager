#include "frb_injection.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace blink_frb {
namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::invalid_argument(message);
}
bool positive(double x) { return std::isfinite(x) && x > 0; }
double scattering_tau(const Event& e, double hz) {
    if (e.scattering_tau_s==0.0) return 0.0;
    const double tau=e.scattering_tau_s*std::pow(hz/e.scattering_ref_hz,e.scattering_index);
    require(positive(tau) && std::isfinite(frbsim::ScatteredBoxCar::duration(e.width_s,tau)),
            "FRB: scattering timescale/duration outside representable range");
    return tau;
}
double number(const std::string& s) {
    std::size_t n=0; const double v=std::stod(s,&n);
    require(n==s.size(),"FRB catalogue: invalid number"); return v;
}
int integer(const std::string& s) {
    std::size_t n=0; const int v=std::stoi(s,&n);
    require(n==s.size(),"FRB catalogue: invalid integer"); return v;
}
long double timestamp(const std::string& s) {
    std::size_t n=0; const auto v=std::stold(s,&n);
    require(n==s.size(),"FRB catalogue: invalid time"); return v;
}
}
void validate(const Event& e) {
    require(!e.id.empty() && !e.observation_id.empty(),"FRB: empty identifier");
    require(std::isfinite(e.ra_deg) && e.ra_deg>=0 && e.ra_deg<360,"FRB: RA outside [0,360)");
    require(std::isfinite(e.dec_deg) && std::abs(e.dec_deg)<=90,"FRB: invalid Dec");
    require(std::isfinite(e.arrival_unix_s),"FRB: invalid arrival time");
    require(std::isfinite(e.dm) && e.dm>=0,"FRB: invalid DM");
    require(positive(e.width_s) && positive(e.peak_jy) && positive(e.units_per_jy),"FRB: invalid width/peak/conversion");
    require(positive(e.ref_hz) || e.ref_hz==std::numeric_limits<double>::infinity(),"FRB: invalid reference frequency");
    require(std::isfinite(e.scattering_tau_s) && e.scattering_tau_s>=0,"FRB: invalid scattering timescale");
    require(positive(e.scattering_ref_hz),"FRB: invalid scattering reference frequency");
    require(std::isfinite(e.scattering_index),"FRB: invalid scattering index");
    require(e.phase_sign==1 || e.phase_sign==-1,"FRB: phase sign must be +1 or -1");
    require(e.tsub>=0 && e.tsub<=4096 && e.fsub>0 && e.fsub<=4096,"FRB: invalid subsampling counts");
    require(std::isfinite(e.peak_jy*e.units_per_jy) && e.peak_jy*e.units_per_jy <= std::numeric_limits<float>::max(),"FRB: amplitude outside float range");
}
Catalogue load_catalogue(const std::string& path) {
    std::ifstream in(path);
    if(!in) throw std::runtime_error("Cannot open FRB catalogue: "+path);
    Catalogue out; std::set<std::string> ids; std::string line; std::size_t lineno=0;
    while(std::getline(in,line)) {
        ++lineno; line=line.substr(0,line.find('#'));
        std::istringstream row(line); std::vector<std::string> fields; std::string field;
        while(row>>field) fields.push_back(field);
        if(fields.empty()) continue;
        try {
            require(fields.size()==13 || fields.size()==15 || fields.size()==16,
                    "expected 13 fields, or 15/16 with scattering_tau_s scattering_ref_mhz [scattering_index]");
            Event e{fields[0],fields[1],number(fields[2]),number(fields[3]),timestamp(fields[4]),
                number(fields[5]),number(fields[6]),number(fields[7]),number(fields[8])*1e6,
                number(fields[9]),integer(fields[10]),integer(fields[11]),integer(fields[12])};
            if(fields.size()>=15) {
                e.scattering_tau_s=number(fields[13]);
                e.scattering_ref_hz=number(fields[14])*1e6;
                if(fields.size()==16) e.scattering_index=number(fields[15]);
            }
            validate(e);
            require(ids.insert(e.id).second,"duplicate injection ID");
            out.push_back(e);
        } catch(const std::exception& ex) {
            throw std::runtime_error(path+":"+std::to_string(lineno)+": "+ex.what());
        }
    }
    if(in.bad()) throw std::runtime_error("Error reading FRB catalogue: "+path);
    return out;
}
Plan make_plan(const Catalogue& catalogue, const Geometry& g) {
    Plan plan;
    if(catalogue.empty()) return plan;
    require(std::isfinite(g.start_unix_s) && positive(g.dt_s) && positive(g.df_hz),"FRB: invalid time/frequency metadata");
    require(std::isfinite(g.phase_ra_deg) && std::isfinite(g.phase_dec_deg) && std::abs(g.phase_dec_deg)<=90,"FRB: invalid phase centre");
    const auto nf=g.frequencies_hz.size();
    require(nf>0 && g.nt>0 && g.nt<=std::numeric_limits<std::size_t>::max()/nf,"FRB: invalid visibility shape");
    for(double f:g.frequencies_hz) require(positive(f-g.df_hz/2),"FRB: invalid channel frequency/width");
    const double end=static_cast<double>(g.nt)*g.dt_s;
    require(std::isfinite(end),"FRB: invalid buffer duration");
    for(const auto& e:catalogue) {
        validate(e);
        if(e.observation_id!="*" && e.observation_id!=g.observation_id) continue;
        const double t0=static_cast<double>(e.arrival_unix_s-g.start_unix_s);
        require(std::isfinite(t0),"FRB: arrival time offset outside double range");
        frbsim::Source source{{t0,e.width_s,e.peak_jy,e.dm,e.ref_hz},e.ra_deg,e.dec_deg};
        const double ra=(e.ra_deg-g.phase_ra_deg)*frbsim::pi/180;
        const double dc=e.dec_deg*frbsim::pi/180, dc0=g.phase_dec_deg*frbsim::pi/180;
        const double n=std::sin(dc)*std::sin(dc0)+std::cos(dc)*std::cos(dc0)*std::cos(ra);
        if(n<=0) continue; // outside the front hemisphere; Python sqrt(n^2) loses this sign
        const auto lmn=source.lmn(g.phase_ra_deg,g.phase_dec_deg);
        for(std::size_t f=0;f<nf;++f) {
            const double fc=g.frequencies_hz[f];
            // Delay monotonic for nonnegative DM. Reject inactive channels before scanning time bins.
            const double lo=t0+source.lightcurve.delay(fc+g.df_hz/2);
            // The dispersion delay and scattering duration need not have the
            // same monotonicity. Bound the end using the largest tau at either edge.
            double hi=t0+source.lightcurve.delay(fc-g.df_hz/2)+e.width_s;
            if(e.scattering_tau_s>0) {
                const double tau_max=std::max(scattering_tau(e,fc-g.df_hz/2),
                                              scattering_tau(e,fc+g.df_hz/2));
                hi+=frbsim::ScatteredBoxCar::tail_efolds*tau_max;
            }
            require(std::isfinite(lo) && std::isfinite(hi),"FRB: nonfinite dispersed arrival");
            if(hi<=0 || lo>=end) continue;
            // Conservative range with a one-bin margin for floating-point rounding.
            const auto first=static_cast<std::size_t>(std::max(0.0,std::floor(lo/g.dt_s)-1));
            const auto last=static_cast<std::size_t>(std::min(static_cast<double>(g.nt),std::ceil(hi/g.dt_s)+1));
            struct Subchannel { double delay, tau; };
            std::vector<Subchannel> subchannels; subchannels.reserve(e.fsub);
            for(int k=0;k<e.fsub;++k) {
                const double fs=fc+((double(k)+0.5)/e.fsub-0.5)*g.df_hz;
                subchannels.push_back({source.lightcurve.delay(fs),scattering_tau(e,fs)});
            }
            for(std::size_t t=first;t<last;++t) {
                double fraction=0;
                const double begin=double(t)*g.dt_s, centre=(double(t)+0.5)*g.dt_s;
                for(const auto& sub:subchannels) {
                    const double d=sub.delay;
                    if(e.tsub==0) {
                        const double x=(begin-t0)-d;
                        fraction+=frbsim::ScatteredBoxCar::average(x,g.dt_s,e.width_s,sub.tau);
                    } else {
                        for(int k=0;k<e.tsub;++k) {
                            const double off=((double(k)+0.5)/e.tsub-0.5)*g.dt_s;
                            const double x=(centre+off-t0)-d;
                            if(sub.tau==0.0) {
                                if(x>=0 && x<=e.width_s) fraction+=1.0/e.tsub;
                            } else {
                                fraction+=frbsim::ScatteredBoxCar::evaluate(x,e.width_s,sub.tau)/e.tsub;
                            }
                        }
                    }
                }
                fraction=std::clamp(fraction/e.fsub,0.0,1.0);
                if(fraction>0) plan[t*nf+f].push_back({lmn,fraction*e.peak_jy*e.units_per_jy,e.phase_sign});
            }
        }
    }
    return plan;
}
void add_slice(std::complex<float>* slice,std::size_t nb,const float* u,const float* v,const float* w,
               const int* flags,double frequency_hz,const std::vector<Contribution>& contributions) {
    // Every iteration owns one baseline. Sources are summed in catalogue order.
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(nb>65536)
    #endif
    for(std::size_t b=0;b<nb;++b) {
        if(flags[b]) continue;
        frbsim::Cx delta{};
        for(const auto& c:contributions) {
            const auto ph=frbsim::phase_term({u[b],v[b],w[b]},c.lmn,frequency_hz);
            delta+=c.amplitude*(c.phase_sign==1?ph:std::conj(ph));
        }
        for(std::size_t p:{0u,3u}) {
            auto& z=slice[4*b+p];
            if(std::isfinite(z.real()) && std::isfinite(z.imag()))
                z=static_cast<std::complex<float>>(static_cast<frbsim::Cx>(z)+delta);
        }
    }
}
} // namespace blink_frb
