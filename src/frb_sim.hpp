#pragma once
// A standalone C++17 port of the pure numerical parts of
// Andywang201605/visvar/mwa/python/vis_frb_sim.py (2026-11-17).
// Deliberately independent of BLINK, CASA, GPU frameworks and NumPy.
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace frbsim {
constexpr double c_m_s = 299792458.0;
constexpr double k_dm = 4.148808e-3; // s GHz^2 / (pc cm^-3)
constexpr double pi = 3.141592653589793238462643383279502884;
using Cx = std::complex<double>;
using UVW = std::array<double, 3>;  // metres
using LMN = std::array<double, 3>;
enum class Correlations { XX_YY, XX_XY_YX_YY, I };
inline std::size_t ncor(Correlations corr) {
    return corr == Correlations::XX_XY_YX_YY ? 4 :
           corr == Correlations::XX_YY ? 2 : 1;
}

struct BoxCarBurst {
    double t0, width, peak, dm;
    double ref_frequency_hz = std::numeric_limits<double>::infinity();

    double delay(double frequency_hz) const {
        // Same expression/order as Python, including infinite ref frequency.
        double d = k_dm * dm * std::pow(frequency_hz / 1e9, -2.0);
        d -= k_dm * dm * std::pow(ref_frequency_hz / 1e9, -2.0);
        return d;
    }
    double evaluate(double time_s, double frequency_hz) const {
        const double x = time_s - t0 - delay(frequency_hz);
        return (x >= 0.0 && x <= width) ? peak : 0.0; // inclusive edges
    }
    // Output [time, channel]. Mean of a rectangular tsub x fsub grid.
    std::vector<double> sample(const std::vector<double>& times_s, double dt_s,
                               const std::vector<double>& frequencies_hz,
                               double df_hz, int tsub=10, int fsub=10) const {
        if (tsub <= 0 || fsub <= 0) throw std::invalid_argument("subsample counts must be positive");
        std::vector<double> out(times_s.size()*frequencies_hz.size());
        for (std::size_t t=0; t<times_s.size(); ++t) {
            for (std::size_t f=0; f<frequencies_hz.size(); ++f) {
                double total = 0.0;
                for (int ti=0; ti<tsub; ++ti) {
                    const double toff = ((static_cast<double>(ti)+0.5)/tsub - 0.5) * dt_s;
                    for (int fi=0; fi<fsub; ++fi) {
                        const double foff = ((static_cast<double>(fi)+0.5)/fsub - 0.5) * df_hz;
                        total += evaluate(times_s[t]+toff, frequencies_hz[f]+foff);
                    }
                }
                out[t*frequencies_hz.size()+f] = total / (tsub*fsub);
            }
        }
        return out;
    }
};

// Unit-height intrinsic boxcar convolved with a normalized, causal
// exponential. x is time since the dispersed leading edge (seconds).
// Values are relative to the intrinsic peak; the fluence is width (seconds).
// Truncate after 25 e-foldings beyond the intrinsic end to keep sparse plans
// finite. The omitted fraction of intrinsic fluence is <= exp(-25).
struct ScatteredBoxCar {
    static constexpr double tail_efolds = 25.0;
    static double duration(double width, double tau) {
        return width + tail_efolds*tau;
    }
    static double evaluate(double x, double width, double tau) {
        if (x < 0.0) return 0.0;
        if (tau == 0.0) return x <= width ? 1.0 : 0.0;
        if (x >= duration(width,tau)) return 0.0;
        if (x <= width) return -std::expm1(-x/tau);
        return -std::expm1(-width/tau)*std::exp(-(x-width)/tau);
    }
    // Mean over [x, x+dt]. Split at the end of the intrinsic boxcar.
    // Integrate each part directly to avoid cancellation between CDFs.
    static double average(double x, double dt, double width, double tau) {
        if (tau == 0.0)
            return std::fmax(0.0,std::fmin(x+dt,width)-std::fmax(x,0.0))/dt;
        const double a=std::fmax(0.0,x), b=std::fmin(x+dt,duration(width,tau));
        if (b <= a) return 0.0;
        double area=0.0;
        if (a < width) {
            const double h=std::fmin(b,width)-a, z=h/tau;
            // 1 - (1-exp(-z))/z; series prevents cancellation for small z.
            const double loss=z<1e-3 ? z*(0.5+z*(-1.0/6+z*(1.0/24+z*(-1.0/120))))
                                     : 1.0+std::expm1(-z)/z;
            area+=h*(-std::expm1(-a/tau)+std::exp(-a/tau)*loss);
        }
        if (b > width) {
            const double start=std::fmax(a,width), h=b-start, z=h/tau;
            const double integ=z>1.0 ? -tau*std::expm1(-z)
                                     : h*(z==0.0 ? 1.0 : -std::expm1(-z)/z);
            area+=-std::expm1(-width/tau)*std::exp(-(start-width)/tau)*integ;
        }
        return area/dt;
    }
};

struct Source {
    BoxCarBurst lightcurve;
    double ra_deg, dec_deg;

    LMN lmn(double phase_ra_deg, double phase_dec_deg) const {
        const double rad = pi/180.0;
        const double ra = ra_deg*rad, dec = dec_deg*rad;
        const double ra0 = phase_ra_deg*rad, dec0 = phase_dec_deg*rad;
        const double dra = ra-ra0;
        const double l = std::cos(dec)*std::sin(dra);
        const double m = std::sin(dec)*std::cos(dec0) - std::cos(dec)*std::sin(dec0)*std::cos(dra);
        const double n = std::sqrt(std::fmax(1.0-l*l-m*m, 0.0));
        return {l,m,n};
    }
    // Output [time, channel, correlation]; convention matches the Python source.
    std::vector<Cx> dynamic_spectrum(const std::vector<double>& times_s, double dt_s,
                       const std::vector<double>& frequencies_hz, double df_hz,
                       Correlations corrs, int tsub=10, int fsub=10) const {
        const auto intensity = lightcurve.sample(times_s,dt_s,frequencies_hz,df_hz,tsub,fsub);
        const std::size_t np = ncor(corrs);
        std::vector<Cx> out(intensity.size()*np, Cx{0.,0.});
        for (std::size_t i=0; i<intensity.size(); ++i) {
            out[i*np] = Cx(intensity[i], 0.0);
            if (corrs == Correlations::XX_YY) out[i*np+1] = Cx(intensity[i],0.0);
            if (corrs == Correlations::XX_XY_YX_YY) out[i*np+3] = Cx(intensity[i],0.0);
        }
        return out;
    }
};

inline Cx phase_term(const UVW& uvw_m, const LMN& lmn, double frequency_hz) {
    const double path = uvw_m[0]*lmn[0]+uvw_m[1]*lmn[1]+uvw_m[2]*(lmn[2]-1.0);
    const double angle = 2.0*pi*path*(frequency_hz/c_m_s);
    return std::exp(Cx(0.0,angle)); // Python uses positive +2j exponent
}

// CASA-like rows: times.size() == uvw_m.size(). Output [row,chan,pol].
inline std::vector<Cx> predict_tabrow(const Source& source, double ra0_deg, double dec0_deg,
                const std::vector<double>& times_s, double dt_s,
                const std::vector<UVW>& uvw_m, const std::vector<double>& frequencies_hz,
                double df_hz, Correlations corrs, int tsub=10, int fsub=10) {
    if (times_s.size()!=uvw_m.size()) throw std::invalid_argument("time and UVW row counts differ");
    const std::size_t nf=frequencies_hz.size(), np=ncor(corrs);
    const LMN coords=source.lmn(ra0_deg,dec0_deg);
    auto out=source.dynamic_spectrum(times_s,dt_s,frequencies_hz,df_hz,corrs,tsub,fsub);
    for(std::size_t row=0; row<times_s.size(); ++row)
        for(std::size_t f=0; f<nf; ++f) {
            const Cx phase=phase_term(uvw_m[row],coords,frequencies_hz[f]);
            for(std::size_t p=0; p<np; ++p)
                out[(row*nf+f)*np+p]*=phase;
        }
    return out;
}

// The Python 'predict_visibilities' has inconsistent baseline/time ordering
// when fixuvw=true: repeated UVWs are baseline-major but its time-index array
// is time-major. Set legacy_fixuvw_order=true only for exact Python comparison.
// For corrected mode, both are time-major.
inline std::vector<Cx> predict_regular(const Source& source, double ra0_deg, double dec0_deg,
                const std::vector<double>& times_s, double dt_s,
                const std::vector<UVW>& uvw_m, const std::vector<double>& frequencies_hz,
                double df_hz, Correlations corrs, int tsub=10, int fsub=10,
                bool fixuvw=true, bool legacy_fixuvw_order=false) {
    const std::size_t nt=times_s.size();
    if(nt==0) throw std::invalid_argument("need at least one time");
    if(!fixuvw && (uvw_m.size()%nt!=0)) throw std::invalid_argument("UVW rows must divide by ntime");
    const std::size_t nbl=fixuvw ? uvw_m.size():uvw_m.size()/nt;
    const std::size_t nf=frequencies_hz.size(), np=ncor(corrs);
    auto src=source.dynamic_spectrum(times_s,dt_s,frequencies_hz,df_hz,corrs,tsub,fsub);
    const LMN coords=source.lmn(ra0_deg,dec0_deg);
    std::vector<Cx> out(nt*nbl*nf*np);
    for(std::size_t row=0;row<nt*nbl;++row) {
        const std::size_t ti = row/nbl;
        const UVW& uvw = fixuvw ? uvw_m[legacy_fixuvw_order ? (row/nt) : (row%nbl)] : uvw_m[row];
        for(std::size_t f=0;f<nf;++f){
            const Cx phase=phase_term(uvw,coords,frequencies_hz[f]);
            for(std::size_t p=0;p<np;++p)
                out[(row*nf+f)*np+p]=src[(ti*nf+f)*np+p]*phase;
        }
    }
    return out;
}

// Float storage matches the actual visibility buffer; no writes to source input.
inline void add_inplace(std::vector<std::complex<float>>& destination,
                        const std::vector<Cx>& injection) {
    if (destination.size()!=injection.size()) throw std::invalid_argument("visibility buffer shape mismatch");
    for(std::size_t i=0;i<destination.size();++i)
        destination[i]=static_cast<std::complex<float>>(static_cast<Cx>(destination[i])+injection[i]);
}
} // namespace frbsim
