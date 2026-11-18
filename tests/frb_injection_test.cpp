#include "frb_injection_astroio.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using namespace blink_frb;
using CF=std::complex<float>;
void check(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
void near(CF a,CF b,const char* msg){check(std::abs(a-b)<2e-5f,msg);}
Event event(){return {"obs","test",92.,-13.,1500000000.025L,0.,0.030,4.,150e6,2.,-1,0,10};}
void test_plan(){
    Geometry g{"obs",1500000000.L,0.020,10000,92.,-13.,5,{150e6,151e6}};
    auto e=event(); auto p=make_plan({e},g);
    check(p.size()==4,"active cells");
    // [25,55]ms overlaps [20,40] and [40,60] by 15ms each.
    for(auto& [key,v]:p){check(key>=2&&key<=5,"time/channel ordering");check(std::abs(v[0].amplitude-6)<1e-6,"exact time overlap");}
    e.tsub=20; e.dm=10; e.arrival_unix_s=g.start_unix_s-0.010L; e.ref_hz=151e6;
    frbsim::BoxCarBurst b{double(e.arrival_unix_s-g.start_unix_s),e.width_s,e.peak_jy,e.dm,e.ref_hz};
    std::vector<double> times;for(unsigned t=0;t<g.nt;++t)times.push_back((t+0.5)*g.dt_s);
    auto ref=b.sample(times,g.dt_s,g.frequencies_hz,g.df_hz,e.tsub,e.fsub);
    p=make_plan({e},g);
    for(size_t i=0;i<ref.size();++i){auto it=p.find(i);double val=it==p.end()?0:it->second[0].amplitude;
        check(std::abs(val-ref[i]*e.units_per_jy)<1e-10,"Python midpoint parity");}
    check(make_plan({e},Geometry{"different",g.start_unix_s,g.dt_s,g.df_hz,g.phase_ra_deg,g.phase_dec_deg,g.nt,g.frequencies_hz}).empty(),"observation filter");
    e.arrival_unix_s+=10000;check(make_plan({e},g).empty(),"inactive time");
    e=event();e.ra_deg=272;e.dec_deg=13;check(make_plan({e},g).empty(),"back hemisphere");
    e=event();auto p1=make_plan({e,e},g);check(p1.begin()->second.size()==2,"overlapping events");
    // Event straddles two buffers; total fluence preserved.
    e=event();e.arrival_unix_s=g.start_unix_s+0.095L; e.width_s=0.020;
    auto a=make_plan({e},g);auto g2=g;g2.start_unix_s+=0.1L;
    auto c=make_plan({e},g2);double sum=0;
    for(auto* q:{&a,&c})for(auto& [key,v]:*q)if(key%2==0)sum+=v[0].amplitude*g.dt_s;
    check(std::abs(sum-e.width_s*e.peak_jy*e.units_per_jy)<1e-7,"buffer-boundary fluence");
}
// Independent reference profile, evaluated by midpoint quadrature below.
long double reference_profile(long double x,long double width,long double tau){
    if(x<0 || x>=width+25*tau)return 0;
    return x<width ? 1-std::exp(-x/tau)
                   : std::exp(-(x-width)/tau)-std::exp(-x/tau);
}
double amplitude(const Plan& p,size_t key){
    const auto it=p.find(key);return it==p.end()?0.:it->second.at(0).amplitude;
}
void test_scattering(){
    using S=frbsim::ScatteredBoxCar;
    const double w=0.03,tau=0.05,dt=0.02;
    check(S::evaluate(-0.01,w,tau)==0,"causal profile");
    check(std::abs(S::evaluate(w,w,tau)-(1-std::exp(-w/tau)))<1e-15,"peak attenuation");
    check(std::abs(S::evaluate(w+tau,w,tau)/S::evaluate(w,w,tau)-std::exp(-1.0))<1e-15,"one e-fold");
    check(S::evaluate(w+26*tau,w,tau)==0,"tail cutoff");
    check(std::abs(S::average(0,1e-12,w,1.)-5e-13)<1e-24,"stable tiny rising integral");
    check(std::abs(S::average(0,dt,w,1e-15)-1)<1e-12,"small tau limit");
    for(double x:{-0.03,-0.005,0.,0.01,0.02,0.03,0.1,1.27,1.30}){
        long double sum=0;constexpr int n=100000;
        for(int k=0;k<n;++k)sum+=reference_profile(x+(k+.5L)*dt/n,w,tau);
        check(std::abs(S::average(x,dt,w,tau)-double(sum/n))<1e-9,"analytic average vs quadrature");
    }
    Geometry g{"obs",1500000000.L,dt,10000,92.,-13.,200,{150e6,200e6}};
    auto e=event();e.scattering_tau_s=tau;e.scattering_ref_hz=150e6;e.scattering_index=-4;
    e.fsub=3;
    auto p=make_plan({e},g);
    const double t0=double(e.arrival_unix_s-g.start_unix_s);
    for(size_t f=0;f<2;++f){
        double fluence=0;
        for(size_t t=0;t<g.nt;++t){
            double val=amplitude(p,t*2+f);fluence+=val*dt;
            long double expected=0; constexpr int n=2000;
            for(int k=0;k<e.fsub;++k){
                const long double hz=g.frequencies_hz[f]+((k+.5L)/e.fsub-.5L)*g.df_hz;
                const long double tau_f=tau*std::pow(hz/e.scattering_ref_hz,-4.L);
                for(int j=0;j<n;++j)expected+=reference_profile((t+(j+.5L)/n)*dt-t0,w,tau_f);
            }
            expected*=e.peak_jy*e.units_per_jy/(e.fsub*n);
            check(std::abs(val-double(expected))<2e-7,"frequency scaling/averaging vs quadrature");
        }
        check(std::abs(fluence-e.peak_jy*e.units_per_jy*w)<1e-10,"fluence conserved per channel");
    }
    check(amplitude(p,20)>amplitude(p,21),"longer low-frequency tail");
    // Midpoint path, including dispersion and a separate scattering reference.
    e.dm=10;e.ref_hz=200e6;e.tsub=20;auto m=make_plan({e},g);
    for(size_t t=0;t<g.nt;++t)for(size_t f=0;f<2;++f){
        long double expected=0;
        for(int k=0;k<e.fsub;++k){
            const long double hz=g.frequencies_hz[f]+((k+.5L)/e.fsub-.5L)*g.df_hz;
            const long double delay=frbsim::k_dm*e.dm*(std::pow(hz/1e9L,-2.L)-std::pow(e.ref_hz/1e9L,-2.L));
            const long double tau_f=tau*std::pow(hz/e.scattering_ref_hz,-4.L);
            for(int j=0;j<e.tsub;++j)expected+=reference_profile((t+(j+.5L)/e.tsub)*dt-t0-delay,w,tau_f);
        }
        expected*=e.peak_jy*e.units_per_jy/(e.fsub*e.tsub);
        check(std::abs(amplitude(m,t*2+f)-double(expected))<1e-10,"dispersed midpoint scattering");
    }
    // A later buffer begins after the unscattered boxcar ended. No tail lost.
    e.dm=0;e.tsub=0;g.start_unix_s=1500000000.L;g.dt_s=0.02;
    const auto full=make_plan({e},g);auto first=g;first.nt=5;
    auto second=g;second.nt=g.nt-first.nt;second.start_unix_s+=0.1L;
    const auto p1=make_plan({e},first),p2=make_plan({e},second);
    check(!p2.empty()&&amplitude(p2,0)>0,"tail survives next buffer");
    for(size_t t=0;t<g.nt;++t)for(size_t f=0;f<2;++f){
        const double val=t<first.nt?amplitude(p1,t*2+f):amplitude(p2,(t-first.nt)*2+f);
        check(std::abs(val-amplitude(full,t*2+f))<2e-8,"split buffer parity");
    }
    auto later=g;later.start_unix_s+=2.L;check(make_plan({e},later).empty(),"inactive after cutoff");
    // Positive indices are supported; the high edge then sets the tau bound.
    e.scattering_index=4;auto positive=make_plan({e},g);
    check(amplitude(positive,21)>amplitude(positive,20),"positive scattering index");
    for(int bad=0;bad<3;++bad){auto b=e;
        if(bad==0)b.scattering_tau_s=-1;
        if(bad==1)b.scattering_ref_hz=std::numeric_limits<double>::infinity();
        if(bad==2)b.scattering_index=std::numeric_limits<double>::quiet_NaN();
        bool threw=false;try{validate(b);}catch(const std::invalid_argument&){threw=true;}
        check(threw,"invalid scattering parameter rejected");
    }
}
void test_slice(){
    CF data[12];std::fill_n(data,12,CF{1,2});float u[]{0,100,120},v[]{0,30,40},w[]{0,20,10};int flags[]{1,0,1};
    const double f=150e6;const frbsim::LMN lmn{0.01,0.02,std::sqrt(1-0.0005)};
    const double phase=-2*frbsim::pi*(100*lmn[0]+30*lmn[1]+20*(lmn[2]-1))*f/frbsim::c_m_s;
    const CF delta=std::polar(3.f,float(phase));
    add_slice(data,3,u,v,w,flags,f,{{lmn,3.,-1}});
    near(data[4],CF(1,2)+delta,"XX negative phase");near(data[7],data[4],"YY");
    for(int i:{0,1,2,3,5,6,8,9,10,11})near(data[i],CF(1,2),"unchanged flag/crosshand");
    data[4]={std::numeric_limits<float>::quiet_NaN(),0};
    add_slice(data,3,u,v,w,flags,f,{{lmn,3.,1}});
    check(std::isnan(data[4].real()),"NaN unchanged");near(data[7],CF(1,2)+delta+std::conj(delta),"positive phase");
}
void test_adapter(bool gpu){
    ObservationInfo info{};info.nAntennas=3;info.nPolarizations=2;info.nFrequencies=4;info.nTimesteps=10;
    info.timeResolution=0.010;info.frequencyResolution=0.005;info.startTime=1500000000;info.id="obs";
    constexpr size_t nb=6, nf=2, nt=5, n=nb*4*nf*nt;
    MemoryBuffer<CF> raw(n);std::fill_n(raw.data(),n,CF{1,2});
    Visibilities vis(std::move(raw),info,2,2);
    MemoryBuffer<float> u(nb),v(nb),w(nb);MemoryBuffer<double> f(nf);
    std::fill_n(u.data(),nb,0.f);std::fill_n(v.data(),nb,0.f);std::fill_n(w.data(),nb,0.f);
    f[0]=150e6;f[1]=151e6;std::vector<int> flags{1,0,1,0,1,1};
    if(gpu)vis.to_gpu();
    check(inject(vis,{},92,-13,u,v,w,f,flags)==0,"empty catalogue");
    const auto e=event();check(inject(vis,{e},92,-13,u,v,w,f,flags)==4,"adapter active count");
    if(gpu){check(vis.on_gpu(),"device residency preserved");vis.to_cpu();}
    for(size_t t=0;t<nt;++t)for(size_t ch=0;ch<nf;++ch)for(size_t b=0;b<nb;++b)for(size_t p=0;p<4;++p){
        const bool active=(t==1||t==2)&&!flags[b]&&(p==0||p==3);
        near(vis.at(t,ch,b)[p],active?CF{7,2}:CF{1,2},"AstroIO indexing/scale");}
    // A second, scattered pulse follows the same visibility adapter path.
    auto scattered=e;scattered.scattering_tau_s=0.050;
    Geometry g{"obs",1500000000.L,0.020,10000,92.,-13.,nt,{150e6,151e6}};
    const auto plan=make_plan({scattered},g);
    check(amplitude(plan,8)>0,"adapter late tail planned");
    if(gpu)vis.to_gpu();
    check(inject(vis,{scattered},92,-13,u,v,w,f,flags)==plan.size(),"adapter scattering active count");
    if(gpu){check(vis.on_gpu(),"scattered device residency preserved");vis.to_cpu();}
    for(size_t t=0;t<nt;++t)for(size_t ch=0;ch<nf;++ch)for(size_t b=0;b<nb;++b)for(size_t p=0;p<4;++p){
        const bool parallel=!flags[b]&&(p==0||p==3);
        const float original=(parallel&&(t==1||t==2))?7.f:1.f;
        const float delta=parallel?float(amplitude(plan,t*nf+ch)):0.f;
        near(vis.at(t,ch,b)[p],CF{original+delta,2},"scattered AstroIO injection");}
    auto bad=flags;bad.pop_back();bool threw=false;
    try{inject(vis,{e},92,-13,u,v,w,f,bad);}catch(const std::invalid_argument&){threw=true;}
    check(threw,"invalid shape rejected");
}
void test_catalogue(const std::string& filename){
    auto cat=load_catalogue(filename);check(cat.size()==3,"catalogue parse");
    check(cat[0].id=="demo"&&cat[0].ref_hz==155e6&&cat[0].tsub==0,"catalogue units");
    check(cat[0].scattering_tau_s==0,"legacy catalogue defaults");
    check(cat[1].scattering_tau_s==0.050&&cat[1].scattering_ref_hz==150e6&&cat[1].scattering_index==-4.,"default scattering index");
    check(cat[2].scattering_tau_s==0.120&&cat[2].scattering_ref_hz==155e6&&cat[2].scattering_index==-4.4,"explicit scattering index");
    auto e=cat[0];e.dm=-1;bool threw=false;try{validate(e);}catch(const std::invalid_argument&){threw=true;}check(threw,"invalid event rejected");
}
int main(int argc,char** argv){try{
    check(argc>=2,"usage: frb_injection_test catalogue [--gpu]");
    test_plan();test_scattering();test_slice();test_catalogue(argv[1]);test_adapter(false);
    if(argc>2 && std::string(argv[2])=="--gpu"){
        #ifdef __GPU__
        test_adapter(true);
        #else
        throw std::runtime_error("GPU test requested in CPU build");
        #endif
    }
    std::cout<<"FRB injection tests passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
