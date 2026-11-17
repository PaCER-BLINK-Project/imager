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
    auto bad=flags;bad.pop_back();bool threw=false;
    try{inject(vis,{e},92,-13,u,v,w,f,bad);}catch(const std::invalid_argument&){threw=true;}
    check(threw,"invalid shape rejected");
}
void test_catalogue(const std::string& filename){
    auto cat=load_catalogue(filename);check(cat.size()==1,"catalogue parse");
    check(cat[0].id=="demo"&&cat[0].ref_hz==155e6&&cat[0].tsub==0,"catalogue units");
    auto e=cat[0];e.dm=-1;bool threw=false;try{validate(e);}catch(const std::invalid_argument&){threw=true;}check(threw,"invalid event rejected");
}
int main(int argc,char** argv){try{
    check(argc>=2,"usage: frb_injection_test catalogue [--gpu]");
    test_plan();test_slice();test_catalogue(argv[1]);test_adapter(false);
    if(argc>2 && std::string(argv[2])=="--gpu"){
        #ifdef __GPU__
        test_adapter(true);
        #else
        throw std::runtime_error("GPU test requested in CPU build");
        #endif
    }
    std::cout<<"FRB injection tests passed\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
