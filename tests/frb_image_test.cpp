// Uses BLINK's existing observation fixture (BLINK_TEST_DATADIR).
#include "pacer_imager.h"
#ifdef __GPU__
#include "gpu/pacer_imager_hip.h"
#endif
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

void run_test(CPacerImager& imager, ObservationInfo info, bool gpu){
    info.nTimesteps=200; // one 20-ms integration for VCS
    const auto n=std::size_t(info.nAntennas)*(info.nAntennas+1)/2*4;
    auto blank=[&](){
        MemoryBuffer<std::complex<float>> buf(n);
        std::fill_n(buf.data(),n,std::complex<float>{0,0});
        Visibilities vis(std::move(buf),info,info.nTimesteps,info.nFrequencies);
        if(gpu)vis.to_gpu();return vis;
    };
    imager.apply_geom_correction=false;imager.apply_cable_correction=false;
    imager.set_frb_injections({}); // override environment and verify no-op
    auto vis0=blank();auto img0=imager.run(vis0);img0.to_cpu();
    // Reuse actual phase centre resolved by grid()/metadata, not a guessed sky location.
    const double ra=std::fmod(img0.ra_deg+360.,360.);
    const double dt=info.timeResolution*info.nTimesteps;
    blink_frb::Event event{"*","centre",ra,img0.dec_deg,static_cast<long double>(info.startTime),
        0.,2*dt,4.,155e6,1.,-1,0,10};
    imager.set_frb_injections({event});
    auto vis1=blank();auto img1=imager.run(vis1);img1.to_cpu();
    const auto side=img1.side_size;const auto centre=(side/2)*side+side/2;
    const auto* orig=img0.at(0,0);const auto* inj=img1.at(0,0);
    const auto delta=inj[centre]-orig[centre];
    if(!std::isfinite(delta.real()) || std::abs(delta.real()-4.)>0.02 || std::abs(delta.imag())>0.02)
        throw std::runtime_error("Phase-centre image test: wrong amplitude");
    for(size_t i=0;i<img1.image_size();++i){
        const auto d=inj[i]-orig[i];
        if(!std::isfinite(d.real()) || d.real()>delta.real()+0.02)
            throw std::runtime_error("Phase-centre image test: invalid or brighter off-centre peak");
    }
    std::cout<<"Phase-centre injection recovered "<<delta.real()<<" visibility units\n";
}
int main(int argc,char** argv){try{
    const char* root=std::getenv("BLINK_TEST_DATADIR");if(!root){std::cerr<<"Set BLINK_TEST_DATADIR for image test\n";return 77;}
    const std::string base=std::string(root)+"/mwa/1276619416/";
    auto fixture=Visibilities::from_fits_file(base+"imager_stages/1s_ch000/input_visibilities.fits",VCS_OBSERVATION_INFO);
    auto info=fixture.obsInfo;
    const std::string meta=base+"20200619163000.metafits",weighting="N";
    std::vector<int> flags{21,25,58,71,80,81,92,101,108,114,119,125};
    const bool gpu=argc>1&&std::string(argv[1])=="--gpu";
    if(gpu){
#ifdef __GPU__
        CPacerImagerHip imager(meta,256,flags,false,Polarization::I,2.f,-1000.,std::numeric_limits<double>::infinity(),weighting.c_str());
        run_test(imager,info,true);
#else
        throw std::runtime_error("GPU image test requested in CPU build");
#endif
    }else{
        CPacerImager imager(meta,256,flags,false,Polarization::I,2.f,-1000.,std::numeric_limits<double>::infinity(),weighting.c_str());
        run_test(imager,info,false);
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
