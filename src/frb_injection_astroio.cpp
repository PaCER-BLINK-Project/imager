#include "frb_injection_astroio.hpp"
#include <gpu_macros.hpp>
#include <algorithm>
#include <stdexcept>
#include <iostream>
namespace blink_frb {
std::size_t inject(Visibilities& vis, const Catalogue& catalogue,
                   double ra, double dec, const MemoryBuffer<float>& u,
                   const MemoryBuffer<float>& v, const MemoryBuffer<float>& w,
                   const MemoryBuffer<double>& frequencies, const std::vector<int>& flags) {
    if(catalogue.empty()) return 0;
    if(std::none_of(catalogue.begin(),catalogue.end(),[&](const Event& e){
        return e.observation_id=="*" || e.observation_id==vis.obsInfo.id;
    })) return 0;
    const auto& info=vis.obsInfo;
    if(info.nPolarizations!=2 || !vis.nIntegrationSteps || !vis.nAveragedChannels ||
       !info.nAntennas || !info.nTimesteps || !vis.nFrequencies ||
       info.nFrequencies%vis.nAveragedChannels || info.nTimesteps%vis.nIntegrationSteps ||
       vis.nFrequencies!=info.nFrequencies/vis.nAveragedChannels)
        throw std::invalid_argument("FRB injection requires 4 products and complete time/frequency integrations");
    const std::size_t nb=std::size_t(info.nAntennas)*(std::size_t(info.nAntennas)+1)/2;
    if(!vis || !u || !v || !w || !frequencies || u.size()!=nb || v.size()!=nb || w.size()!=nb ||
       flags.size()!=nb || frequencies.size()<vis.nFrequencies || u.on_gpu() || v.on_gpu() || w.on_gpu() || frequencies.on_gpu() ||
       vis.matrix_size()!=4*nb || static_cast<const MemoryBuffer<std::complex<float>>&>(vis).size()<vis.size())
        throw std::invalid_argument("FRB injection: unexpected AstroIO layout/UVW");
    Geometry g{info.id,static_cast<long double>(info.startTime),info.timeResolution*vis.nIntegrationSteps,
       info.frequencyResolution*vis.nAveragedChannels*1e6,ra,dec,vis.integration_intervals(),
       std::vector<double>(frequencies.data(),frequencies.data()+vis.nFrequencies)};
    const auto plan=make_plan(catalogue,g);
    if(plan.empty()) return 0;
    // DEVICE data stays on the device. Transfer only active [time,channel] slices.
    // Synchronous transfers also order this work after GPU phase corrections.
    std::vector<std::complex<float>> staging;
    const bool device=vis.on_gpu();
    if(device) {
        #ifdef __GPU__
        gpuDeviceSynchronize();
        staging.resize(4*nb);
        #else
        throw std::runtime_error("FRB injection: device data in CPU-only compilation");
        #endif
    }
    for(const auto& [cell, contributions]:plan) {
        const auto t=cell/vis.nFrequencies, f=cell%vis.nFrequencies;
        auto* ptr=vis.at(static_cast<unsigned>(t),static_cast<unsigned>(f),0u);
        auto* host=ptr;
        #ifdef __GPU__
        if(device) {
            gpuMemcpy(staging.data(),ptr,staging.size()*sizeof(staging[0]),gpuMemcpyDeviceToHost);
            host=staging.data();
        }
        #endif
        add_slice(host,nb,u.data(),v.data(),w.data(),flags.data(),g.frequencies_hz[f],contributions);
        #ifdef __GPU__
        if(device) gpuMemcpy(ptr,staging.data(),staging.size()*sizeof(staging[0]),gpuMemcpyHostToDevice);
        #endif
    }
    std::clog << "FRB injection: obs=" << info.id << " start_unix=" << info.startTime
              << " coarse_channel=" << info.coarseChannel << " active_cells=" << plan.size() << '\n';
    return plan.size();
}
} // namespace blink_frb
