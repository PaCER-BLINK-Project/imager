#ifndef _GRIDDING_GPU_H__
#define _GRIDDING_GPU_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>
#include <astroio.hpp>
#include <memory_buffer.hpp>




void gridding_gpu(Visibilities& xcorr, int time_step, int fine_channel,
      MemoryBuffer<float>& u_gpu,  MemoryBuffer<float>& v_gpu, 
      int* antenna_flags, float* antenna_weights,
       MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, MemoryBuffer<float>& grids_counters,
      MemoryBuffer<std::complex<float>>& grids);

#endif 
