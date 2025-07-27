#ifndef _GRIDDING_GPU_H__
#define _GRIDDING_GPU_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>
#include <astroio.hpp>
#include <memory_buffer.hpp>




void gridding_gpu(const Visibilities& xcorr,
      const MemoryBuffer<float>& u_gpu,  const MemoryBuffer<float>& v_gpu, 
      const MemoryBuffer<int>& antenna_flags, const MemoryBuffer<float>& antenna_weights,
      const MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, MemoryBuffer<float>& grids_counters,
      MemoryBuffer<std::complex<float>>& grids);

#endif 
