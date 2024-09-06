#ifndef _GRIDDING_GPU_H__
#define _GRIDDING_GPU_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>
#include <astroio.hpp>
#include <bg_fits.h>
#include <memory_buffer.hpp>




void gridding_gpu(Visibilities& xcorr, int time_step, int fine_channel,
      CBgFits& fits_vis_u, CBgFits& fits_vis_v,
      int* antenna_flags, float* antenna_weights,
       MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, MemoryBuffer<float>& grids_counters_buffer,
      MemoryBuffer<std::complex<float>>& grids_buffer);

#endif 
