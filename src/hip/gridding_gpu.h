#ifndef _GRIDDING_GPU_H__
#define _GRIDDING_GPU_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>
#include <astroio.hpp>
#include <bg_fits.h>
#include <memory_buffer.hpp>

// same as above but using xcorr structure instead :
__global__ void gridding_imaging_cuda_xcorr( int xySize, // size of the correlation matrix
                                      int n_ant,
                                      float *u_cuda, float *v_cuda, 
                                      int* antenna_flags, float* antenna_weights,
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda,
                                      float *vis_cuda,  
                                      float *uv_grid_counter_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda);




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

__device__ int calculate_pos(float u,
                             float v,
                             double delta_u_cuda,
                             double delta_v_cuda,
                             double wavelength_cuda,
                             double min_uv_cuda,
                             int n_pixels_cuda,
                             int center_x_cuda,
                             int center_y_cuda,
                             int is_odd_x_cuda,
                             int is_odd_y_cuda,
                             int uv_sign_cuda);



void gridding_gpu(Visibilities& xcorr, int time_step, int fine_channel,
      CBgFits& fits_vis_u, CBgFits& fits_vis_v,
      int* antenna_flags, float* antenna_weights,
       MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, MemoryBuffer<float>& grids_counters_buffer,
      MemoryBuffer<std::complex<float>>& grids_buffer);

#endif 
