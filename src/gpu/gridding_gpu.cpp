/* 
Header file for both gridding + cuFFT, implemented in cuda 
- gridding_imaging

References: 
https://stackoverflow.com/questions/17489017/can-we-declare-a-variable-of-type-cufftcomplex-in-side-a-kernel
*/

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>
#include <astroio.hpp>
#include <memory_buffer.hpp>

#include "pacer_imager_hip_defines.h"
#include "../gridding.hpp"


__device__ inline int wrap_index(int i, int side){
    if(i >= 0) return i % side;
    else return (side + i);
}


__device__ int calculate_pos(float u,
                             float v,
                             double delta_u,
                             double delta_v,
                             double wavelength,
                             double min_uv,
                             int n_pixels,
                             int uv_sign) // will be called with +1 and -1 
{
   // Operation 1: uv_lambda()
   double u_lambda = (u)/(wavelength); 
   double v_lambda = (v)/(wavelength); 

   // Calculating distance between the two antennas 
   double uv_distance = sqrt(u_lambda*u_lambda + v_lambda*v_lambda);

   if( uv_distance > min_uv )
   {            
      // (For all the rows of the Correlation Matrix)
      // Operation 2: uv_index()
      double u_pix, v_pix;
      u_pix = round(u_lambda/delta_u); 
      v_pix = round(v_lambda/delta_v);
      int u_index = wrap_index(uv_sign*u_pix, n_pixels); 
      int v_index = wrap_index(uv_sign*v_pix, n_pixels);
      return (n_pixels*v_index) + u_index; 
   }

   // same as else :   
   return -1;
}



__global__ void gridding_kernel(const float *visibilities, unsigned int n_baselines, unsigned int n_frequencies, unsigned int n_intervals,
                                      int n_ant,
                                      const float *u, const float *v, 
                                      const int* antenna_flags, const float* antenna_weights,
                                      const double *frequencies, int image_size, double delta_u, double delta_v, 
                                      int n_pixels,
                                      float *uv_grid_counter, double min_uv, Polarization pol,
                                      gpufftComplex *m_in_buffer) {

   unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
   unsigned int grid_size = gridDim.x * blockDim.x;
   unsigned int total_baselines = n_baselines * n_frequencies * n_intervals;
   const int n_pols_prod = 4;

   for (; i < total_baselines; i += grid_size){
      unsigned int baseline = i % n_baselines;
      unsigned int m_idx = i / n_baselines;
      unsigned int fine_channel = m_idx % n_frequencies;

      float re {0}, im {0};

      if(pol == Polarization::XX){
         re = visibilities[i * 2 * n_pols_prod];
         im = visibilities[i * 2 * n_pols_prod + 1];
      }else if(pol == Polarization::YY){
         re = visibilities[i * 2 * n_pols_prod + 6];
         im = visibilities[i * 2 * n_pols_prod + 7];

      
      }else {
         // Stokes I
         re = (visibilities[i * 2 * n_pols_prod] + visibilities[i * 2 * n_pols_prod + 6]) / 2.0f;
         im = (visibilities[i * 2 * n_pols_prod + 1] + visibilities[i * 2 * n_pols_prod + 7]) / 2.0f;
      }
      
      unsigned int a1 {static_cast<unsigned int>(-0.5 + sqrt(0.25 + 2*baseline))};
      unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
      if(a1 == a2) continue;

      // Checking for NaN values 
      if( !isnan(re) && !isnan(im) && antenna_flags[a2]<=0 && antenna_flags[a1]<=0 ) {
         int pos = calculate_pos( u[baseline], v[baseline], delta_u, delta_v, VEL_LIGHT / frequencies[fine_channel], min_uv, n_pixels,  +1 );
         if(pos>=0 && pos<image_size) {
            // Allocating in uv_grid       
            // WARNING: this might not give us the exact count for all time steps because of nan values, but it is a tradeoff for memory
            // Maybe over time, if not uniform weighting          
            if(m_idx / n_frequencies == 0) atomicAdd(&uv_grid_counter[image_size * m_idx + pos],1);

            // Allocating inside m_in_buffer as well 
            atomicAdd(&m_in_buffer[image_size * m_idx + pos].x,re);
            atomicAdd(&m_in_buffer[image_size * m_idx + pos].y,im);
         }   

         int pos2 = calculate_pos(u[baseline], v[baseline], delta_u, delta_v, VEL_LIGHT / frequencies[fine_channel], min_uv, n_pixels, -1 );
         if(pos2>=0 && pos2<image_size)
         {
            if(m_idx / n_frequencies == 0) atomicAdd(&uv_grid_counter[image_size * m_idx + pos2],1);
            // Allocating inside m_in_buffer as well 
            atomicAdd(&m_in_buffer[image_size * m_idx + pos2].x,re);
            atomicAdd(&m_in_buffer[image_size * m_idx + pos2].y,-im);
         }        
      }
   }
}





void gridding_gpu(const Visibilities& xcorr,
      const MemoryBuffer<float>& u_gpu,  const MemoryBuffer<float>& v_gpu, 
      const MemoryBuffer<int>& antenna_flags, const MemoryBuffer<float>& antenna_weights,
      const MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, Polarization pol, MemoryBuffer<float>& grids_counters,
      MemoryBuffer<std::complex<float>>& grids){

  int n_ant = xcorr.obsInfo.nAntennas;
  int image_size {n_pixels * n_pixels}; 

   size_t n_images {xcorr.integration_intervals() * xcorr.nFrequencies};
   size_t buffer_size {image_size * n_images};
      
   int n_baselines = ((xcorr.obsInfo.nAntennas + 1) * xcorr.obsInfo.nAntennas) / 2;
   struct gpuDeviceProp_t props;
   int gpu_id = -1;
   gpuGetDevice(&gpu_id);
   gpuGetDeviceProperties(&props, gpu_id);
   unsigned int n_blocks = props.multiProcessorCount * 2;
   gridding_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<const float*>(xcorr.data()), n_baselines, xcorr.nFrequencies, xcorr.integration_intervals(),
      n_ant, u_gpu.data(), v_gpu.data(), antenna_flags.data(), antenna_weights.data(), frequencies.data(), image_size,
      delta_u, delta_v, n_pixels, grids_counters.data(), min_uv, pol, (gpufftComplex*) grids.data());
   gpuGetLastError();
   gpuDeviceSynchronize();
}