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
#include <bg_fits.h>

#include "pacer_imager_hip_defines.h"

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

 
     
      // TODO : understand why this is giving wrong image with a black stripe in the centre !
      // This may optimise this code in the future (remove if-s) if it also produces the same results
      // WARNING : does not support ODD image sizes (add ASSERT !!!)
      // TODO : maybe copysign(int,int) is required - doesn't copysign use if-s too ?
      // int x_grid = round(u_index + copysignf( 1.0, float(center_x-u_index))*center_x);
      // int y_grid = round(v_index + copysignf( 1.0, float(center_y-v_index))*center_y);
      

      // Operation 4: Assignment of (re,im)vis to uv_grid
      // Position for assignment 
      return (n_pixels*v_index) + u_index; 
   }

   // same as else :   
   return -1;
}



__global__ void gridding_kernel(float *visibilities, unsigned int n_baselines, unsigned int n_frequencies, unsigned int n_intervals,
                                      int n_ant,
                                      float *u, float *v, 
                                      int* antenna_flags, float* antenna_weights,
                                      double *frequencies, int image_size, double delta_u, double delta_v, 
                                      int n_pixels,
                                      float *uv_grid_counter, double min_uv, 
                                      gpufftComplex *m_in_buffer) {

   unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
   unsigned int grid_size = gridDim.x * blockDim.x;
   unsigned int total_baselines = n_baselines * n_frequencies * n_intervals;
   const int n_pols_prod = 4;

   for (; i < total_baselines; i += grid_size){
      unsigned int baseline = i % n_baselines;
      unsigned int m_idx = i / n_baselines;
      unsigned int fine_channel = m_idx % n_frequencies;

      double re = visibilities[i * 2 * n_pols_prod];
      double im = visibilities[i * 2 * n_pols_prod + 1];
      
      unsigned int a1 {static_cast<unsigned int>(-0.5 + sqrt(0.25 + 2*baseline))};
      unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
      if(a1 == a2) continue;

      // Checking for NaN values 
      if( !isnan(re) && !isnan(im) && antenna_flags[a2]<=0 && antenna_flags[a1]<=0 ) {
         int pos = calculate_pos( u[a1 * n_ant + a2], v[a1 * n_ant + a2], delta_u, delta_v, VEL_LIGHT / frequencies[fine_channel], min_uv, n_pixels,  +1 );
         if(pos>=0 && pos<image_size) {
            // Allocating in uv_grid                
            atomicAdd(&uv_grid_counter[image_size * m_idx + pos],1);

            // Allocating inside m_in_buffer as well 
            atomicAdd(&m_in_buffer[image_size * m_idx + pos].x,re);
            atomicAdd(&m_in_buffer[image_size * m_idx + pos].y,im);
         }   

         int pos2 = calculate_pos(u[a1 * n_ant + a2], v[a1 * n_ant + a2], delta_u, delta_v, VEL_LIGHT / frequencies[fine_channel], min_uv, n_pixels, -1 );
         if(pos2>=0 && pos2<image_size)
         {
            atomicAdd(&uv_grid_counter[image_size * m_idx + pos2],1);
            // Allocating inside m_in_buffer as well 
            atomicAdd(&m_in_buffer[image_size * m_idx + pos2].x,re);
            atomicAdd(&m_in_buffer[image_size * m_idx + pos2].y,-im);
         }        
      }
   }
}





void gridding_gpu(Visibilities& xcorr, int time_step, int fine_channel,
      CBgFits& fits_vis_u, CBgFits& fits_vis_v,
      int* antenna_flags, float* antenna_weights,
      MemoryBuffer<double>& frequencies,
      double delta_u, double delta_v,
      int n_pixels, double min_uv, MemoryBuffer<float>& grids_counters,
      MemoryBuffer<std::complex<float>>& grids){
  std::cout << "Running 'gridding' on GPU.." << std::endl;

  int n_ant = xcorr.obsInfo.nAntennas;
  int corr_size = n_ant * n_ant;
  int image_size {n_pixels * n_pixels}; 

  float *u_cpu = fits_vis_u.get_data();
  float *v_cpu = fits_vis_v.get_data();
   float *u_gpu, *v_gpu;
   gpuMalloc((void**)&u_gpu, corr_size*sizeof(float));
   gpuMalloc((void**)&v_gpu, corr_size*sizeof(float));
   
   gpuMemcpy((float*)u_gpu, (float*)u_cpu, sizeof(float)*corr_size, gpuMemcpyHostToDevice); 
   gpuMemcpy((float*)v_gpu, (float*)v_cpu, sizeof(float)*corr_size, gpuMemcpyHostToDevice);
  


   size_t n_images {xcorr.integration_intervals() * xcorr.nFrequencies};
   size_t buffer_size {image_size * n_images};
   
   gpuMemset(grids_counters.data(), 0, n_images * image_size * sizeof(float));
   gpuMemset(grids.data(), 0, n_images * image_size * sizeof(std::complex<float>));
   
   int n_baselines = (xcorr.obsInfo.nAntennas + 1) * (xcorr.obsInfo.nAntennas / 2);
   struct gpuDeviceProp_t props;
   int gpu_id = -1;
   gpuGetDevice(&gpu_id);
   gpuGetDeviceProperties(&props, gpu_id);
   unsigned int n_blocks = props.multiProcessorCount * 2;

   frequencies.to_gpu();
   gridding_kernel<<<n_blocks, NTHREADS>>>(reinterpret_cast<float*>(xcorr.data()), n_baselines, xcorr.nFrequencies, xcorr.integration_intervals(),
      n_ant, u_gpu, v_gpu, antenna_flags, antenna_weights, frequencies.data(), image_size,
      delta_u, delta_v, n_pixels, grids_counters.data(), min_uv, (gpufftComplex*) grids.data());
   gpuGetLastError();
   gpuFree(u_gpu); 
   gpuFree(v_gpu);
}