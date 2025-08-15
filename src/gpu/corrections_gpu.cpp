#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <gpu_macros.hpp>
#include "pacer_imager_hip_defines.h"


__global__ void apply_cable_corrections(float *visibilities, unsigned int n_baselines, unsigned int n_frequencies, unsigned int n_intervals, double *cable_lengths,
      double* frequencies, double speed_of_light) {
   
   unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
   unsigned int grid_size = gridDim.x * blockDim.x;
   unsigned int total_baselines = n_baselines * n_frequencies * n_intervals;
   const int n_pols_prod = 4;

   for (; i < total_baselines; i += grid_size){
      unsigned int baseline = i % n_baselines;
      unsigned int m_idx = i / n_baselines;
      unsigned int fine_channel = m_idx % n_frequencies;

      for(int pol_idx {0}; pol_idx < n_pols_prod; pol_idx++){
         double re = visibilities[i * 2 * n_pols_prod + 2*pol_idx];
         double im = visibilities[i * 2 * n_pols_prod + 2*pol_idx + 1];

         unsigned int a1 {static_cast<unsigned int>(-0.5 + sqrt(0.25 + 2*baseline))};
         unsigned int a2 {baseline - ((a1 + 1) * a1)/2};

         double cableDeltaLen = (cable_lengths[a2] - cable_lengths[a1]);
         double angle = -2.0*M_PI*cableDeltaLen*frequencies[fine_channel] / speed_of_light;
         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;

         visibilities[i * 2 * n_pols_prod+ 2*pol_idx] = re_prim;
         visibilities[i * 2 * n_pols_prod + 2*pol_idx + 1] = im_prim;
      }
   }
}

__global__ void apply_geometric_corrections(float *visibilities, unsigned int n_baselines, unsigned int n_frequencies, unsigned int n_intervals, unsigned int n_ant, float *w,
      double* frequencies, double speed_of_light) {
   
   unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
   unsigned int grid_size = gridDim.x * blockDim.x;
   unsigned int total_baselines = n_baselines * n_frequencies * n_intervals;
   const int n_pols_prod = 4;

   for (; i < total_baselines; i += grid_size){
      unsigned int baseline = i % n_baselines;
      unsigned int m_idx = i / n_baselines;
      unsigned int fine_channel = m_idx % n_frequencies;
      
      for(int pol_idx {0}; pol_idx < n_pols_prod; pol_idx++){
         double re = visibilities[i * 2 * n_pols_prod + 2*pol_idx];
         double im = visibilities[i * 2 * n_pols_prod + 2*pol_idx + 1];
         
         unsigned int a1 {static_cast<unsigned int>(-0.5 + sqrt(0.25 + 2*baseline))};
         unsigned int a2 {baseline - ((a1 + 1) * a1)/2};
         // TODO this will have to change once we save w on a lower triangular matrix basis
         double angle = 2.0 * M_PI *  w[a1 * n_ant + a2] * frequencies[fine_channel] / speed_of_light;
         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;
         
         visibilities[i * 2 * n_pols_prod + 2*pol_idx] = re_prim;
         visibilities[i * 2 * n_pols_prod + 2*pol_idx + 1] = im_prim;
      }
   }
}


void apply_cable_lengths_corrections_gpu(Visibilities &xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   cable_lengths.to_gpu();
   frequencies.to_gpu();
   int n_baselines = (xcorr.obsInfo.nAntennas + 1) * (xcorr.obsInfo.nAntennas / 2);
   struct gpuDeviceProp_t props;
   int gpu_id = -1;
   gpuGetDevice(&gpu_id);
   gpuGetDeviceProperties(&props, gpu_id);
   unsigned int n_blocks = props.multiProcessorCount * 2;
   apply_cable_corrections<<<n_blocks, NTHREADS>>>(reinterpret_cast<float*>(xcorr.data()), n_baselines, xcorr.nFrequencies, xcorr.integration_intervals(), cable_lengths.data(), frequencies.data(), SPEED_OF_LIGHT);
   gpuCheckLastError();
}


void apply_geometric_corrections_gpu(Visibilities &xcorr, float *w_gpu, MemoryBuffer<double>& frequencies){
   xcorr.to_gpu();
   frequencies.to_gpu();
   int n_baselines = (xcorr.obsInfo.nAntennas + 1) * (xcorr.obsInfo.nAntennas / 2);
   struct gpuDeviceProp_t props;
   int gpu_id = -1;
   gpuGetDevice(&gpu_id);
   gpuGetDeviceProperties(&props, gpu_id);
   unsigned int n_blocks = props.multiProcessorCount * 2;
   apply_geometric_corrections<<<n_blocks, NTHREADS>>>(reinterpret_cast<float*>(xcorr.data()), n_baselines, xcorr.nFrequencies, xcorr.integration_intervals(), xcorr.obsInfo.nAntennas, w_gpu, frequencies.data(), SPEED_OF_LIGHT);
   gpuCheckLastError();
}
