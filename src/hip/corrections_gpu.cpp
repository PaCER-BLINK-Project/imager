#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <gpu_macros.hpp>
#include "gridding_imaging_cuda.h"
#include "pacer_imager_hip_defines.h"

void apply_cable_lengths_corrections_gpu(Visibilities &xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   cable_lengths.to_gpu();
   frequencies.to_gpu();

   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         double frequency_hz = frequencies[fine_channel];
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_cable_corrections<<<nBlocks, NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, cable_lengths.data(), frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
      }
   }
}


void apply_geometric_corrections_gpu(Visibilities &xcorr, float *w_gpu, MemoryBuffer<double>& frequencies){
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   frequencies.to_gpu();
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         double frequency_hz = frequencies[fine_channel];
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_geometric_corrections<<<nBlocks,NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, w_gpu, frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
      }
   }
}