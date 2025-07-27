#include "pacer_imager_hip.h"
#include "pacer_imager_hip_defines.h"

#include "gridding_gpu.h"

#define __GPU__
#include <gpu_macros.hpp>
#include <gpu_fft.hpp>

#include "../pacer_common.h"
#include "gpu_utils.h"
#include <mystring.h>
#include <mydate.h>
#include <exception>
#include "corrections_gpu.h"

#include "../utils.h"

CPacerImagerHip::CPacerImagerHip(const std::string metadata_file, const std::vector<int>& flagged_antennas) : CPacerImager(metadata_file, flagged_antennas) {}

void CPacerImagerHip::UpdateAntennaFlags(int n_ant) {
   if(!antenna_flags_gpu){
       antenna_flags_gpu.allocate(n_ant);      
      antenna_weights_gpu.allocate(n_ant);
      int n_flagged=0;   
      for(int ant=0;ant<n_ant;ant++){
         InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant]; 
         antenna_flags_gpu[ant] = ant1_info.flag;
         antenna_weights_gpu[ant] = 1.00;
         if( ant1_info.flag ){
            antenna_weights_gpu[ant] = 0.00;
            n_flagged++;
         }
      }
      antenna_flags_gpu.to_gpu();
      antenna_weights_gpu.to_gpu();
   PRINTF_INFO("INFO : CPacerImagerHip::UpdateAntennaFlags there are %d flagged antennas passed to GPU imager\n", n_flagged);
   }
}


// TODO : 
//     - do more cleanup of this function as this is nearly "RAW" copy paste from Gayatri's code:
//     - optimise (same in CPU version) uv_grid_counter, uv_grid_real, uv_grid_imag to become member variables.
Images CPacerImagerHip::gridding_imaging(Visibilities& xcorr,
                                     double delta_u, double delta_v,
                                     int    n_pixels,
                                     double min_uv /*=-1000*/,    // minimum UV 
                                     const char* weighting)
{
  std::cout << "Running 'gridding_imaging' on GPU.." << std::endl;

  int n_ant = xcorr.obsInfo.nAntennas;
  int xySize = n_ant * n_ant;
  int image_size {n_pixels * n_pixels}; 

  // update antenna flags before gridding which uses these flags or weights:
  UpdateAntennaFlags( n_ant );

   size_t n_images{xcorr.integration_intervals() * xcorr.nFrequencies};
   size_t buffer_size {image_size * n_images};
   if(!grids_counters) grids_counters.allocate(image_size * xcorr.nFrequencies, true);
   if(!grids) grids.allocate(buffer_size, true);
   MemoryBuffer<std::complex<float>> images_buffer(buffer_size, true);
  
  

   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   if(!u_gpu) {
      u_gpu.allocate(xySize, true);
      v_gpu.allocate(xySize, true);
   }
   gpuMemcpy(u_gpu.data(),  u_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);
   gpuMemcpy(v_gpu.data(),  v_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);

   gridding_gpu(xcorr, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, frequencies_gpu,
      delta_u, delta_v, n_pixels, min_uv, grids_counters, grids);

   gpuEvent_t start, stop;
   gpuEventCreate(&start);
   gpuEventCreate(&stop);
   float elapsed;
       

    if( !m_FFTPlan ){
       int n[2]; 
       n[0] = n_pixels; 
       n[1] = n_pixels;
       gpuEventRecord(start);
       gpufftPlanMany((gpufftHandle*)(&m_FFTPlan), 2, n, NULL, 1, image_size, NULL, 1, image_size, GPUFFT_C2C, n_images );
       gpuEventRecord(stop);
       gpuEventSynchronize(stop);
       gpuEventElapsedTime(&elapsed, start, stop);
       std::cout << "gpufftPlanMany took " << elapsed << "ms" << std::endl;
    }
     gpuEventRecord(start);
     gpufftExecC2C(m_FFTPlan, (gpufftComplex*) grids.data(), (gpufftComplex*) images_buffer.data(), GPUFFT_BACKWARD);
     gpuEventRecord(stop);
     gpuEventSynchronize(stop);
     gpuEventElapsedTime(&elapsed, start, stop);
     std::cout << "gpufftExecC2C took " << elapsed << "ms" << std::endl;
     if(!fnorm) fnorm.allocate(xcorr.nFrequencies, true);
     vector_sum_gpu(grids_counters.data(), image_size, xcorr.nFrequencies, fnorm);
     fft_shift_and_norm_gpu( (gpufftComplex*) images_buffer.data(), n_pixels, n_pixels, n_images, fnorm );
     Images imgs {std::move(images_buffer), xcorr.obsInfo, xcorr.nIntegrationSteps, xcorr.nAveragedChannels, static_cast<unsigned int>(n_pixels)};
      gpuEventDestroy(start);
     gpuEventDestroy(stop);
     if(CImagerParameters::averageImages){
          return image_averaging_gpu(imgs);
     }else{
          return imgs;
     }
}


void CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, MemoryBuffer<float>& w_cpu, MemoryBuffer<double>& frequencies){
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   if(!w_gpu) w_gpu.allocate(xySize, true);
   gpuMemcpy(w_gpu.data(), w_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);
   apply_geometric_corrections_gpu(xcorr, w_gpu.data(), frequencies_gpu);
}


void CPacerImagerHip::ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   apply_cable_lengths_corrections_gpu(xcorr, cable_lengths, frequencies_gpu);
}
