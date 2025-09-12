#include "pacer_imager_hip.h"
#include "pacer_imager_hip_defines.h"
#include "gridding_gpu.h"
#include <gpu_macros.hpp>
#include <gpu_fft.hpp>

#include "../pacer_common.h"
#include "gpu_utils.h"
#include <mystring.h>
#include <mydate.h>
#include <exception>
#include "corrections_gpu.h"

#include "../utils.h"

CPacerImagerHip::CPacerImagerHip(const std::string metadata_file, int n_pixels, const std::vector<int>& flagged_antennas,
   bool average_images, Polarization pol_to_image, float oversampling_factor, double min_uv, const char *weighting) : CPacerImager(
      metadata_file, n_pixels, flagged_antennas, average_images, pol_to_image, oversampling_factor, min_uv, weighting) {}

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


void CPacerImagerHip::gridding(Visibilities& xcorr){
  std::cout << "Running 'gridding' on GPU.." << std::endl;

  int n_ant = xcorr.obsInfo.nAntennas;
  int xySize = n_ant * n_ant;
  int image_size {n_pixels * n_pixels};
   size_t n_images {n_gridded_channels * n_gridded_intervals};
   size_t buffer_size {image_size * n_images};

  // update antenna flags before gridding which uses these flags or weights:
  UpdateAntennaFlags( n_ant );
  if(!u_gpu) {
      u_gpu.allocate(xySize, true);
      v_gpu.allocate(xySize, true);
   }
   gpuMemcpy(u_gpu.data(),  u_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);
   gpuMemcpy(v_gpu.data(),  v_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);

   if(!grids_counters) {
      grids_counters.allocate(image_size * xcorr.nFrequencies, true);
      gpuMemset(grids_counters.data(), 0, grids_counters.size() * sizeof(float));
   }
   if(!grids) {
      grids.allocate(buffer_size, true);
      gpuMemset(grids.data(), 0, grids.size() * sizeof(std::complex<float>));
   }
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   

   gridding_gpu(xcorr, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, frequencies_gpu,
      delta_u, delta_v, n_pixels, min_uv, pol_to_image, grids_counters, grids);   
}



Images CPacerImagerHip::image(ObservationInfo& obs_info){
   if(!grids_counters || !grids) throw std::runtime_error {"Grids are not initialised and cannot be imaged."};
   int image_size {n_pixels * n_pixels};
   size_t n_images {n_gridded_channels * n_gridded_intervals};
   size_t buffer_size {image_size * n_images};

   MemoryBuffer<std::complex<float>> images_buffer(buffer_size, true);
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
    gpuEventDestroy(start);
     gpuEventDestroy(stop);
     
     if(!fnorm) fnorm.allocate(n_gridded_channels, true);
     vector_sum_gpu(grids_counters.data(), image_size, n_gridded_channels, fnorm);
     fft_shift_and_norm_gpu( (gpufftComplex*) images_buffer.data(), n_pixels, n_pixels, n_images, fnorm );
     // reset grids for the next round of imaging
     gpuMemset(grids_counters.data(), 0, grids_counters.size() * sizeof(float));
     gpuMemset(grids.data(), 0, grids.size() * sizeof(std::complex<float>));
   
     Images images {std::move(images_buffer), obs_info,static_cast<unsigned int>(n_gridded_intervals), 
        static_cast<unsigned int>(n_gridded_channels), static_cast<unsigned int>(n_pixels)};
     
    images.ra_deg = m_MetaData.raHrs*15.00;
    images.dec_deg = m_MetaData.decDegs;
    // MAX(u) , pixscale in degrees is later used in WCS FITS keywords CDELT1,2
    images.pixscale_ra = (1.00/(oversampling_factor*u_max))*(180.00/M_PI);
    // MAX(v) , pixscale in degrees is later used in WCS FITS keywords CDELT1,2
    images.pixscale_dec = (1.00/(oversampling_factor*v_max))*(180.00/M_PI);
     if(average_images){
          return image_averaging_gpu(images);
     }else{
          return images;
     }

}


void CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, MemoryBuffer<float>& w_cpu, MemoryBuffer<double>& frequencies){
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
    // TODO: improve the following
   if(!w_gpu) w_gpu.allocate(xySize, true);
   gpuMemcpy(w_gpu.data(), w_cpu.data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);
   apply_geometric_corrections_gpu(xcorr, w_gpu.data(), frequencies_gpu);
}


void CPacerImagerHip::ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
    if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
    gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice); 
    apply_cable_lengths_corrections_gpu(xcorr, cable_lengths, frequencies_gpu);
}
