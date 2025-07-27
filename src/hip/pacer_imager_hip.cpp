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

void memdump(char *ptr, size_t nbytes, std::string filename);

namespace {
       void compare_xcorr_to_fits_file(Visibilities& xcorr, std::string filename){
        auto vis2 = Visibilities::from_fits_file(filename, xcorr.obsInfo);
        size_t fine_channel {0}, int_time {0};
        size_t n_nans {0};
        size_t total {0};
        for(size_t a1 {0}; a1 < xcorr.obsInfo.nAntennas; a1++){
            for(size_t a2 {0}; a2 < a1; a2++){
                std::complex<float> *p1 = xcorr.at(int_time, fine_channel, a1, a2);
                std::complex<float> *p2 = vis2.at(int_time, fine_channel, a1, a2);
                for(size_t p {0}; p < 4; p++){
                    total++;
                    if(isnan(p1->real()) && isnan(p2->real()) && isnan(p2->imag()) && isnan(p1->imag())){
                        n_nans++;
                        continue;
                    }
                    if(*p1 != *p2){
                        std::cerr << "xcorr differs from " << filename << "!!!!" << std::endl;
                        std::cerr << "[a1 = " << a1 << ", a2 = " << a2 << "] p1 = " << *p1 << ", p2 = " << *p2 << std::endl;
                        exit(1);
                    }
                }
            }
        }
        std::cout << "OKK comparison with " << filename << std::endl;
        std::cout << "Percentage NaNs: " << (static_cast<double>(n_nans) / total * 100.0) << std::endl;
    }
}


CPacerImagerHip::CPacerImagerHip()
: CPacerImager(),
  m_FFTPlan(0), vis_gpu(NULL), cable_lengths_gpu(NULL), cable_lengths_cpu(NULL),
  antenna_flags_gpu(NULL), antenna_weights_gpu(NULL), antenna_flags_cpu(NULL), antenna_weights_cpu(NULL)
{

}

CPacerImagerHip::~CPacerImagerHip()
{
   CleanGPUMemory();
}



void CPacerImagerHip::CleanGPUMemory()
{   
   if( vis_gpu )
   {
      (gpuFree( vis_gpu)); 
      vis_gpu = NULL;
   }

   if( cable_lengths_gpu ){
      (gpuFree( cable_lengths_gpu )); 
      cable_lengths_gpu = NULL;
   }
   
   if( cable_lengths_cpu ){
      delete [] cable_lengths_cpu;
      cable_lengths_cpu = NULL;
   }
   

   
   // antenna flags and weights :
   if( antenna_flags_gpu ){
      gpuFree( antenna_flags_gpu );
      antenna_flags_gpu = NULL;
   }
   if( antenna_weights_gpu ){
      gpuFree( antenna_weights_gpu );
      antenna_weights_gpu = NULL;
   }
   if( antenna_flags_cpu ){
      delete [] antenna_flags_cpu;
      antenna_flags_cpu = NULL;
   }
   if( antenna_weights_cpu ){
      delete [] antenna_weights_cpu;
      antenna_weights_cpu = NULL;
   }
   


// TODO : why it is commented out - does it cause memory leak ???    
//   if( m_in_buffer_gpu )
//   {
//      (gpuFree( m_in_buffer_gpu)); 
//      m_in_buffer_gpu = NULL;
//   }
//   if( m_out_buffer_gpu )
//   {
//      (gpuFree( m_out_buffer_gpu)); 
//      m_out_buffer_gpu = NULL;
//   }
}


void CPacerImagerHip::UpdateAntennaFlags( int n_ant )
{
   if(!antenna_flags_cpu){
      antenna_flags_cpu = new int[n_ant];      
   }

   if(!antenna_weights_cpu){
      antenna_weights_cpu = new float[n_ant];
   }
   
   // antenna flags and weights:
   if( !antenna_flags_gpu )
   {
      (gpuMalloc((void**)&antenna_flags_gpu, n_ant*sizeof(int)));
      (gpuMemset((int*)antenna_flags_gpu, 0, n_ant*sizeof(int)));
   }

   if( !antenna_weights_gpu )
   {
      (gpuMalloc((void**)&antenna_weights_gpu, n_ant*sizeof(float)));
      (gpuMemset((float*)antenna_weights_gpu, 0, n_ant*sizeof(float)));
   }


   int n_flagged=0;   
   for(int ant=0;ant<n_ant;ant++){
      InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant]; 
      
      antenna_flags_cpu[ant] = ant1_info.flag;
      antenna_weights_cpu[ant] = 1.00;
      if( ant1_info.flag ){
         antenna_weights_cpu[ant] = 0.00;
         n_flagged++;
      }
   }
   
   PRINTF_INFO("INFO : CPacerImagerHip::UpdateAntennaFlags there are %d flagged antennas passed to GPU imager\n",n_flagged);

   // copy to device:
   (gpuMemcpy((int*)antenna_flags_gpu, (int*)antenna_flags_cpu, sizeof(int)*n_ant, gpuMemcpyHostToDevice));
   (gpuMemcpy((float*)antenna_weights_gpu, (float*)antenna_weights_cpu, sizeof(float)*n_ant, gpuMemcpyHostToDevice));

}

template <typename T>
inline void compare_buffers(MemoryBuffer<T>& a, MemoryBuffer<T>&b ){
   for(size_t i {0}; i < a.size(); i++){
      if(std::abs(a[i] - b[i]) >= 1e-4) {
         std::cout << "Elements differ at position " << i << ": a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
         throw std::exception{};
      }
   }
}

template <>
inline void compare_buffers(MemoryBuffer<std::complex<float>>& a, MemoryBuffer<std::complex<float>>&b ){
   for(size_t i {0}; i < a.size(); i++){
      if(std::abs(a[i].real() - b[i].real()) >= 1e-4 || std::abs(a[i].imag() - b[i].imag()) >= 1e-4) {
         std::cout << "Elements differ at position " << i << ": a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
         //throw std::exception{};
      }
   }
}



// TODO : 
//     - do more cleanup of this function as this is nearly "RAW" copy paste from Gayatri's code:
//     - optimise (same in CPU version) uv_grid_counter, uv_grid_real, uv_grid_imag to become member variables.
Images CPacerImagerHip::gridding_imaging(Visibilities& xcorr, 
                                     int time_step, 
                                     int fine_channel,
                                     CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                                     double delta_u, double delta_v,
                                     int    n_pixels,
                                     double min_uv /*=-1000*/,    // minimum UV 
                                     const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
                                     const char* szBaseOutFitsName /*=NULL*/
                )
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
   gpuMemcpy(u_gpu.data(),  fits_vis_u.get_data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);
   gpuMemcpy(v_gpu.data(),  fits_vis_v.get_data(), sizeof(float)*xySize, gpuMemcpyHostToDevice);

   gridding_gpu(xcorr, time_step, fine_channel, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, frequencies_gpu,
      delta_u, delta_v, n_pixels, min_uv, grids_counters, grids);

   // auto ref_grids_counters = MemoryBuffer<float>::from_dump("/scratch/director2183/cdipietrantonio/cpu_stages_dumps/grids_counters.bin");
   // auto ref_grids =  MemoryBuffer<std::complex<float>>::from_dump("/scratch/director2183/cdipietrantonio/cpu_stages_dumps/grids.bin");
   // grids_counters.to_cpu();
   // grids.to_cpu();
   // std::cout << "Comparing counters..." << std::endl;
   // compare_buffers(ref_grids_counters, grids_counters);
   // std::cout << "Comparing grids..." << std::endl;
   // compare_buffers(ref_grids, grids);
   // grids_counters.to_gpu();
   // grids.to_gpu();
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


void CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_w, MemoryBuffer<double>& frequencies){
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   if(!w_gpu) w_gpu.allocate(xySize, true);
   gpuMemcpy(w_gpu.data(), fits_vis_w.get_data(), sizeof(float)*xySize,  gpuMemcpyHostToDevice);
   apply_geometric_corrections_gpu(xcorr, w_gpu.data(), frequencies_gpu);
}


void CPacerImagerHip::ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   if(!frequencies_gpu) frequencies_gpu.allocate(xcorr.nFrequencies, true);
   gpuMemcpy(frequencies_gpu.data(), frequencies.data(), frequencies.size() * sizeof(double), gpuMemcpyHostToDevice);
   apply_cable_lengths_corrections_gpu(xcorr, cable_lengths, frequencies_gpu);
}
