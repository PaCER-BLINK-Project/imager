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

#include "corrections_gpu.h"

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
  u_gpu(NULL), v_gpu(NULL), w_gpu(NULL),
  m_AllocatedXYSize(0), m_AllocatedImageSize(0),
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


   if( w_gpu )
   {
      (gpuFree( w_gpu)); 
      w_gpu = NULL;
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
   MemoryBuffer<float> grids_counters_buffer(buffer_size, false, true);
   MemoryBuffer<std::complex<float>> grids_buffer(buffer_size, false,  true);
   MemoryBuffer<std::complex<float>> images_buffer(buffer_size, false, true);
  
  

   MemoryBuffer<double> frequencies {xcorr.nFrequencies, false, false};
   for(size_t fine_channel {0}; fine_channel < xcorr.nFrequencies; fine_channel++)
      frequencies[fine_channel] = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
   
   gridding_gpu(xcorr, time_step, fine_channel, fits_vis_u, fits_vis_v, antenna_flags_gpu, antenna_weights_gpu, frequencies,
      delta_u, delta_v, n_pixels, min_uv, grids_counters_buffer, grids_buffer);

    if( !m_FFTPlan ){
       int n[2]; 
       n[0] = n_pixels; 
       n[1] = n_pixels; 

       // START: gpufftPLanMany() 
       std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

       // DONE : make plan -> m_FFTPlan and do it once in the lifetime of the program !!!
       //        it takes ~200ms and is not required every time as it is always the same !!!
       gpufftPlanMany((gpufftHandle*)(&m_FFTPlan), 2, n, NULL, 1, image_size, NULL, 1, image_size, GPUFFT_C2C, n_images );
    }
    gpufftExecC2C(m_FFTPlan, (gpufftComplex*) grids_buffer.data(), (gpufftComplex*) images_buffer.data(), GPUFFT_FORWARD);
    MemoryBuffer<float> fnorm {n_images, false, true};
    sum_gpu_atomicadd(grids_counters_buffer.data(), image_size, n_images, fnorm);
    fft_shift_and_norm_gpu( (gpufftComplex*) images_buffer.data(), n_pixels, n_pixels, n_images, fnorm );

    return {std::move(images_buffer), xcorr.obsInfo, xcorr.nIntegrationSteps, xcorr.nAveragedChannels, static_cast<unsigned int>(n_pixels)};
}


bool CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_w, MemoryBuffer<double>& frequencies){
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   if(!w_gpu){
      gpuMalloc((void**)&w_gpu, xySize*sizeof(float));
   }
   gpuMemcpy(w_gpu, fits_vis_w.get_data(), sizeof(float)*xySize,  gpuMemcpyHostToDevice);
   apply_geometric_corrections_gpu(xcorr, w_gpu, frequencies);
   return true;
}


bool CPacerImagerHip::ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   apply_cable_lengths_corrections_gpu(xcorr, cable_lengths, frequencies);
   return true;
}
