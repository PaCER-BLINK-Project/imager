 #include "pacer_imager_hip.h"
#include "pacer_imager_hip_defines.h"

#include "gridding_imaging_cuda.h"

#define __GPU__
#include <gpu_macros.hpp>
#include <gpu_fft.hpp>

#include "../pacer_common.h"
#include "gpu_utils.h"
#include <mystring.h>
#include <mydate.h>

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

void CPacerImagerHip::AllocGPUMemory( int corr_size, int image_size ) 
{
   m_AllocatedXYSize = corr_size;
   m_AllocatedImageSize = image_size;
   int n_ant = sqrt(corr_size);
   printf("DEBUG :  CPacerImagerHip::AllocGPUMemory( %d , %d ) -> n_ant = %d\n",corr_size,image_size,n_ant);

   // Memory for GPU input variables: 
   if( !u_gpu )
   {
      (gpuMalloc((void**)&u_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)u_gpu, 0, corr_size*sizeof(float)));
   }
   if( !v_gpu )
   {
      (gpuMalloc((void**)&v_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)v_gpu, 0, corr_size*sizeof(float)));
   }
   
  
   
}


void CPacerImagerHip::CleanGPUMemory()
{   
   if( vis_gpu )
   {
      (gpuFree( vis_gpu)); 
      vis_gpu = NULL;
   }

   if( u_gpu )
   {
      (gpuFree( u_gpu)); 
      u_gpu = NULL;
   }

   if( v_gpu )
   {
      (gpuFree( v_gpu)); 
      v_gpu = NULL;
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
Images CPacerImagerHip::gridding_imaging( Visibilities& xcorr, 
                                     int time_step, 
                                     int fine_channel,
                                     CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                                     double delta_u, double delta_v,
                                     int    n_pixels,
                                     double min_uv /*=-1000*/,    // minimum UV 
                                     const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
                                     const char* szBaseOutFitsName /*=NULL*/, 
                                     bool do_gridding,
                                     bool do_dirty_image,
                                     const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
                                     const char* in_fits_file_uv_im,  /*=""*/ // gridded visibilities can be provided externally
                                     bool bSaveIntermediate /*=false*/ , 
                                     bool bSaveImaginary /*=true*/
                )
{
  std::cout << "Running 'gridding_imaging' on GPU.." << std::endl;


  PACER_PROFILER_START
  bool bStatisticsCalculated = false;
  if( m_bPrintImageStatistics ){ // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
     fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );

     // V : 
     fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );

     // W : 
     fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );

     bStatisticsCalculated = true;
  }

  // Input size: u, v and w 
  int n_ant = xcorr.obsInfo.nAntennas;
  int xySize = n_ant * n_ant;


  int image_size {n_pixels * n_pixels}; 

  // In order to include conjugates at (-u,-v) UV point in gridding
  u_min = -u_max;
  v_min = -v_max;


  float *u_cpu = fits_vis_u.get_data();
  float *v_cpu = fits_vis_v.get_data();

  // Allocate only if not-allocated 

  // TODO : warning GPU UV grid is not initialised to ZEROs :
  // For the UV matrices, for now
  AllocGPUMemory(xySize, image_size ); //  out_image_real.get_data() );

   size_t n_images{xcorr.integration_intervals() * xcorr.nFrequencies};
   size_t buffer_size {image_size * n_images};
   MemoryBuffer<float> grids_counters_buffer(buffer_size, false, true);
   MemoryBuffer<std::complex<float>> grids_buffer(buffer_size, false,  true);
   MemoryBuffer<std::complex<float>> images_buffer(buffer_size, false, true);
   gpuMemset(grids_counters_buffer.data(), 0, n_images * image_size*sizeof(float));
   gpuMemset(grids_buffer.data(), 0, n_images * image_size*sizeof(std::complex<float>));
   
   
   
 
  // Step 3: Copy contents from CPU to GPU [input variables]
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy()
  (gpuMemcpy((float*)u_gpu, (float*)u_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)v_gpu, (float*)v_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice));
  
  
  // update antenna flags before gridding which uses these flags or weights:
  UpdateAntennaFlags( n_ant );

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
         double wavelength_m = VEL_LIGHT / frequency_hz;
         int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         std::complex<float>* current_grid = grids_buffer.data() + time_step * xcorr.nFrequencies * image_size + fine_channel * image_size;
         float* current_counter = grids_counters_buffer.data() + time_step * xcorr.nFrequencies * image_size + fine_channel * image_size;

         gridding_imaging_cuda_xcorr<<<nBlocks,NTHREADS>>>( xySize, n_ant, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, wavelength_m, image_size, delta_u, delta_v, n_pixels, vis_local_gpu, current_counter, min_uv, (gpufftComplex*)current_grid); 
         PRINTF_DEBUG("\n GRIDDING CHECK: Step 4 Calls to kernel");
         
         // Gives the error in the kernel! 
         gpuGetLastError();
         gpuDeviceSynchronize();
      }
   }
   
   grids_counters_buffer.to_cpu();
   grids_buffer.to_cpu();
   CBgFits reference_grid, reference_grid_counter;
    reference_grid.ReadFits("/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/uv_grid_real_8192x8192.fits", 0, 1, 1 );
    reference_grid_counter.ReadFits("/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/uv_grid_counter_8192x8192.fits", 0, 1, 1 );

    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(grids_counters_buffer[i] != reference_grid_counter.getXY(i % n_pixels, i / n_pixels)){
            std::cerr << "Error!! Counters are not the same at position " << i << ": " << grids_counters_buffer[i] << " != " << reference_grid_counter.getXY(i % n_pixels, i / n_pixels) << std::endl;
            break;
        }
    }
    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(std::abs(grids_buffer[i].real() - reference_grid.getXY(i % n_pixels, i / n_pixels)) > 1e-4){
            std::cerr << "Error!! Grids are not the same at position " << i << ": " << grids_buffer[i].real() << " != " << reference_grid.getXY(i % n_pixels, i / n_pixels) << std::endl;
            exit(1);
        }
    }
   grids_counters_buffer.to_gpu();
   grids_buffer.to_gpu();
  // TODO: 2024-06-22 : DIVIDE m_in_buffer_gpu and uv_grid_real_gpu, uv_grid_imag_gpu by uv_grid_counter_gpu for uniform and other weightings to really work
  //            are uv_grid_imag_gpu uv_grid_real_gpu really needed ???


  // TEMPORARY - check if UV Grid values are already truncated straight after after the kernel call
  // (gpuMemcpy((float*)uv_grid_real_cpu, (float*)uv_grid_real_gpu, sizeof(float)*image_size, gpuMemcpyDeviceToHost)); 
  // printf("\nDEBUG : GPU gridding (4,0) = %.20f [just after kernel call gridding_imaging_cuda]\n",m_uv_grid_real->getXY(4,0));  


  // uv_grid_counter_xSize = width
  // uv_grid_counter_ySize = height
  // size = image_size: (width x height)

  // Checking Execution time for cuFFT
    if( !m_FFTPlan ){
       int n[2]; 
       n[0] = n_pixels; 
       n[1] = n_pixels; 

       // START: gpufftPLanMany() 
       std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

       // DONE : make plan -> m_FFTPlan and do it once in the lifetime of the program !!!
       //        it takes ~200ms and is not required every time as it is always the same !!!
       gpufftPlanMany((gpufftHandle*)(&m_FFTPlan), 2, n, NULL, 1, image_size, NULL, 1, image_size, GPUFFT_C2C, n_images );    
       std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
       PRINTF_BENCHMARK("BENCHMARK : gpufftPlanMany executed and took %.6f seconds.\n",time_span.count());
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    gpufftExecC2C(m_FFTPlan, (gpufftComplex*) grids_buffer.data(), (gpufftComplex*) images_buffer.data(), GPUFFT_FORWARD);
    gpuDeviceSynchronize();

  // TODO: CPU->GPU : calculate this sum on GPU, can it be done in the gridding kernel itself ???
  // TODO: fix this
   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         std::complex<float>* current_image = images_buffer.data() + time_step * xcorr.nFrequencies * image_size + fine_channel * image_size;
         float* current_counter = grids_counters_buffer.data() + time_step * xcorr.nFrequencies * image_size + fine_channel * image_size;

         double fnorm = 1.00/sum_gpu_atomicadd( current_counter, image_size );
         // double fnorm = 1.00/sum_gpu_parallel_reduce( uv_grid_counter_gpu, image_size );

         
         // TODO : for now keeping it as it was as there is no clear advantage of doing normalsation on GPU 
         // apply normalisation constant on GPU :
         // int nBlocksImage = (image_size + NTHREADS -1)/NTHREADS;
         // mult_by_const<<<nBlocksImage,NTHREADS>>>( (gpufftComplex*)m_out_buffer_gpu, image_size, fnorm );
         
         // FFT shift together with multiplication by fnorm (normalisation)
         // bool fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, float fnorm=1.00 );
         fft_shift_and_norm_gpu( (gpufftComplex*) current_image, n_pixels, n_pixels, fnorm );
      }
   }
    return {std::move(images_buffer), xcorr.obsInfo, xcorr.nIntegrationSteps, xcorr.nAveragedChannels, static_cast<unsigned int>(n_pixels)};
}


bool CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_w, MemoryBuffer<double>& frequencies){
   std::cout << "Applying geometric corrections on GPU.." << std::endl;
   
   xcorr.to_cpu();
   ::compare_xcorr_to_fits_file(xcorr, "/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/01_before_geo_corrections.fits");
   
   
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   if(!w_gpu){
      gpuMalloc((void**)&w_gpu, xySize*sizeof(float));
   }
   gpuMemcpy(w_gpu, fits_vis_w.get_data(), sizeof(float)*xySize,  gpuMemcpyHostToDevice);
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_geometric_corrections<<<nBlocks,NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, w_gpu, frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
         gpuDeviceSynchronize();
      }
   }
   xcorr.to_cpu();
   ::compare_xcorr_to_fits_file(xcorr, "/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/01_after_geo_corrections.fits");
   
   printf("DEBUG : after call of apply_geometric_corrections kernel\n");
   return true;
}


bool CPacerImagerHip::ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   std::cout << "Applying cable corrections on GPU.." << std::endl;
   xcorr.to_cpu();
   ::compare_xcorr_to_fits_file(xcorr, "/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/02_before_cable_corrections.fits");
    

   if(!xcorr.on_gpu()) xcorr.to_gpu();
   cable_lengths.to_gpu();

   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_cable_corrections<<<nBlocks, NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, cable_lengths.data(), frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
         gpuDeviceSynchronize();
      }
   }
   xcorr.to_cpu();
   ::compare_xcorr_to_fits_file(xcorr, "/software/projects/director2183/cdipietrantonio/test-data/mwa/1276619416/imager_stages/1s_ch000/02_after_cable_corrections.fits");
    xcorr.to_gpu();
   return true;
}
