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



CPacerImagerHip::CPacerImagerHip()
: CPacerImager(),
  u_gpu(NULL), v_gpu(NULL), w_gpu(NULL),
  m_AllocatedXYSize(0), m_AllocatedImageSize(0),
  m_FFTPlan(0), vis_gpu(NULL), cable_lengths_gpu(NULL), cable_lengths_cpu(NULL),
  antenna_flags_gpu(NULL), antenna_weights_gpu(NULL), antenna_flags_cpu(NULL), antenna_weights_cpu(NULL), m_out_data(NULL)
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
   
   // gridded visibilities :      
   if( !uv_grid_counter_gpu )
   {
      (gpuMalloc((void**)&uv_grid_counter_gpu, image_size*sizeof(float)));
      (gpuMemset((float*)uv_grid_counter_gpu, 0, image_size*sizeof(float)));
   }

   if( !uv_grid_counter_cpu )
   {
      uv_grid_counter_cpu = new float[image_size];
      memset(uv_grid_counter_cpu,'\0',image_size*sizeof(float));
   }
   
   // input and output to cuFFT :
   if( !m_in_buffer_gpu )
   {
      // m_in_buffer 
      (gpuMalloc((void**) &m_in_buffer_gpu, image_size*sizeof(gpufftComplex)));
      (gpuMemset((gpufftComplex*)m_in_buffer_gpu, 0, image_size*sizeof(gpufftComplex)));
   }   
   if( !m_out_buffer_gpu )
   {      
      (gpuMalloc((void**) &m_out_buffer_gpu, image_size*sizeof(gpufftComplex))); 
      (gpuMemset((gpufftComplex*)m_out_buffer_gpu, 0, image_size*sizeof(gpufftComplex)));
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

   if( uv_grid_counter_gpu )
   {
      (gpuFree( uv_grid_counter_gpu)); 
      uv_grid_counter_gpu = NULL;
   }

   if( uv_grid_counter_cpu ){
      delete [] uv_grid_counter_cpu;
      uv_grid_counter_cpu = NULL;
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
void CPacerImagerHip::gridding_imaging( Visibilities& xcorr, 
                                     int time_step, 
                                     int fine_channel,
                                     CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                                     double delta_u, double delta_v,
                                     double frequency_mhz,
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
  int u_xSize = fits_vis_u.GetXSize();
  int u_ySize = fits_vis_u.GetYSize();
  int xySize = u_xSize*u_ySize;

  int vis_real_xSize = n_ant; 
  int vis_real_ySize = n_ant; 


  int image_size {n_pixels * n_pixels}; 
  int vis_real_size = (vis_real_xSize*vis_real_ySize);
  int vis_imag_size = (vis_real_xSize*vis_real_ySize);

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n SIZE CHECK xySize = %d", xySize);
     printf("\n SIZE CHECK image_size = %d", image_size);
     printf("\n SIZE CHECK vis_real_size = %d",vis_real_size); 
     printf("\n SIZE CHECK vis_imag_size = %d", vis_imag_size); 
  }

  // In order to include conjugates at (-u,-v) UV point in gridding
  u_min = -u_max;
  v_min = -v_max;


//   int center_x = int(n_pixels/2);
//   int center_y = int(n_pixels/2);
//   int is_odd_x = 0 , is_odd_y = 0;
//   if( (n_pixels % 2) == 1 ){
//      is_odd_x = 1;
//   }
//   if( (n_pixels % 2) == 1 ){
//      is_odd_y = 1;
//   }

  // initialise values in the UV GRID to zeros: 
  // 2024-06-14 : this code is not need as it is only used to initialise GPU arrays with zeros later in the code (cudaMemcpy), but these arrays are already initialised with zeros in 
  //              AllocGPUMemory function using cudaMemset
  // m_uv_grid_real->SetZeroValue();
  // m_uv_grid_imag->SetZeroValue();
  // m_uv_grid_counter->SetZeroValue();

  // Setting the initial values of out_image_real/out_image_imag 
  // MS (2024-06-12) : not required as these are overwritten later in a loop where normalisation factor fnorm is applied:
  // out_image_real.SetValue( 0.00 );
  // out_image_imag.SetValue( 0.00 );

  // Step 1: Declare GPU(Device) and CPU(Host) Variables 
  // CPU input variables 
  float *u_cpu = fits_vis_u.get_data();
  float *v_cpu = fits_vis_v.get_data();
  float *w_cpu = fits_vis_w.get_data();

  // Allocate only if not-allocated 

  // TODO : warning GPU UV grid is not initialised to ZEROs :
  AllocGPUMemory(xySize, image_size ); //  out_image_real.get_data() );
 
  // Step 3: Copy contents from CPU to GPU [input variables]
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy()
  (gpuMemcpy((float*)u_gpu, (float*)u_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)v_gpu, (float*)v_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice));
  
  // TODO : COPY xcorr strucuttre here:
  VISIBILITY_TYPE* vis_local_gpu = NULL;

  
  // update antenna flags before gridding which uses these flags or weights:
  UpdateAntennaFlags( n_ant );

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
         double wavelength_m = VEL_LIGHT / frequency_hz;
      int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
      gridding_imaging_cuda_xcorr<<<nBlocks,NTHREADS>>>( xySize, n_ant, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, wavelength_m, image_size, delta_u, delta_v, n_pixels, vis_local_gpu, uv_grid_counter_gpu, min_uv, (gpufftComplex*)m_in_buffer_gpu ); 
      PRINTF_DEBUG("\n GRIDDING CHECK: Step 4 Calls to kernel");
      
      // Gives the error in the kernel! 
      gpuGetLastError();
      gpuDeviceSynchronize();
      }
   }
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
     gpufftPlan2d((gpufftHandle*)(&m_FFTPlan), n_pixels, n_pixels, GPUFFT_C2C);
     PRINTF_INFO("INFO : gpufftPlan2d created\n");
  }
  gpufftExecC2C(((gpufftHandle)m_FFTPlan), (gpufftComplex*)m_in_buffer_gpu, (gpufftComplex*)m_out_buffer_gpu, GPUFFT_FORWARD);


  // Step 5: Copy contents from GPU variables to CPU variables
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy() 
  clock_t start_time5 = clock();

  // TODO: CPU->GPU : calculate this sum on GPU, can it be done in the gridding kernel itself ???
  // double fnorm = 1.00/m_uv_grid_counter->Sum();
  double fnorm = 1.00/sum_gpu_atomicadd( uv_grid_counter_gpu, image_size );
  // double fnorm = 1.00/sum_gpu_parallel_reduce( uv_grid_counter_gpu, image_size );

  
  // TODO : for now keeping it as it was as there is no clear advantage of doing normalsation on GPU 
  // apply normalisation constant on GPU :
  // int nBlocksImage = (image_size + NTHREADS -1)/NTHREADS;
  // mult_by_const<<<nBlocksImage,NTHREADS>>>( (gpufftComplex*)m_out_buffer_gpu, image_size, fnorm );
  
  // FFT shift together with multiplication by fnorm (normalisation)
  // bool fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, float fnorm=1.00 );
  fft_shift_and_norm_gpu( (gpufftComplex*)m_out_buffer_gpu, n_pixels, n_pixels, fnorm );
  
  // End of gpuMemcpy() GPU to CPU 
  clock_t end_time5 = clock();
  double duration_sec5 = double(end_time5-start_time5)/CLOCKS_PER_SEC;
  double duration_ms5 = duration_sec5*1000;
  printf("\n ** CLOCK gpuMemcpy() GPU to CPU took : %.6f [seconds], %.3f [ms]\n",duration_sec5,duration_ms5);
  PRINTF_DEBUG("\n GRIDDING CHECK: Step 5 GPU to CPU copied"); 

}


bool CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_w){
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
      }
   }
   printf("DEBUG : after call of apply_geometric_corrections kernel\n");
   return true;
}


bool CPacerImagerHip::ApplyCableCorrections( Visibilities& xcorr){
   if(!xcorr.on_gpu()) xcorr.to_gpu();

   if(!cable_lengths_cpu){
      cable_lengths_cpu = new float[xcorr.obsInfo.nAntennas];
      for(int ant=0; ant< xcorr.obsInfo.nAntennas; ant++){
         InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant]; 
         cable_lengths_cpu[ant] = ant1_info.cableLenDelta;
      }
   }
   if(!cable_lengths_gpu){
      gpuMalloc(&cable_lengths_gpu, xcorr.obsInfo.nAntennas*sizeof(float));
      gpuMemcpy(cable_lengths_gpu, cable_lengths_cpu, sizeof(float)*xcorr.obsInfo.nAntennas, gpuMemcpyHostToDevice);      
   }
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_cable_corrections<<<nBlocks, NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, cable_lengths_gpu, frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
      }
   }
   return true;
}
