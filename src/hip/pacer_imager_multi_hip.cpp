#include "pacer_imager_multi_hip.h"
#include "pacer_imager_hip_defines.h"

// #include "gridding_imaging_cuda.h"
#include "gridding_multi_image_cuda.h"


#define __GPU__
#include <gpu_macros.hpp>
#include <gpu_fft.hpp>


// In order to use high resolution clock 
#include <ctime>
#include <ratio>
#include <chrono>

#include "../pacer_common.h"
#include <array_config_common.h>
#include <mystring.h>
#include <mydate.h>

CPacerImagerMultiFreqHip::CPacerImagerMultiFreqHip()
: CPacerImagerHip(), m_nStreams(15), m_bGridBlocks(false), m_GpuStreams(NULL), m_nCrossCorrBlockSize(0), m_nAutoCorrBlockSize(0),
  m_mapping_array_gpu(NULL), m_mapping_array_cpu(NULL), m_inputs_gpu(NULL)
{
   m_in_cross_correlations_gpu = NULL;
   m_in_auto_correlations_gpu  = NULL;
}

void CPacerImagerMultiFreqHip::SetNBlocks( int N )
{ 
   m_N = N;
   m_UVCounterSum.assign( m_N, 1.00 ); // set normalisation factors for all channels to 1 - later updated when UV counters are calculated   
   PRINTF_DEBUG("DEBUG : number of blocks(layers/channels) set to %d -> normalisation factors initialised to 1.00, m_UVCounterSum.size() = %d\n",N,int(m_UVCounterSum.size()));
}

void CPacerImagerMultiFreqHip::SetGridBlocks( bool bGridBlocks )
{ 
   m_bGridBlocks = bGridBlocks; 
   PRINTF_DEBUG("DEBUG : number of CUDA blocks set to %d\n",m_bGridBlocks);
}

void CPacerImagerMultiFreqHip::InitStreams()
{
   CleanStreams();

   m_GpuStreams = new gpuStream_t[m_nStreams];
   for(int i = 0; i < m_nStreams; i++)
   {
      gpuStreamCreate(&m_GpuStreams[i]); 
   }
}

void CPacerImagerMultiFreqHip::CleanStreams()
{
   if( m_GpuStreams ){
      for(int i=0;i<m_nStreams;i++)
      {
         gpuStreamDestroy( m_GpuStreams[i] );
      }
      delete [] m_GpuStreams;
   }
}

CPacerImagerMultiFreqHip::~CPacerImagerMultiFreqHip()
{
   CleanGPUMemory();

   CleanStreams();
}

bool CPacerImagerMultiFreqHip::CheckMemoryBuffers()
{
   if( !m_in_cross_correlations_gpu ){
       printf("ERROR : m_in_cross_correlations_gpu buffer not allocated\n");
       return false;
   }
   if( !m_in_auto_correlations_gpu ){
       printf("ERROR : m_in_auto_correlations_gpu buffer not allocated\n");
       return false;
   }

   if( !u_gpu )
   {
      printf("ERROR : u_gpu buffer not allocated\n");
      return false;
   }
 
   if( !v_gpu )
   {
      printf("ERROR : v_gpu buffer not allocated\n");
      return false;
   }
 
   if( !uv_grid_counter_gpu )
   {
      printf("ERROR : uv_grid_counter_gpu buffer not allocated\n");
      return false;
   }

   if( !uv_grid_counter_cpu )
   {
      printf("ERROR : uv_grid_counter_cpu buffer not allocated\n");
      return false;
   }

   if( !m_in_buffer_gpu )
   {
       printf("ERROR : m_in_buffer_gpu buffer not allocated\n");
       return false;
   }

   if( !m_out_buffer_gpu )
   {
      printf("ERROR : m_out_buffer_gpu buffer not allocated\n");
      return false;
   }

   if( !m_mapping_array_gpu ){
      printf("ERROR : m_mapping_array_gpu not allocated\n");
      return false;
   }
 
   if( !m_mapping_array_cpu ){
      printf("ERROR : m_mapping_array_cpu not allocated\n");
      return false;
   }

   if( !m_inputs_gpu ){
      printf("ERROR : m_inputs_gpu not allocated\n");
      return false;
   }
 
   return true;
}

bool CPacerImagerMultiFreqHip::ZeroMemoryBuffers()
{
   if( !CheckMemoryBuffers() ){
      return false;
   }

   // clean memory buffers required for imaging:
   if( !m_ImagerParameters.m_bConstantUVW )
   {
      PRINTF_DEBUG("ZEROing memory for UVW coordinates, XYsize = %d, m_N = %d\n",m_AllocatedXYSize,m_N);
      // if UVW is calculated for very new timestamp it has to be zeroed too:
      (gpuMemset((float*)u_gpu, 0, m_AllocatedXYSize*sizeof(float)*1));
      (gpuMemset((float*)v_gpu, 0, m_AllocatedXYSize*sizeof(float)*1));

      PRINTF_DEBUG("ZEROing memory for UV counter coordinates, XYsize = %d, m_N = %d\n",m_AllocatedImageSize,m_N);
      (gpuMemset((float*)uv_grid_counter_gpu, 0, m_AllocatedImageSize*sizeof(float)*m_N));
   }

   // UV grid used in gridding has to be always cleaned up :
   (gpuMemset((gpufftComplex*)m_in_buffer_gpu, 0, m_AllocatedImageSize*sizeof(gpufftComplex)*m_N));   

   return true;
}

void CPacerImagerMultiFreqHip::AllocInVisBuffer( int size_cross_corr , int size_auto_corr )
{
   if( !m_in_cross_correlations_gpu )
   {
      (gpuMalloc((void**) &m_in_cross_correlations_gpu, size_cross_corr*sizeof(gpufftComplex)));
      (gpuMemset((gpufftComplex*)m_in_cross_correlations_gpu, 0, size_cross_corr*sizeof(gpufftComplex)));
      printf("CPacerImagerMultiFreqHip::AllocInVisBuffer : allocated and initilised to zero %d bytes for cross correlations\n",int(size_cross_corr*sizeof(gpufftComplex)));

      m_nCrossCorrBlockSize = size_cross_corr;
      printf("DEBUG : allocated buffer for %d cross-correlations\n",m_nCrossCorrBlockSize);
   }   

   if( !m_in_auto_correlations_gpu )
   {
      (gpuMalloc((void**) &m_in_auto_correlations_gpu, size_auto_corr*sizeof(float)));
      (gpuMemset((float*)m_in_auto_correlations_gpu, 0, size_auto_corr*sizeof(float)));
      printf("CPacerImagerMultiFreqHip::AllocInVisBuffer : allocated and initilised to zero %d bytes for auto correlations\n",int(size_auto_corr*sizeof(gpufftComplex)));

      m_nAutoCorrBlockSize = size_auto_corr;
      printf("DEBUG : allocated buffer for %d auto-correlations\n",m_nAutoCorrBlockSize);
   }  
}

// Clean GPU Memory 
void CPacerImagerMultiFreqHip::CleanGPUMemory()
{
   if( m_in_cross_correlations_gpu ){
      (gpuFree( m_in_cross_correlations_gpu ));
      m_in_cross_correlations_gpu = NULL;
   }

   if( m_in_auto_correlations_gpu ){
      (gpuFree( m_in_auto_correlations_gpu ));
      m_in_auto_correlations_gpu = NULL;
   }

   if( m_mapping_array_gpu ){
      (gpuFree( m_mapping_array_gpu ));
      m_mapping_array_gpu = NULL;
   }

   if( m_mapping_array_cpu ){
      delete [] m_mapping_array_cpu;
      m_mapping_array_cpu = NULL;
   }

   if( m_inputs_gpu ){
      (gpuFree( m_inputs_gpu ));
      m_inputs_gpu = NULL;
   }

   // call standard clean from the parent class :
   CPacerImagerHip::CleanGPUMemory();
}


// Not all arrays are allocated as in CPacerImagerHip
void CPacerImagerMultiFreqHip::AllocGPUMemory( int corr_size, int image_size  ) 
{
   m_AllocatedXYSize = corr_size;
   m_AllocatedImageSize = image_size;

   if( !u_gpu )
   {
//      if( m_ImagerParameters.m_bConstantUVW ){
      // regardless if constant UVW or not -> no need to recalculate for all time blocks
      // WARNING : this is because UVW are stored in meters and are later divided by wavelength
      // -> this means for m_bConstantUVW we do not need separate UV array for different frequency channes (wavelengths)
      //    we just keep them in meters and then convert to meters/wavelenghts !
      (gpuMalloc((void**)&u_gpu,corr_size*sizeof(float)*1));
      (gpuMemset((float*)u_gpu, 0, corr_size*sizeof(float)*1));
//      }else{
//         (gpuMalloc((void**)&u_gpu, corr_size*sizeof(float)*m_N));
//         (gpuMemset((float*)u_gpu, 0, corr_size*sizeof(float)*m_N));
//      }
   }
   if( !v_gpu )
   {
//       if( m_ImagerParameters.m_bConstantUVW ){
      // regardless if constant UVW or not -> no need to recalculate for all time blocks
      // WARNING : this is because UVW are stored in meters and are later divided by wavelength
      // -> this means for m_bConstantUVW we do not need separate UV array for different frequency channes (wavelengths)
      //    we just keep them in meters and then convert to meters/wavelenghts !
      (gpuMalloc((void**)&v_gpu, corr_size*sizeof(float)*1));
      (gpuMemset((float*)v_gpu, 0, corr_size*sizeof(float)*1));
//      }else{
//         (gpuMalloc((void**)&v_gpu, corr_size*sizeof(float)*m_N));
//         (gpuMemset((float*)v_gpu, 0, corr_size*sizeof(float)*m_N));
//      }
   }
   if( !uv_grid_counter_gpu )
   {
      //if( m_ImagerParameters.m_bConstantUVW ){         
      //   (gpuMalloc((void**)&uv_grid_counter_gpu, image_size*sizeof(float)*N));
      //   (gpuMemset((float*)uv_grid_counter_gpu, 0, image_size*sizeof(float)*N));
      //}else{
      // WARNING : counter depends on wavelength as it operates in UV in units of wavelength (after division)
      (gpuMalloc((void**)&uv_grid_counter_gpu, image_size*sizeof(float)*m_N));
      (gpuMemset((float*)uv_grid_counter_gpu, 0, image_size*sizeof(float)*m_N));
      //}
   }
   
   if( !uv_grid_counter_cpu )
   {
      // WARNING : counter depends on wavelength as it operates in UV in units of wavelength (after division)
      // TODO : may just be a temporary code of for testing only 
      uv_grid_counter_cpu = new float[image_size*m_N];
      memset( uv_grid_counter_cpu, '\0', image_size*sizeof(float)*m_N);
   }
   
   // input and output to cuFFT :
   if( !m_in_buffer_gpu )
   {
      // m_in_buffer 
      (gpuMalloc((void**) &m_in_buffer_gpu, image_size*sizeof(gpufftComplex)*m_N));
      (gpuMemset((gpufftComplex*)m_in_buffer_gpu, 0, image_size*sizeof(gpufftComplex)*m_N));
   }   
   if( !m_out_buffer_gpu )
   {      
      (gpuMalloc((void**) &m_out_buffer_gpu, image_size*sizeof(gpufftComplex)*m_N)); 
      (gpuMemset((gpufftComplex*)m_out_buffer_gpu, 0, image_size*sizeof(gpufftComplex)*m_N));
   }
}

bool CPacerImagerMultiFreqHip::InitialiseAntennaMapping( std::vector< std::vector<int> >& antenna_mapping_array, std::vector<InputMapping>& inputs_cpu )
{
   // Memory for GPU input variables: 
   int n_inputs = antenna_mapping_array.size();
   printf("DEBUG : double check of inputs %d vs. %d -> OK = %d\n",n_inputs,int(inputs_cpu.size()),(n_inputs==inputs_cpu.size()));
   
   if( n_inputs != (2*m_nAntennas) ){
       printf("ERROR : number of inputs in the mapping array = %d != %d from the input data (for %d antennas)\n",n_inputs,(2*m_nAntennas),m_nAntennas);
       return false;
   } 

   int size = n_inputs*n_inputs;
   if( !m_mapping_array_gpu ){
      // just single - should also be some fast memory 
      // TODO : use texture memory ???
      (gpuMalloc((void**)&m_mapping_array_gpu,size*sizeof(int)*1));
      (gpuMemset((int*)m_mapping_array_gpu, 0, size*sizeof(int)*1));
   }
   if( !m_mapping_array_cpu ){
       m_mapping_array_cpu = new int[size];
       memset(m_mapping_array_cpu,'\0',size*sizeof(int));
   }
   if( !m_inputs_gpu ){
     (gpuMalloc((void**)&m_inputs_gpu,inputs_cpu.size()*sizeof(InputMapping)));
     // copy inputs to GPU memory :
     (gpuMemcpy((InputMapping*)m_inputs_gpu, (InputMapping*)(&(inputs_cpu[0])), sizeof(InputMapping)*inputs_cpu.size(), gpuMemcpyHostToDevice));
   }

   if( !CheckMemoryBuffers() ){
      printf("ERROR : at least one of the memory buffers has not been allocated !!!\n");
      return false;
   }

   printf("DEBUG : %d x %d vs. input array dimensions = %d x %d\n",n_inputs,n_inputs,int(antenna_mapping_array.size()),int(antenna_mapping_array[0].size()));
   if( antenna_mapping_array.size() != n_inputs || antenna_mapping_array[0].size()!=n_inputs ){
      printf("ERROR : in code, dimensions of the mapping arrays are not compatible !!!\n");
      return false;
   }
 
   // just copy one to one :
   for(int y=0;y<n_inputs;y++){
      for(int x=0;x<n_inputs;x++){
         int pos = y*n_inputs + x;
         m_mapping_array_cpu[pos] = antenna_mapping_array[y][x];         

         // TODO :
         // int pos2 = x*n_inputs+y;
         // m_mapping_array_cpu[pos2] = antenna_mapping_array[y][x];
      }
   } 
   // TODO : should I make this array m_mapping_array_cpu symetric -> fill
   // values -1 with the value from the other side of the diagonal ???
   // to avoid if-s in the gridding kernel ???
   // TODO : later optimisation once I am 100% code works ok !!!

/*   int index=0;
   for(int pos=0;pos<size;pos++){
     m_mapping_array_cpu[pos] = -1;
   }
   for(int y=0;y<n_inputs;y++){
      for(int x=0;x<n_inputs;x++){
         if( y < x ){
            int pos = y*n_inputs + x;
            m_mapping_array_cpu[pos] = antenna_mapping_array[y][x];
// this won't be OK because we only have 1 half of the matrix !!!
            if( index != antenna_mapping_array[y][x] ){
                printf("ERROR in CPacerImagerMultiFreqHip::InitialiseAntennaMapping code pos = %d != %d index\n",pos,index);
            }
//            printf("%d ",m_mapping_array_cpu[pos]);
            index++;
         }
      }
//      printf("\n");
   }*/

   // copy mapping matrix from CPU to GPU : 
   (gpuMemcpy((int*)m_mapping_array_gpu, (int*)m_mapping_array_cpu, sizeof(int)*size, gpuMemcpyHostToDevice));


   if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_ALL ){
      SaveDebugFiles( antenna_mapping_array, inputs_cpu, size );
   }


   return true;   
}

bool CPacerImagerMultiFreqHip::SaveGriddedVisibilities()
{
   bool ret = true;
   CBgFits tmp_fits( m_ImagerParameters.m_ImageSize, m_ImagerParameters.m_ImageSize );
   float* tmp_data = tmp_fits.get_data();

   gpufftComplex* tmp_buffer = new gpufftComplex[m_AllocatedImageSize*m_N];
   (gpuMemcpy(tmp_buffer, m_in_buffer_gpu,sizeof(gpufftComplex)*m_AllocatedImageSize*m_N, gpuMemcpyDeviceToHost));   

   for(int freq_channel=0; freq_channel<m_N; freq_channel++)
   {
      gpufftComplex* ptr_image_b_cpu = tmp_buffer + (freq_channel*m_AllocatedImageSize); // pointer to memory where gridded visibilities will be inserted

      for(int pos=0;pos<m_AllocatedImageSize;pos++){
         tmp_data[pos] = ptr_image_b_cpu[pos].x;
      }
      char szOutFits[128];
      sprintf(szOutFits,"gpu_uvgrid_real_ch%03d.fits",freq_channel);
      if( tmp_fits.WriteFits(szOutFits) ){
          printf("ERROR : could not write UV grid control file %s (channel %d)\n",szOutFits,freq_channel);          
          ret = false;
      }else{
          printf("DEBUG : control file file %s (channel %d) saved OK\n",szOutFits,freq_channel);
      }
   }
   
   delete [] tmp_buffer;   

   return ret;
}

void CPacerImagerMultiFreqHip::SaveDebugFiles( std::vector< std::vector<int> >& antenna_mapping_array, std::vector<InputMapping>& inputs_cpu, int size )
{
   printf("DEEP_DEBUG : debugging antenna_mapping_array -> GPU memory\n");
   if( m_mapping_array_gpu ){
      int* tmp_buffer = new int[size];
      (gpuMemcpy((int*)tmp_buffer, (int*)m_mapping_array_gpu, sizeof(int)*size, gpuMemcpyDeviceToHost));
      CBgFits tmp_fits( antenna_mapping_array.size(), antenna_mapping_array[0].size() );
      float* tmp_data = tmp_fits.get_data();
      for(int i=0;i<size;i++){
         tmp_data[i] = tmp_buffer[i];
//       if( tmp_buffer[i] <= -1  ){
//          tmp_data[i] = 0.00;
//       }
      }
      tmp_fits.WriteFits("debug_mapping_gpu.fits");

      delete [] tmp_buffer;
   }

   // save inputs 
   // CPU :
   FILE* out_f_cpu = fopen("debug_inputs_cpu.txt","w");
   for(int i=0;i<inputs_cpu.size();i++){
      InputMapping& input = inputs_cpu[i];

      fprintf(out_f_cpu,"%d %d %s %c %d\n",input.input,input.antenna,input.szAntName.c_str(),input.pol,input.flag);
   }
   fclose(out_f_cpu);

   // GPU :
   InputMapping* tmp_inputs = new InputMapping[inputs_cpu.size()];
   (gpuMemcpy((InputMapping*)tmp_inputs, (InputMapping*)m_inputs_gpu, sizeof(InputMapping)*inputs_cpu.size(),gpuMemcpyDeviceToHost));   
   FILE* out_f_gpu = fopen("debug_inputs_gpu.txt","w");
   for(int i=0;i<inputs_cpu.size();i++){
      InputMapping& input = tmp_inputs[i];

      fprintf(out_f_gpu,"%d %d %s %c %d\n",input.input,input.antenna,"",input.pol,input.flag);
   }
   fclose(out_f_gpu);
//   delete tmp_inputs;

}

bool CPacerImagerMultiFreqHip::InitialiseUVW()
{
   if( !u_gpu || !v_gpu || !uv_grid_counter_gpu || !uv_grid_counter_cpu ){
       printf("ERROR in CPacerImagerMultiFreqHip::InitialiseUVW() - gpu memory for UV and counter has not been allocated !!!\n");
       return false;
   }

   int xySize = m_U.GetXSize()*m_U.GetYSize();
   (gpuMemcpy((float*)u_gpu, (float*)m_U.get_data(), sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
   (gpuMemcpy((float*)v_gpu, (float*)m_V.get_data(), sizeof(float)*xySize, gpuMemcpyHostToDevice)); 

   
   // TODO : calcualte wavelength for each frequency channel 
   // calculate_counter<<<nBlocks,NTHREADS>>>(xySize, u_gpu, v_gpu, wavelength, image_size, delta_u, delta_v, n_pixels, center_x, center_y, is_odd_x, is_odd_y, vis_real_gpu, vis_imag_gpu, uv_grid_counter_single_gpu, min_uv );

   double channel_bw_mhz = m_ImagerParameters.m_fBandwidthMHz/m_N;
   double start_freq_mhz = m_ImagerParameters.m_fCenterFrequencyMHz - m_ImagerParameters.m_fBandwidthMHz/2.00;
   double first_channel_center_freq_mhz = start_freq_mhz + channel_bw_mhz/2.00;

   double FoV_radians = m_ImagerParameters.m_ImageFOV_degrees*M_PI/180; 
   double delta_u = 1.00/(FoV_radians);
   double delta_v = 1.00/(FoV_radians);

   int is_odd_x = 0;
   int is_odd_y = 0;
   if( (m_ImagerParameters.m_ImageSize % 2) != 0 ){
      is_odd_x = 1;
      is_odd_y = 1;
   }
   int center_x = int(m_ImagerParameters.m_ImageSize/2);
   int center_y = int(m_ImagerParameters.m_ImageSize/2);
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

   // 
   m_UVCounterSum.clear();

   // temporary to calculate Sum in a standard way 
   // TODO : move calculation as parallel reduction (loop over warp memory in a counter_sum kernel)
   CBgFits counter_fits("",0,0);
   
   // This happens once in a lifetime of the program so does not have to be superoptimal (loop is ok here)
   PACER_PROFILER_START
   printf("DEBUG : not using blocks (loop), xySize = %d, imageSize = %d\n",m_AllocatedXYSize,m_AllocatedImageSize);
   for(int b=0; b<m_N; b++)
   {
      float* ptr_counter_gpu = uv_grid_counter_gpu + (b*m_AllocatedImageSize);
      double freqMHz = first_channel_center_freq_mhz + b*channel_bw_mhz;
      double wavelength = SPEED_OF_LIGHT/(freqMHz*1e6); // c/freq_in_Hz

      if( CPacerImager::m_ImagerDebugLevel >= IMAGER_ALL_MSG_LEVEL ){
         printf("DEBUG : channel %d : pointer to GPU counter = %ld , wavelength = %.6f [m] ( frequency  = %.6f [MHz] )\n",b,(long int)ptr_counter_gpu,wavelength,freqMHz);
      }

      // implement calculate_counter_novis - same as calculate_counter but without visibility parameter !!!
      calculate_counter_novis<<<nBlocks,NTHREADS>>>( m_AllocatedXYSize, u_gpu, v_gpu, wavelength, m_AllocatedImageSize, delta_u, delta_v, m_ImagerParameters.m_ImageSize, 
                                                     center_x, center_y, is_odd_x, is_odd_y, ptr_counter_gpu, m_ImagerParameters.m_fMinUV );
   }
   gpuDeviceSynchronize();
   PACER_PROFILER_END("calculate_counter_novis took")

   // 2nd loop to keep the 1st one purly GPU kernel calls 
   // TODO : temporary : counter is on CPU to be able to save and calculat Sum -> later paralllel reduction on GPU
   for(int b=0; b<m_N; b++)
   {
      float* ptr_counter_gpu = uv_grid_counter_gpu + (b*m_AllocatedImageSize);

      // temporary : counter is on CPU to be able to save and calculat Sum -> later paralllel reduction on GPU 
      float* ptr_counter_cpu = uv_grid_counter_cpu + (b*m_AllocatedImageSize);      
      (gpuMemcpy(ptr_counter_cpu, ptr_counter_gpu, sizeof(float)*m_AllocatedImageSize, gpuMemcpyDeviceToHost));
      counter_fits.SetData( m_ImagerParameters.m_ImageSize, m_ImagerParameters.m_ImageSize, ptr_counter_cpu );
      double counter_sum = counter_fits.Sum();
      m_UVCounterSum.push_back( counter_sum );

      if( CPacerImager::m_ImagerDebugLevel >= IMAGER_ALL_MSG_LEVEL ){
         printf("DEBUG (CPacerImagerMultiFreqHip::InitialiseUVW) : normalisation factor calculated for channel %d to be %.8f\n",b,counter_sum);
      }

      if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_ALL ){
         // for debugging save counter file :
         char szCounterFits[128];
         sprintf(szCounterFits,"counter_channel%03d.fits",b);
         if( counter_fits.WriteFits( szCounterFits ) ){
             printf("ERROR : could not save counter fits file %s (for channel %d)\n",szCounterFits,b);
         }else{
             printf("DEBUG : image %s (freq. channel %d) saved OK\n",szCounterFits,b);
         }
      }
   }

   PRINTF_DEBUG("CPacerImagerMultiFreqHip::InitialiseUVW : m_UVCounterSum.size() = %d , m_N = %d\n",int(m_UVCounterSum.size()),m_N);

   return true;
}

void CPacerImagerMultiFreqHip::CopyInputVisibilities( float* pBaselinesAutoBufferGPU, 
                                                      gpufftComplex* pBaselinesCrossBuffer )
{
   (gpuMemcpy((float*)m_in_auto_correlations_gpu, pBaselinesAutoBufferGPU, m_nAutoCorrBlockSize*sizeof(float), gpuMemcpyHostToDevice));    
   (gpuMemcpy((gpufftComplex*)m_in_cross_correlations_gpu, pBaselinesCrossBuffer, m_nCrossCorrBlockSize*sizeof(gpufftComplex), gpuMemcpyHostToDevice));
}

void CPacerImagerMultiFreqHip::gridding_imaging_multi_freq( 
                                             int n_fine_channels, double bw_mhz, // description on input buffers (output of the correlator)
                                             // CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, - use m_U, m_V, m_W
                                             const char* weighting /*=""*/,   // weighting : U for uniform (others not implemented)
                                             bool do_gridding      /*=true*/, // if gridding is performed
                                             bool do_dirty_image   /*=true*/  // if cu/hip FFT is called
                                            )
{
    PRINTF_DEBUG("DEBUG : CPacerImagerMultiFreqHip::gridding_imaging_multi_freq start (m_bGridBlocks = %d)\n",m_bGridBlocks);
    InitStreams(); // only if not already initialised
    
    if( !CheckMemoryBuffers() ){     
       return;
    }

    // zero UV grid memory buffer and UV when required (non constant UV -> MWA case )
    ZeroMemoryBuffers(); // zero UV grid etc

    double channel_bw_mhz = m_ImagerParameters.m_fBandwidthMHz/m_N;
    double start_freq_mhz = m_ImagerParameters.m_fCenterFrequencyMHz - m_ImagerParameters.m_fBandwidthMHz/2.00;
    double first_channel_center_freq_mhz = start_freq_mhz + channel_bw_mhz/2.00;
    int    n_pixels = m_ImagerParameters.m_ImageSize;
    int    image_size = (n_pixels*n_pixels);
//    double min_uv   = m_ImagerParameters.m_fMinUV;

    double FoV_radians = m_ImagerParameters.m_ImageFOV_degrees*M_PI/180; 
    double delta_u = 1.00/(FoV_radians);
    double delta_v = 1.00/(FoV_radians);

    int nBlocks = (m_AllocatedXYSize + NTHREADS -1)/NTHREADS;
    int cudaBlocks = 1;   

    high_resolution_clock::time_point gridding_imaging_start = high_resolution_clock::now();

    if( m_bGridBlocks ){       
       cudaBlocks = n_fine_channels;
       PRINTF_INFO("INFO : executing gridding_imaging_lfile_vis_blocks (%d, %d, 1)\n",nBlocks,cudaBlocks);

       PACER_PROFILER_START
       dim3 blocks( nBlocks, n_fine_channels, 1 );
       if( m_ImagerParameters.m_bConstantUVW ){
            float* ptr_auto_correlations = NULL; 
            if( m_bIncludeAutos ){
               // if auto-correlation pointer is set, gridding_imaging_lfile_vis will use in the ant1=ant2 case
               ptr_auto_correlations = m_in_auto_correlations_gpu;
               PRINTF_DEBUG("Including auto-correlations\n");
            }

           // gridding kernel using direct output of the correlator :
           gridding_imaging_lfile_vis_blocks<<<nBlocks,NTHREADS>>>( m_AllocatedXYSize, 
                                                                    first_channel_center_freq_mhz, channel_bw_mhz,
                                                                    0, 'X', 0, 'X',
                                                                    u_gpu, v_gpu, 
                                                                    delta_u, delta_v, 
                                                                    m_in_cross_correlations_gpu, m_nCrossCorrBlockSize, 
                                                                    ptr_auto_correlations, m_nAutoCorrBlockSize,
                                                                    (gpufftComplex*)m_in_buffer_gpu, image_size, n_pixels, 
                                                                    m_nAntennas, n_fine_channels, m_inputs_gpu, (2*m_nAntennas), m_mapping_array_gpu 
                                                                  );

       }else{
           printf("ERROR : BLOCKs version with non-constant UVW not implemented yet !!!\n");
       }
       PACER_PROFILER_END("griding with blocks took")
    }else{
       // 
       PRINTF_INFO("DEBUG : CPacerImagerMultiFreqHip::gridding_imaging_multi_freq not using blocks (loop)\n");
       PACER_PROFILER_START
       for(int freq_channel=0; freq_channel<n_fine_channels; freq_channel++)
       {
           gpufftComplex* ptr_image_b_gpu = (gpufftComplex*)m_in_buffer_gpu + (freq_channel*image_size); // pointer to memory where gridded visibilities will be inserted

           double freqMHz = first_channel_center_freq_mhz + freq_channel*channel_bw_mhz;
           double wavelength = SPEED_OF_LIGHT/(freqMHz*1e6); // c/freq_in_Hz

           if( m_ImagerParameters.m_bConstantUVW ){
               // gridding kernel using direct output of the correlator : 
               float* ptr_auto_correlations = NULL; 
               if( m_bIncludeAutos ){
                  // if auto-correlation pointer is set, gridding_imaging_lfile_vis will use in the ant1=ant2 case
                  ptr_auto_correlations = m_in_auto_correlations_gpu;
                  PRINTF_DEBUG("Including auto-correlations\n");
               }

               gridding_imaging_lfile_vis<<<nBlocks,NTHREADS, 0, m_GpuStreams[freq_channel % m_nStreams]>>>( m_AllocatedXYSize, freq_channel, wavelength,
                                                                                                             0, 'X', 0, 'X',
                                                                                                             u_gpu, v_gpu, 
                                                                                                             delta_u, delta_v, 
                                                                                                             m_in_cross_correlations_gpu, m_nCrossCorrBlockSize, 
                                                                                                             ptr_auto_correlations, m_nAutoCorrBlockSize,
                                                                                                             ptr_image_b_gpu, image_size, n_pixels, 
                                                                                                             m_nAntennas, n_fine_channels, m_inputs_gpu, (2*m_nAntennas), m_mapping_array_gpu 
                                                                                                           );
           }else{ 
              // float* ptr_counter_gpu = uv_grid_counter_gpu + (freq_channel*image_size);
              printf("ERROR : non-blocks version with non-constant UVW not implemented yet !!!\n");              
           }
       }
       gpuDeviceSynchronize();
       PACER_PROFILER_END("griding with streams took")
    }

    if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_ALL ){
       // save gridded visibilities :
       SaveGriddedVisibilities();
    }

    if( !m_FFTPlan ){
       int n[2]; 
       n[0] = n_pixels; 
       n[1] = n_pixels; 

       // START: gpufftPLanMany() 
       std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

       // DONE : make plan -> m_FFTPlan and do it once in the lifetime of the program !!!
       //        it takes ~200ms and is not required every time as it is always the same !!!
       gpufftPlanMany((gpufftHandle*)(&m_FFTPlan), 2, n, NULL, 1, image_size, NULL, 1, image_size, GPUFFT_C2C, n_fine_channels );    
       std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
       std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
       PRINTF_BENCHMARK("BENCHMARK : gpufftPlanMany executed and took %.6f seconds. PARAMETERS ( N_PIXELS , N_BLOCKS , N_STREAMS , N_CHANNELS ) = ( %d , %d , %d , %d ) \n",time_span.count(),n_pixels,cudaBlocks,m_nStreams,n_fine_channels);
    }

    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    gpufftExecC2C(((gpufftHandle)m_FFTPlan), (gpufftComplex*)m_in_buffer_gpu, (gpufftComplex*)m_out_buffer_gpu, GPUFFT_BACKWARD);
    gpuDeviceSynchronize();

    // END: gpufftPlanMany() 
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    PRINTF_BENCHMARK("BENCHMARK : gpufftExecC2C took: %.6f seconds. PARAMETERS ( N_PIXELS , N_BLOCKS , N_STREAMS , N_CHANNELS ) = ( %d , %d , %d , %d ) \n",time_span.count(),n_pixels,cudaBlocks,m_nStreams,n_fine_channels);
    
    high_resolution_clock::time_point gridding_imaging_end = high_resolution_clock::now();
    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(gridding_imaging_end-gridding_imaging_start);
    PRINTF_BENCHMARK("BENCHMARK : full gridding_imaging took %.3f ms\n",time_span.count()*1000.00);
}

int CPacerImagerMultiFreqHip::CopyImagesGpu2Cpu()
{
   if(!m_out_buffer_cpu)
   {
      m_out_buffer_cpu = (gpufftComplex*)malloc(sizeof(gpufftComplex)*m_AllocatedImageSize*m_N);
      if( !m_out_buffer_cpu )
      {
         printf("ERROR : while allocating Host (Output) memory size  ...\n");
         return -1;
      }
   }  

   (gpuMemcpy(m_out_buffer_cpu, m_out_buffer_gpu, sizeof(gpufftComplex)*m_AllocatedImageSize*m_N, gpuMemcpyDeviceToHost));

   return m_N;
}

bool CPacerImagerMultiFreqHip::SaveSkyImage( int ch, int timestamp, double unixtime, int image_size, const char* szOutDir, int iSaveImaginaryFITS )
{
   CBgFits real_image( image_size, image_size ), imag_image( image_size, image_size );

   if( GetShiftedImage(ch,real_image,imag_image) ){
      char szOutFitsReal[1048],szOutFitsImag[1048],szOutFitsBase[1024],szPrefix[64];
      sprintf(szPrefix,"ch%05d/gpumulti_dirty_image_ch%03d_",ch,ch);
      get_filename_base( unixtime, szOutFitsBase, szOutDir, szPrefix );

      sprintf(szOutFitsReal,"%s_real.fits",szOutFitsBase);
      printf("DEBUG : saving channel %d image to FITS file %s\n",ch,szOutFitsReal);
      // real_image.WriteFits( szOutFitsReal );
      ((CPacerImager*)(this))->SaveSkyImage( szOutFitsReal, &real_image, unixtime );
      

      if( iSaveImaginaryFITS > 0 && (timestamp % iSaveImaginaryFITS)==0 ){
         sprintf(szOutFitsImag,"%s_imag.fits",szOutFitsBase);
         printf("DEBUG : saving channel %d image to FITS file %s\n",ch,szOutFitsImag);
         // imag_image.WriteFits( szOutFitsImag );
         ((CPacerImager*)(this))->SaveSkyImage( szOutFitsImag, &imag_image, unixtime );
      }

      return true;
   }
   
   return false; 
}

bool CPacerImagerMultiFreqHip::GetShiftedImage( int freq_channel, CBgFits& real_image_out, CBgFits& imag_image_out )
{
   CBgFits tmp_image_real( real_image_out.GetXSize(), real_image_out.GetYSize() ), tmp_image_imag( real_image_out.GetXSize(), real_image_out.GetYSize() );  
   double fnorm = 1.00;
   if( freq_channel < m_UVCounterSum.size() ){
      fnorm = 1.00/(m_UVCounterSum[freq_channel]);
      printf("DEBUG : normalisation factor for channel = %d is %.8f\n",freq_channel,fnorm);
   }else{
      printf("ERROR : m_UVCounterSum not initialised to handle at least %d frequency channels, m_UVCounterSum.size() = %d\n",(freq_channel+1),int(m_UVCounterSum.size()));
      return false;
   }
   
   // copy data from gpufftComplex -> CBgFits(float) :
   gpufftComplex* ptr_channel_image = (gpufftComplex*)m_out_buffer_cpu + m_AllocatedImageSize*freq_channel;
   float* ptr_data_real = tmp_image_real.get_data();
   float* ptr_data_imag = tmp_image_imag.get_data();
   for(int i=0;i<m_AllocatedImageSize;i++){
      ptr_data_real[i] = ptr_channel_image[i].x*fnorm;
      ptr_data_imag[i] = ptr_channel_image[i].y*fnorm;
   }

   // FFT Shift :
   fft_shift( tmp_image_real, real_image_out );     
   fft_shift( tmp_image_imag, imag_image_out );     

   return true;
}
