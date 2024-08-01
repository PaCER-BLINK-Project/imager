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

/*void __cuda_check_error(cudaError_t err, const char *file, int line)
{
	if(err != cudaSuccess){
        fprintf(stderr, "CUDA error (%s:%d): %s\n", file, line, cudaGetErrorString(err));
        exit(1);
    }
}*/


CPacerImagerHip::CPacerImagerHip()
: CPacerImager(),
  u_gpu(NULL), v_gpu(NULL), w_gpu(NULL), vis_real_gpu(NULL), vis_imag_gpu(NULL), uv_grid_real_gpu(NULL), uv_grid_imag_gpu(NULL), uv_grid_counter_gpu(NULL), uv_grid_counter_cpu(NULL),
  m_in_buffer_gpu(NULL), m_out_buffer_gpu(NULL), m_AllocatedXYSize(0), m_AllocatedImageSize(0), m_out_buffer_cpu(NULL),
  m_FFTPlan(0), vis_gpu(NULL), cable_lengths_gpu(NULL), cable_lengths_cpu(NULL), test_data_real_gpu(NULL), test_data_imag_gpu(NULL),
  antenna_flags_gpu(NULL), antenna_weights_gpu(NULL), antenna_flags_cpu(NULL), antenna_weights_cpu(NULL), m_out_data(NULL)
{

}

CPacerImagerHip::~CPacerImagerHip()
{
   CleanGPUMemory();
}

void CPacerImagerHip::AllocGPUMemoryForXCorr( Visibilities* p_xcorr )
{
   if( p_xcorr ){
      if( p_xcorr->on_gpu() ){
         PRINTF_INFO("INFO : data already on GPU -> no need to create another buffer and copy\n");
      }else{
         // only allocate local memory buffer if xcorr structure keeps it on CPU 
         // otherwise (data is on GPU already) do not create another copy 
         if( !vis_gpu ){
            // int vis_size_bytes = sizeof(std::complex<VISIBILITY_TYPE>) * p_xcorr->size();
            int vis_size_bytes = sizeof(std::complex<VISIBILITY_TYPE>) * p_xcorr->matrix_size(); // 2024-04-07 : only allocate memory for a single correlation matrix in this version
            (gpuMalloc((void**)&vis_gpu, vis_size_bytes));
            (gpuMemset((VISIBILITY_TYPE*)vis_gpu, 0, vis_size_bytes));
         }
      }
   }else{
      printf("INFO : p_xcorr = NULL -> no need for vis_gpu allocation\n");
   }
}

void CPacerImagerHip::AllocGPUMemory( int corr_size, int image_size ) 
{
   m_AllocatedXYSize = corr_size;
   m_AllocatedImageSize = image_size;
   int n_ant = sqrt(corr_size);
   printf("DEBUG :  CPacerImagerHip::AllocGPUMemory( %d , %d ) -> n_ant = %d\n",corr_size,image_size,n_ant);

   // Memory for GPU input variables: 
   if( !vis_real_gpu )
   {
      gpuMalloc((void**)&vis_real_gpu, corr_size*sizeof(float));
      gpuMemset((float*)vis_real_gpu, 0, corr_size*sizeof(float));
   }
   if( !vis_imag_gpu )
   {
      (gpuMalloc((void**)&vis_imag_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)vis_imag_gpu, 0, corr_size*sizeof(float)));
   }
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
   if( !w_gpu )
   {
      (gpuMalloc((void**)&w_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)w_gpu, 0, corr_size*sizeof(float)));
   }

   // gridded visibilities :   
   if( !uv_grid_real_gpu )
   {
      // Memory for GPU output variables:  
      (gpuMalloc((void**)&uv_grid_real_gpu, image_size*sizeof(float)));
      (gpuMemset((float*)uv_grid_real_gpu, 0, image_size*sizeof(float)));      
   }
   
   if( !uv_grid_imag_gpu )
   {
      (gpuMalloc((void**)&uv_grid_imag_gpu, image_size*sizeof(float)));
      (gpuMemset((float*)uv_grid_imag_gpu, 0, image_size*sizeof(float)));
   }
   
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
   
   // test buffers :
   if( !test_data_real_gpu ){
      (gpuMalloc((void**)&test_data_real_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)test_data_real_gpu, 0, corr_size*sizeof(float)));
   }
   if( !test_data_imag_gpu ){
      (gpuMalloc((void**)&test_data_imag_gpu, corr_size*sizeof(float)));
      (gpuMemset((float*)test_data_imag_gpu, 0, corr_size*sizeof(float)));
   }
   
   if( !m_out_data ){
      m_out_data = (gpufftComplex*)malloc(sizeof(gpufftComplex) * image_size);
   }
   
}


void CPacerImagerHip::CleanGPUMemory()
{
   if( vis_real_gpu )
   {
      (gpuFree( vis_real_gpu)); 
      vis_real_gpu = NULL;
   }
  
   if( vis_imag_gpu )
   {
      (gpuFree( vis_imag_gpu)); 
      vis_imag_gpu = NULL;
   }
   
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

   if( uv_grid_real_gpu )
   {
      (gpuFree( uv_grid_real_gpu)); 
      uv_grid_real_gpu = NULL;
   }

   if( uv_grid_imag_gpu )
   {
      (gpuFree( uv_grid_imag_gpu)); 
      uv_grid_imag_gpu = NULL;
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
   
   // clean test buffers:
   if( test_data_real_gpu ){
      (gpuFree( test_data_real_gpu ));
      test_data_real_gpu = NULL;
   }
   if( test_data_imag_gpu ){
      (gpuFree( test_data_imag_gpu ));
      test_data_imag_gpu = NULL;
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
   
   if( m_out_data )
   {
      free( (gpufftComplex*)m_out_data);
      m_out_data = NULL;
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

void CPacerImagerHip::InitCableLengts( int n_ant )
{
   if( !cable_lengths_gpu ){
      // allocate cpu library :
      if( !cable_lengths_cpu ){
         cable_lengths_cpu = new float[n_ant];
      }

      // allocate gpu memory :
      (gpuMalloc((void**)&cable_lengths_gpu, n_ant*sizeof(float)));
      (gpuMemset((float*)cable_lengths_gpu, 0, n_ant*sizeof(float)));

      
      for(int ant=0;ant<n_ant;ant++){
         InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant]; 
         cable_lengths_cpu[ant] = ant1_info.cableLenDelta;
      }
      
      // copy to device:
      (gpuMemcpy((float*)cable_lengths_gpu, (float*)cable_lengths_cpu, sizeof(float)*n_ant, gpuMemcpyHostToDevice));      
   }
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

// Saving intermediate / test FITS files and printing statistics:
void CPacerImagerHip::SaveTestFitsFilesAndShowStat( int n_pixels, 
                                                    const char* weighting,
                                                    const char* szBaseOutFitsName, 
                                                    bool bSaveIntermediate, 
                                                    bool bSaveImaginary 
                                                  )
{
   if( CPacerImager::m_SaveFilesLevel > 0 ){
       // out_image_real and out_image_imag 
       // 2024-06-23 : these images are not saved in GPU version. In CPU version there were after 2D FFT and normalisation but before FFTShift. 
       //              However, in GPU version 2D FFT is followed by FFTShift+normalisation in one kernel (no DeviceToHost copy in between) -> cannot save this product
       // CBgFits out_image_real( m_uv_grid_real->GetXSize(), m_uv_grid_real->GetYSize() ), out_image_imag( m_uv_grid_real->GetXSize(), m_uv_grid_real->GetYSize() ); 
   
       // CPU Output variables
       // MS (2024-06-14) I leave these for know, but in fact they should only be used/executed if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG )
       //                 because otherwise we do not need to have all these data copied from GPU to CPU
       // TODO : add if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ) -> will require more testing
       int uv_grid_counter_xSize = m_uv_grid_counter->GetXSize();
       int uv_grid_counter_ySize = m_uv_grid_counter->GetYSize();
       int image_size = (uv_grid_counter_xSize*uv_grid_counter_ySize); 

       float *uv_grid_real_cpu = m_uv_grid_real->get_data();
       float *uv_grid_imag_cpu = m_uv_grid_imag->get_data();
       float *uv_grid_counter_cpu = m_uv_grid_counter->get_data();
 
       // 2024-06-22: TODO : move saving FITS files to a separate function:
       // Saving gridding() output files 
       if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG )
       {
          (gpuMemcpy((float*)uv_grid_counter_cpu, (float*)uv_grid_counter_gpu, sizeof(float)*image_size, gpuMemcpyDeviceToHost)); 
          (gpuMemcpy((float*)uv_grid_real_cpu, (float*)uv_grid_real_gpu, sizeof(float)*image_size, gpuMemcpyDeviceToHost)); 
          (gpuMemcpy((float*)uv_grid_imag_cpu, (float*)uv_grid_imag_gpu, sizeof(float)*image_size, gpuMemcpyDeviceToHost)); 
          PRINTF_DEBUG("\nDEBUG : GPU gridding (4,0) = %.20f [just after gpuMemcpy]\n",m_uv_grid_real->getXY(4,0));
  
          // WARNING : this is only required for debugging (saving intermediate test files) - hence moved inside this if 
          // TODO: Uniform weighting (Not implemented in GPU version). Divie UV grid real/imag by counter before FFT-2D
          if( strcmp(weighting, "U" ) == 0 )
          {
             m_uv_grid_real->Divide( *m_uv_grid_counter );
             m_uv_grid_imag->Divide( *m_uv_grid_counter );
          }      
  
          char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
          sprintf(uv_grid_re_name,"%s/uv_grid_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
          sprintf(uv_grid_im_name,"%s/uv_grid_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
          sprintf(uv_grid_counter_name,"%s/uv_grid_counter_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
    
          PRINTF_DEBUG("\nDEBUG : GPU gridding (4,0) = %.20f [just before saving]\n",m_uv_grid_real->getXY(4,0));
          if( m_uv_grid_real->WriteFits( uv_grid_re_name ) ){
             printf("ERROR : could not write output file %s\n",uv_grid_re_name);
          }else{
             PRINTF_INFO("INFO : saved file %s\n",uv_grid_re_name);
          }

          if( m_uv_grid_imag->WriteFits( uv_grid_im_name ) ){
             printf("ERROR : could not write output file %s\n",uv_grid_im_name);
          }else{
             PRINTF_INFO("INFO : saved file %s\n",uv_grid_im_name);
          }
  
          if( m_uv_grid_counter->WriteFits( uv_grid_counter_name ) ){
             printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
          }else{
             PRINTF_INFO("INFO : saved file %s\n",uv_grid_counter_name);
          }
       }

       // CPU Variable
       // TODO : this will only be done if saving intermediate / debug FITS files is required, otherwise stays on GPU :   
       (gpuMemcpy( (gpufftComplex*)m_out_data, m_out_buffer_gpu, sizeof(gpufftComplex)*image_size, gpuMemcpyDeviceToHost));

       // 2024-06-22: TODO : move saving FITS files to a separate function:

       // 2024-06-22 : fft shift is executed on GPU together with normalisation (multiplication by fnorm)
       // TODO : once fft_shift is implemneted on GPU this code may be move inside if( bSaveIntermediate ){ 
       //        then this operation will be performed only if really required for saving intermediate files 
       //        Once fft_shift is performed on GPU multipilcation by fnorm can also be done there. 
       //        Then final output will be copied from GPU to CPU only when required by CPacerImager::m_SaveFilesLevel parameter
       // 
       // float pointers to 1D Arrays 
       // float* out_data_real = out_image_real.get_data();
       // float* out_data_imag = out_image_imag.get_data();

       // DONE : next part to move to GPU
       // Assigning back 
       /*for(int i = 0; i < image_size; i++) 
       {
          out_data_real[i] = ((gpufftComplex*)m_out_data)[i].x*fnorm; // was *fnorm - now on GPU
          out_data_imag[i] = ((gpufftComplex*)m_out_data)[i].y*fnorm; // was *fnorm - now on GPU
       } */  

       // 2024-06-22 : temporary, TODO : this if-s "level" should be adjusted perhaps >= SAVE_FILES_FINAL, also get rid of m_pSkyImageReal / m_pSkyImageImag
       // 2024-06-23 : these images are not saved in GPU version. In CPU version there were after 2D FFT and normalisation but before FFTShift. 
       //              However, in GPU version 2D FFT is followed by FFTShift+normalisation in one kernel (no DeviceToHost copy in between) -> cannot save this product 
       /*if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_FINAL ){
          for(int i = 0; i < image_size; i++){
            out_data_real[i] = ((gpufftComplex*)m_out_data)[i].x;
            out_data_imag[i] = ((gpufftComplex*)m_out_data)[i].y;
          }
        } 

        if( bSaveIntermediate ){ // I will keep this if - assuming it's always TRUE, but there is still control using , if bSaveIntermediate=false it has priority over m_SaveFilesLevel
           if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
              sprintf(outDirtyImageReal,"%s/dirty_test_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),uv_grid_counter_xSize,uv_grid_counter_ySize);
              sprintf(outDirtyImageImag,"%s/dirty_test_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),uv_grid_counter_xSize,uv_grid_counter_ySize);
   
              out_image_real.WriteFits( outDirtyImageReal );
              out_image_imag.WriteFits( outDirtyImageImag );
           }
         }*/
   
        // 2022-04-02 : test change to use member variable for final image (have to be careful with threads and to not use this class as global variable):
        // calculate and save FFT-shifted image :
        // CBgFits out_image_real2( out_image_real.GetXSize(), out_image_real.GetYSize() ), out_image_imag2( out_image_real.GetXSize(), out_image_real.GetYSize() );
        AllocOutPutImages( m_uv_grid_real->GetXSize(), m_uv_grid_real->GetYSize()  );

        char outDirtyImageReal[1024],outDirtyImageImag[1024];   
        if( !m_pSkyImageReal || !m_pSkyImageImag )
        {
           printf("ERROR in code : internal image buffers not allocated -> cannot continue\n");
           return;
        }

        // 2022-04-02 : test change to use member variable for final image (have to be careful with threads and to not use this class as global variable):
        // TODO : CPU -> GPU 
        // fft_shift( out_image_real, *m_pSkyImageReal );
        // fft_shift( out_image_imag, *m_pSkyImageImag );
   
        int rest = 1; // just so that by default it is !=0 -> image not saved 
        if( CPacerImager::m_SaveControlImageEveryNth > 0 )
        {
           rest = (m_SkyImageCounter % CPacerImager::m_SaveControlImageEveryNth);
           if( rest == 0 )
           {
              PRINTF_INFO("INFO : saving %d-th control sky image\n",m_SkyImageCounter);
           }
        }

         if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_FINAL || rest==0 )
         {   
            // 2024-06-22 : temporary code m_pSkyImageReal and m_pSkyImageImag are no longer needed as normalised and FFT shifted image is already in out_image_real / out_image_imag
            float* sky_data_real = m_pSkyImageReal->get_data();
            float* sky_data_imag = m_pSkyImageImag->get_data();
            for(int i = 0; i < image_size; i++){
               sky_data_real[i] = ((gpufftComplex*)m_out_data)[i].x; // was in CPU FFTSHIFT version : out_data_real[i];
               sky_data_imag[i] = ((gpufftComplex*)m_out_data)[i].y; // was in CPU FFTSHIFT version : out_data_imag[i];
            }
   
            if( szBaseOutFitsName && strlen(szBaseOutFitsName) ){
               sprintf(outDirtyImageReal,"%s/%s_real.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
            }else{
               // sprintf(outDirtyImageReal,"dirty_test_real_fftshift_%dx%d.fits",width,height);
               // const char* get_filename(  time_t ut_time , char* out_buffer, int usec=0, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix="", const char* formater="%.2u%.2u%.2uT%.2u%.2u%.2u" );
               get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageReal, m_ImagerParameters.m_szOutputDirectory.c_str(), "dirty_image_", "_real" ); // uxtime=0 -> it will be taken as current system time
            }
            SaveSkyImage( outDirtyImageReal , m_pSkyImageReal );
            PRINTF_DEBUG("Saved read file to %s\n",outDirtyImageReal);
   
            if( bSaveImaginary ){
               if( szBaseOutFitsName && strlen(szBaseOutFitsName) )
               {
                  sprintf(outDirtyImageImag,"%s/%s_imag.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
               }else{
                  // sprintf(outDirtyImageImag,"dirty_test_imag_fftshift_%dx%d.fits",width,height);
                  get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageImag, m_ImagerParameters.m_szOutputDirectory.c_str(), "dirty_image_", "_imag" );
               }

               m_pSkyImageImag->SetFileName( outDirtyImageImag );      
               m_pSkyImageImag->WriteFits( outDirtyImageImag );
               PRINTF_DEBUG("Saved imaginary file to %s\n",outDirtyImageImag);
            }
         }
   
         if( CPacerImager::m_bPrintImageStatistics ){
            double mean, rms, minval, maxval, median, iqr, rmsiqr;
            int cnt;
            int radius = int( sqrt( m_pSkyImageReal->GetXSize()*m_pSkyImageReal->GetXSize() + m_pSkyImageReal->GetYSize()*m_pSkyImageReal->GetYSize() ) ) + 10;
            // m_SkyImageReal.GetStat( mean, rms, minval, maxval );
            m_pSkyImageReal->GetStatRadiusAll( mean, rms, minval, maxval, median, iqr, rmsiqr, cnt, radius, true );
            printf("STAT : full image %s statistics in radius = %d around the center using %d pixels : mean = %.6f , rms = %.6f, minval = %.6f, maxval = %.6f, median = %.6f, rms_iqr = %.6f\n",outDirtyImageReal,radius,cnt,mean, rms, minval, maxval, median, rmsiqr );
      
      
            // TODO : this will be parameterised to specified requested windown in the image to get RMS value from:
            double mean_window, rms_window, minval_window, maxval_window, median_window, iqr_window, rmsiqr_window;
            radius = 10; // TODO : make it use parameter and also position in the image 
            m_pSkyImageReal->GetStatRadiusAll( mean_window, rms_window, minval_window, maxval_window, median_window, iqr_window, rmsiqr_window, cnt, radius, true );
            printf("STAT : statistics of %s in radius = %d around the center using %d pixels : mean = %.6f , rms = %.6f, minval = %.6f, maxval = %.6f, median = %.6f, rms_iqr = %.6f\n",outDirtyImageReal,radius,cnt,mean_window, rms_window, minval_window, maxval_window, median_window, rmsiqr_window );
         }
   
         // TODO : re-grid to SKY COORDINATES !!!
         // convert cos(alpha) to alpha - see notes !!!
         // how to do it ???
   }else{ 
      printf("WARNING : saving temporary FITS files is disabled. Hence, statistics is not printed either. Enable with option -V 100 (or lower)\n");
   }
}

// TODO : 
//     - do more cleanup of this function as this is nearly "RAW" copy paste from Gayatri's code:
//     - optimise (same in CPU version) uv_grid_counter, uv_grid_real, uv_grid_imag to become member variables.
void CPacerImagerHip::gridding_imaging( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
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
  PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f (20230929)\n",min_uv);

  // allocates data structures for gridded visibilities:
  AllocGriddedVis( n_pixels, n_pixels );


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
  int u_xSize = fits_vis_u.GetXSize();
  int u_ySize = fits_vis_u.GetYSize();

  int vis_real_xSize =  fits_vis_real.GetXSize(); 
  int vis_real_ySize =  fits_vis_real.GetYSize(); 

  int vis_imag_xSize =  fits_vis_imag.GetXSize(); 
  int vis_imag_ySize =  fits_vis_real.GetYSize(); 



  // Output size: uv_grid_real, uv_grid_imag, uv_grid_counter 
  int uv_grid_counter_xSize = m_uv_grid_counter->GetXSize();
  int uv_grid_counter_ySize = m_uv_grid_counter->GetYSize();

  int xySize = (u_xSize*u_ySize);
  int image_size = (uv_grid_counter_xSize*uv_grid_counter_ySize); 
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

  double frequency_hz = frequency_mhz*1e6;
  double wavelength_m = VEL_LIGHT / frequency_hz;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL)
  {  
     PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg]\n",pixscale_zenith_deg,pixscale_deg_version2);
     m_PixscaleAtZenith = pixscale_deg_version2; 
     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  int center_x = int(n_pixels/2);
  int center_y = int(n_pixels/2);
  int is_odd_x = 0 , is_odd_y = 0;
  if( (n_pixels % 2) == 1 ){
     is_odd_x = 1;
  }
  if( (n_pixels % 2) == 1 ){
     is_odd_y = 1;
  }

  // initialise values in the UV GRID to zeros: 
  // 2024-06-14 : this code is not need as it is only used to initialise GPU arrays with zeros later in the code (cudaMemcpy), but these arrays are already initialised with zeros in 
  //              AllocGPUMemory function using cudaMemset
  // m_uv_grid_real->SetZeroValue();
  // m_uv_grid_imag->SetZeroValue();
  // m_uv_grid_counter->SetZeroValue();

  // Setting the initial values of out_image_real/out_image_imag 
  // MS (2024-06-12) : not required as these are overwritten later in a loop where normalisation factor fnorm is applied:
  // out_image_real.SetZeroValue();
  // out_image_imag.SetZeroValue();

  // Step 1: Declare GPU(Device) and CPU(Host) Variables 
  // CPU input variables 
  float *u_cpu = fits_vis_u.get_data();
  float *v_cpu = fits_vis_v.get_data();
  float *w_cpu = fits_vis_w.get_data();
  float *vis_real_cpu = fits_vis_real.get_data();
  float *vis_imag_cpu = fits_vis_imag.get_data();

  // GPU Variables: (Corresponding GPU Variables) declared in class 

  PACER_PROFILER_RESTART
  // Allocate only if not-allocated 
  // TODO : warning GPU UV grid is not initialised to ZEROs :
  AllocGPUMemory(xySize, image_size ); //  out_image_real.get_data() );
  
  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n GRIDDING CHECK: Step 2 Memory allocated for GPU variables");
     printf("\n GRIDDING CHECK: Step 2 GPU Variables initialised: AllocGPUMemory() called");
  }
  // End of gpuMalloc() 
  PACER_PROFILER_END("GPU Memory allocations took")

  // Step 3: Copy contents from CPU to GPU [input variables]
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy()
  PACER_PROFILER_RESTART
  (gpuMemcpy((float*)u_gpu, (float*)u_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)v_gpu, (float*)v_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)w_gpu, (float*)w_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice));
  (gpuMemcpy((float*)vis_real_gpu, (float*)vis_real_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)vis_imag_gpu, (float*)vis_imag_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemset((float*)uv_grid_counter_gpu, 0, image_size*sizeof(float))); // UV grid has to be set to zero
  (gpuMemset((float*)uv_grid_real_gpu, 0, image_size*sizeof(float)));    // UV grid has to be set to zero
  (gpuMemset((float*)uv_grid_imag_gpu, 0, image_size*sizeof(float)));    // UV grid has to be set to zero
  (gpuMemset((gpufftComplex*)m_in_buffer_gpu, 0, image_size*sizeof(gpufftComplex))); // UV grid has to be set to zero
  PACER_PROFILER_SHOW("GPU Memory copy host to device took")
  
  int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n GRIDDING CHECK: Step 3 CPU to GPU copied"); 

     printf("\n GRIDDING CHECK: NTHREADS = %d", NTHREADS);
     printf("\n GRIDDING CHECK: nBlocks = %d", nBlocks);
  }
  
  // update antenna flags before gridding which uses these flags or weights:
  int n_ant = fits_vis_u.GetXSize();
  UpdateAntennaFlags( n_ant );

  // Step 4: Call to GPU kernel
  
  // Start of kernel call 
  PACER_PROFILER_START_TIMER1 // just to measure gridding and imaging (excluding memory allocation and copying)
  PACER_PROFILER_RESTART
  gridding_imaging_cuda<<<nBlocks,NTHREADS>>>( xySize, n_ant, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, wavelength_m, image_size, delta_u, delta_v, n_pixels, center_x, center_y, is_odd_x, is_odd_y, vis_real_gpu, vis_imag_gpu, uv_grid_counter_gpu, uv_grid_real_gpu, uv_grid_imag_gpu, min_uv, (gpufftComplex*)m_in_buffer_gpu ); 
  PRINTF_DEBUG("\n GRIDDING CHECK: Step 4 Calls to kernel");
  
  // Gives the error in the kernel! 
  (gpuGetLastError());
  (gpuDeviceSynchronize());
  PACER_PROFILER_SHOW("Gridding kernel call took")
  
  // TODO: 2024-06-22 : DIVIDE m_in_buffer_gpu and uv_grid_real_gpu, uv_grid_imag_gpu by uv_grid_counter_gpu for uniform and other weightings to really work
  //            are uv_grid_imag_gpu uv_grid_real_gpu really needed ???

  // Checking Execution time for cuFFT 
  PACER_PROFILER_RESTART
  if( !m_FFTPlan ){
     gpufftPlan2d((gpufftHandle*)(&m_FFTPlan), uv_grid_counter_xSize, uv_grid_counter_ySize, GPUFFT_C2C);
     PRINTF_INFO("INFO : gpufftPlan2d created\n");
  }
  gpufftExecC2C(((gpufftHandle)m_FFTPlan), (gpufftComplex*)m_in_buffer_gpu, (gpufftComplex*)m_out_buffer_gpu, GPUFFT_FORWARD);

  // End of cuFFT 
  PACER_PROFILER_SHOW("cuFFT execution took")
  PACER_PROFILER_SHOW_TIMER1("gridding and cuFFT (total) took")

  // Step 5: Copy contents from GPU variables to CPU variables
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy() 
  PACER_PROFILER_RESTART
  
  // TODO: CPU->GPU : calculate this sum on GPU, can it be done in the gridding kernel itself ???
  // double fnorm = 1.00/m_uv_grid_counter->Sum();
  // TODO : ??? this can also be done in gridding_imaging_cuda<<<nBlocks,NTHREADS>>> , right, will not reduce efficiency if another atomicAdd is used for the sum ???
  double fnorm = 1.00/sum_gpu_atomicadd( uv_grid_counter_gpu, image_size );
  // double fnorm = 1.00/sum_gpu_parallel_reduce( uv_grid_counter_gpu, image_size );
  
  
  // TODO : for now keeping it as it was as there is no clear advantage of doing normalsation on GPU 
  // apply normalisation constant on GPU :
  // int nBlocksImage = (image_size + NTHREADS -1)/NTHREADS;
  // mult_by_const<<<nBlocksImage,NTHREADS>>>( (gpufftComplex*)m_out_buffer_gpu, image_size, fnorm );
  
  // FFT shift together with multiplication by fnorm (normalisation)
  // bool fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, float fnorm=1.00 );
  fft_shift_and_norm_gpu( (gpufftComplex*)m_out_buffer_gpu, m_uv_grid_real->GetXSize(), m_uv_grid_real->GetYSize(), fnorm );

  // End of gpuMemcpy() GPU to CPU 
  PACER_PROFILER_SHOW("GPU memory copy from device to host took")
  
  // save test FITS files and print STAT (if file debug level is high enough)
  SaveTestFitsFilesAndShowStat( n_pixels, weighting, szBaseOutFitsName, bSaveIntermediate, bSaveImaginary );
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
  PRINTF_DEBUG("DEBUG : gridding_imaging_cuda_version : min_uv = %.4f\n",min_uv);

  // allocates data structures for gridded visibilities:
  AllocGriddedVis( n_pixels, n_pixels );


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
  int u_xSize = fits_vis_u.GetXSize();
  int u_ySize = fits_vis_u.GetYSize();

  int n_ant = xcorr.obsInfo.nAntennas;
  int vis_real_xSize = n_ant; 
  int vis_real_ySize = n_ant; 

  // Output size: uv_grid_real, uv_grid_imag, uv_grid_counter 
  int uv_grid_counter_xSize = m_uv_grid_counter->GetXSize();
  int uv_grid_counter_ySize = m_uv_grid_counter->GetYSize();

  int xySize = (u_xSize*u_ySize);
  int image_size = (uv_grid_counter_xSize*uv_grid_counter_ySize); 
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

  double frequency_hz = frequency_mhz*1e6;
  double wavelength_m = VEL_LIGHT / frequency_hz;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL)
  {  
     PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg]\n",pixscale_zenith_deg,pixscale_deg_version2);
     m_PixscaleAtZenith = pixscale_deg_version2; 
     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  int center_x = int(n_pixels/2);
  int center_y = int(n_pixels/2);
  int is_odd_x = 0 , is_odd_y = 0;
  if( (n_pixels % 2) == 1 ){
     is_odd_x = 1;
  }
  if( (n_pixels % 2) == 1 ){
     is_odd_y = 1;
  }

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
//  float *vis_real_cpu = fits_vis_real.get_data();
//  float *vis_imag_cpu = fits_vis_imag.get_data();

  // CPU Output variables
  // MS (2024-06-14) I leave these for know, but in fact they should only be used/executed if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG )
  //                 because otherwise we do not need to have all these data copied from GPU to CPU
  // TODO : add if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ) -> will require more testing
  float *uv_grid_real_cpu = m_uv_grid_real->get_data();
  float *uv_grid_imag_cpu = m_uv_grid_imag->get_data();
  float *uv_grid_counter_cpu = m_uv_grid_counter->get_data();

  // GPU Variables: (Corresponding GPU Variables) declared in class 

  // Calculating gpuMalloc() 
  clock_t start_time2 = clock();

  // Step 2: Allocate memory for GPU variables
  // Memory for GPU input variables: 
  // (gpuMalloc((float**)&vis_real_gpu, xySize*sizeof(float)));
  // (gpuMalloc((float**)&vis_imag_gpu, xySize*sizeof(float)));
  // (gpuMalloc((float**)&u_gpu, xySize*sizeof(float)));
  // (gpuMalloc((float**)&v_gpu, xySize*sizeof(float)));

  // Allocate only if not-allocated 

  // TODO : warning GPU UV grid is not initialised to ZEROs :
  AllocGPUMemory(xySize, image_size ); //  out_image_real.get_data() );
 
  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n GRIDDING CHECK: Step 2 Memory allocated for GPU variables");
     printf("\n GRIDDING CHECK: Step 2 GPU Variables initialised: AllocGPUMemory() called");
  }

  // End of gpuMalloc() 
  clock_t end_time2 = clock();
  double duration_sec2 = double(end_time2-start_time2)/CLOCKS_PER_SEC;
  double duration_ms2 = duration_sec2*1000;
  printf("\n ** CLOCK gpuMalloc() took : %.6f [seconds], %.3f [ms]\n",duration_sec2,duration_ms2);

  // Step 3: Copy contents from CPU to GPU [input variables]
  // gpuMemcpy(destination, source, size, HostToDevice)

  // Start of gpuMemcpy()
  clock_t start_time3 = clock();

  (gpuMemcpy((float*)u_gpu, (float*)u_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)v_gpu, (float*)v_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemcpy((float*)w_gpu, (float*)w_cpu, sizeof(float)*xySize, gpuMemcpyHostToDevice)); 
  (gpuMemset((float*)uv_grid_counter_gpu, 0, image_size*sizeof(float))); // UV grid has to be set to zero
  (gpuMemset((float*)uv_grid_real_gpu, 0, image_size*sizeof(float)));    // UV grid has to be set to zero
  (gpuMemset((float*)uv_grid_imag_gpu, 0, image_size*sizeof(float)));    // UV grid has to be set to zero
  (gpuMemset((gpufftComplex*)m_in_buffer_gpu, 0, image_size*sizeof(gpufftComplex))); // UV grid has to be set to zero

  // TODO : COPY xcorr strucuttre here:
  VISIBILITY_TYPE* vis_local_gpu = NULL;
 
  // vis_gpu=NULL -> data is already on GPU (otherwise error in code !)
  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n GRIDDING CHECK: Step 2 getting GPU pointer to visibilities (CORRECTED)\n");
  }
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
     // save correlation matrix before any correcitons are done:
     CBgFits tmp_real(n_ant,n_ant),tmp_imag(n_ant,n_ant);
     ConvertXCorr2Fits( xcorr, tmp_real, tmp_imag, time_step, fine_channel, "gridding_imaging_gpu_precorr" );
  }
  
  if( xcorr.on_gpu() ){
     // xcorr.data() is on GPU -> setting the pointed
     
     // WARNING : this code does not use time/frequency:
     // vis_local_gpu = (VISIBILITY_TYPE*)xcorr.data();         
     // NEW CODE - 2024-04-07 :
     vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
  }else{
     // data in xcorr is in CPU 
     // this is only allocated when data in xcorr structure is on CPU then to run GPU imager (which it is) we need to copy from CPU to GPU :a
     if( !vis_gpu ){
        // if not already allocated -> allocate required GPU memory
        AllocGPUMemoryForXCorr( &xcorr );
     }
     // set the pointer and copy data from Host (in xcorr) to device:
     vis_local_gpu = vis_gpu;    
     VISIBILITY_TYPE* corr_matrix_ptr = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
//     (gpuMemcpy((VISIBILITY_TYPE*)vis_local_gpu, (VISIBILITY_TYPE*)xcorr.data(), sizeof(std::complex<VISIBILITY_TYPE>) * xcorr.size(), gpuMemcpyHostToDevice));
     (gpuMemcpy((VISIBILITY_TYPE*)vis_local_gpu, (VISIBILITY_TYPE*)corr_matrix_ptr, sizeof(std::complex<VISIBILITY_TYPE>) * xcorr.matrix_size(), gpuMemcpyHostToDevice));

     // TODO: frequency / time is not acconted for here
     // printf("ERROR in CODE : frequency / time is not accounted in this version of the code (hybrid - imagerGPU and blinkpipelineCPU)\n");     
  }
  
  int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
  if( m_ImagerParameters.m_bApplyGeomCorr ){ // turned off for testing
     // Apply geometric correction here :
     apply_geometric_corrections<<<nBlocks,NTHREADS>>>( xySize, n_ant, vis_local_gpu, w_gpu, frequency_hz, SPEED_OF_LIGHT );
     printf("DEBUG : after call of apply_geometric_corrections kernel\n");
     
     if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){ // TODO : remove in realese version:
        // Gives the error in the kernel! 
        (gpuGetLastError());
        (gpuDeviceSynchronize());
     
        // __global__ void vis2corrmatrix( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *vis_corrmatrix_re_cuda, float* vis_corrmatrix_im_cuda );
        vis2corrmatrix<<<nBlocks,NTHREADS>>>( xySize, n_ant, vis_local_gpu, test_data_real_gpu, test_data_imag_gpu );


        CBgFits tmp_real(n_ant,n_ant),tmp_imag(n_ant,n_ant);
        (gpuMemcpy((float*)tmp_real.get_data(), (float*)test_data_real_gpu, sizeof(float)*xySize, gpuMemcpyDeviceToHost));
        (gpuMemcpy((float*)tmp_imag.get_data(), (float*)test_data_imag_gpu, sizeof(float)*xySize, gpuMemcpyDeviceToHost));
        
        tmp_real.WriteFits("vis_re_geom_corr_gpu.fits");
        tmp_imag.WriteFits("vis_im_geom_corr_gpu.fits");
     }
  }

  if( m_ImagerParameters.m_bApplyCableCorr ){ // turned off for testing 
     // init cable length buffer for cuda :
     InitCableLengts( n_ant );
  
     // Apply cable correction here :     
     apply_cable_corrections<<<nBlocks,NTHREADS>>>( xySize, n_ant, vis_local_gpu, cable_lengths_gpu, frequency_hz, SPEED_OF_LIGHT );
     printf("DEBUG : after call of apply_cable_corrections kernel\n");

     if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){ // TODO : remove in realese version:
        // Gives the error in the kernel! 
        (gpuGetLastError());
        (gpuDeviceSynchronize());
     
        // __global__ void vis2corrmatrix( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *vis_corrmatrix_re_cuda, float* vis_corrmatrix_im_cuda );
        vis2corrmatrix<<<nBlocks,NTHREADS>>>( xySize, n_ant, vis_local_gpu, test_data_real_gpu, test_data_imag_gpu );


        // WARNING : strange code with extra memcpy as otherwise it crashes on hipMemcpy ???!!!
/*        printf("DEBUG : values %d x %d = %d vs. %d\n",n_ant,n_ant,(n_ant*n_ant),xySize);
        CBgFits tmp_real(n_ant,n_ant),tmp_imag(n_ant,n_ant);
        // CBgFits* p_tmp_real= new CBgFits(n_ant,n_ant);
        float* test_ptr = new float[n_ant*n_ant];
        //(hipMemcpy(test_ptr, (float*)test_data_real_gpu, sizeof(float)*n_ant*n_ant, hipMemcpyDeviceToHost));
        (hipMemcpy((float*)test_ptr, (float*)test_data_real_gpu, sizeof(float)*n_ant*n_ant, hipMemcpyDeviceToHost));
        memcpy(tmp_real.get_data(),test_ptr,sizeof(float)*n_ant*n_ant);
        //(hipMemcpy((float*)p_tmp_real->get_data(), (float*)test_data_real_gpu, sizeof(float)*xySize, hipMemcpyDeviceToHost));
        (hipMemcpy((float*)test_ptr, (float*)test_data_imag_gpu, sizeof(float)*n_ant*n_ant, hipMemcpyDeviceToHost));
        memcpy(tmp_imag.get_data(),test_ptr,sizeof(float)*n_ant*n_ant);
        printf("DEBUG : crashed ???\n");
        delete [] test_ptr;
*/        
        CBgFits tmp_real(n_ant,n_ant),tmp_imag(n_ant,n_ant);
        (gpuMemcpy((float*)tmp_real.get_data(), (float*)test_data_real_gpu, sizeof(float)*n_ant*n_ant, gpuMemcpyDeviceToHost));
        (gpuMemcpy((float*)tmp_imag.get_data(), (float*)test_data_imag_gpu, sizeof(float)*n_ant*n_ant, gpuMemcpyDeviceToHost));
                
        tmp_real.WriteFits("vis_re_cable_corr_gpu.fits");
        tmp_imag.WriteFits("vis_im_cable_corr_gpu.fits");
        printf("DEBUG : saved vis_re_cable_corr_gpu.fits and vis_im_cable_corr_gpu.fits (CPacerImagerHip::gridding_imaging)\n");
     }

  }
  
  // update antenna flags before gridding which uses these flags or weights:
  UpdateAntennaFlags( n_ant );

  clock_t end_time3 = clock();
  double duration_sec3 = double(end_time3-start_time3)/CLOCKS_PER_SEC;
  double duration_ms3 = duration_sec3*1000;
  printf("\n ** CLOCK gpuMemcpy() CPU to GPU took : %.6f [seconds], %.3f [ms]\n",duration_sec3,duration_ms3);

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n GRIDDING CHECK: Step 3 CPU to GPU copied"); 

     printf("\n GRIDDING CHECK: NTHREADS = %d", NTHREADS);
     printf("\n GRIDDING CHECK: nBlocks = %d", nBlocks);
  }

  // Step 4: Call to GPU kernel
  
  // Start of kernel call 
  clock_t start_time4 = clock();

  gridding_imaging_cuda_xcorr<<<nBlocks,NTHREADS>>>( xySize, n_ant, u_gpu, v_gpu, antenna_flags_gpu, antenna_weights_gpu, wavelength_m, image_size, delta_u, delta_v, n_pixels, center_x, center_y, is_odd_x, is_odd_y, vis_local_gpu, uv_grid_counter_gpu, uv_grid_real_gpu, uv_grid_imag_gpu, min_uv, (gpufftComplex*)m_in_buffer_gpu ); 
  PRINTF_DEBUG("\n GRIDDING CHECK: Step 4 Calls to kernel");
  
  // End of kernel call 
  clock_t end_time4 = clock();
  double duration_sec4 = double(end_time4-start_time4)/CLOCKS_PER_SEC;
  double duration_ms4 = duration_sec4*1000;
  printf("\n ** CLOCK kernel call took : %.6f [seconds], %.3f [ms]\n",duration_sec4,duration_ms4);

  // Gives the error in the kernel! 
  (gpuGetLastError());
  (gpuDeviceSynchronize());
  
  // TODO: 2024-06-22 : DIVIDE m_in_buffer_gpu and uv_grid_real_gpu, uv_grid_imag_gpu by uv_grid_counter_gpu for uniform and other weightings to really work
  //            are uv_grid_imag_gpu uv_grid_real_gpu really needed ???


  // TEMPORARY - check if UV Grid values are already truncated straight after after the kernel call
  // (gpuMemcpy((float*)uv_grid_real_cpu, (float*)uv_grid_real_gpu, sizeof(float)*image_size, gpuMemcpyDeviceToHost)); 
  // printf("\nDEBUG : GPU gridding (4,0) = %.20f [just after kernel call gridding_imaging_cuda]\n",m_uv_grid_real->getXY(4,0));  


  // uv_grid_counter_xSize = width
  // uv_grid_counter_ySize = height
  // size = image_size: (width x height)

  // Checking Execution time for cuFFT 
  clock_t start_time6 = clock();

  if( !m_FFTPlan ){
     gpufftPlan2d((gpufftHandle*)(&m_FFTPlan), uv_grid_counter_xSize, uv_grid_counter_ySize, GPUFFT_C2C);
     PRINTF_INFO("INFO : gpufftPlan2d created\n");
  }
  gpufftExecC2C(((gpufftHandle)m_FFTPlan), (gpufftComplex*)m_in_buffer_gpu, (gpufftComplex*)m_out_buffer_gpu, GPUFFT_FORWARD);

  // End of cuFFT 
  clock_t end_time6 = clock();
  double duration_sec6 = double(end_time6-start_time6)/CLOCKS_PER_SEC;
  double duration_ms6 = duration_sec6*1000;
  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
     printf("\n ** CLOCK cuFFT() took : %.6f [seconds], %.3f [ms]\n",duration_sec6,duration_ms6);
     printf("\n Imaging CHECK: cuFFT completed: \n "); 
  }

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
  fft_shift_and_norm_gpu( (gpufftComplex*)m_out_buffer_gpu, m_uv_grid_real->GetXSize(), m_uv_grid_real->GetYSize(), fnorm );
  
  // End of gpuMemcpy() GPU to CPU 
  clock_t end_time5 = clock();
  double duration_sec5 = double(end_time5-start_time5)/CLOCKS_PER_SEC;
  double duration_ms5 = duration_sec5*1000;
  printf("\n ** CLOCK gpuMemcpy() GPU to CPU took : %.6f [seconds], %.3f [ms]\n",duration_sec5,duration_ms5);
  PRINTF_DEBUG("\n GRIDDING CHECK: Step 5 GPU to CPU copied"); 

  // save test FITS files and print STAT (if file debug level is high enough)
  SaveTestFitsFilesAndShowStat( n_pixels, weighting, szBaseOutFitsName, bSaveIntermediate, bSaveImaginary );
}

bool CPacerImagerHip::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz, int time_step, int fine_channel )
{
   printf("DEBUG: CPacerImagerHip::ApplyGeometricCorrections is a void function just NOT TO DO THIS IN CPU, but using kernel : apply_geometric_corrections (called in gridding_imaging)\n");
   return false;
}

bool CPacerImagerHip::ApplyCableCorrections( Visibilities& xcorr, double frequency_mhz, int time_step, int fine_channel )
{
   printf("DEBUG: CPacerImagerHip::ApplyCableCorrections is a void function just NOT TO DO THIS IN CPU, but using kernel : apply_cable_corrections (called in gridding_imaging)\n");
   return false;
}

//
// MS : let's keep this function commented and use it in case I need to test if "make test" works ok with HIP version using CPU gridding and FFTW
//
// TEST : use the same version as in CPacerImager which should pass both tests now !!!
// void CPacerImagerHip::gridding_imaging( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
//                                     double delta_u, double delta_v,
//                                     double frequency_mhz,
//                                     int    n_pixels,
//                                     double min_uv /*=-1000*/,    // minimum UV 
//                                     const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
//                                     const char* szBaseOutFitsName /*=NULL*/, 
//                                     bool do_gridding,
//                                     bool do_dirty_image,
//                                     const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
//                                     const char* in_fits_file_uv_im  /*=""*/ // gridded visibilities can be provided externally
//                )
/*{
  // TODO: make it member variable to optimise and do not allocate every time !!!
  CBgFits uv_grid_counter( n_pixels, n_pixels ),uv_grid_real( n_pixels, n_pixels ) , uv_grid_imag( n_pixels, n_pixels );  
  
  if( do_gridding ){
     gridding_fast( fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, uv_grid_real, uv_grid_imag, uv_grid_counter, delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting );
  }else{
     if( strlen(in_fits_file_uv_re) && strlen(in_fits_file_uv_im) ){
        uv_grid_counter.SetValue(1.00);
        
        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_re);
        if( uv_grid_real.ReadFits( in_fits_file_uv_re, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_re);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_re);
        }

        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_im);
        if( uv_grid_real.ReadFits( in_fits_file_uv_im, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_im);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_im);
        }

     }else{
        printf("ERROR : when gridding is disabled (-g 0) options -r and -i with REAL and IMAG FITS file names must be specified -> cannot continue !\n");
        exit(-1);
     }
  }

  if( do_dirty_image ){
     // dirty image :  
     PRINTF_INFO("PROGRESS : executing dirty image\n");
     dirty_image( uv_grid_real, uv_grid_imag, uv_grid_counter, 
                  true, // do not save intermediate FITS files
                  szBaseOutFitsName, // output filename template
                  true,   // save imaginary image
                  false    // do FFTunshift true when gridding() , false for gridding_fast()
                );             
  }


}                

*/
