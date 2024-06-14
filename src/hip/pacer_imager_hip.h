#ifndef _PACER_IMAGER_GPU_H__
#define _PACER_IMAGER_GPU_H__

#include <typeinfo>

#include "../pacer_imager.h"

class CPacerImagerHip : public CPacerImager 
{
protected :
   // CUDA / HIP FFT plan:
   long int m_FFTPlan; // initialised on the first usage , WARNING : on Setonix was void* , int is too small and causes CORE DUMP CRASH !!!

   // member variables for cuFFT :
   void* m_in_buffer_gpu;  // buffer for gridded visibilities which are passed to cuFFT/hipFFT

   void* m_out_buffer_gpu; // output of cu/hip FFT - GPU memory (Device)
   void* m_out_buffer_cpu; // output of cu/fip FFT - Host memory (Host)
   void* m_out_data; // temporary CPU buffer to apply normalisation. TODO : apply normalisation in a GPU kernal


   // Additional GPU Input variables 
   // GPU Input variables 
   float *u_gpu; 
   float *v_gpu;
   float *w_gpu;
   float *vis_real_gpu; 
   float *vis_imag_gpu; 
   VISIBILITY_TYPE* vis_gpu;
   float *cable_lengths_gpu;
   float *cable_lengths_cpu;
  
   // GPU Output variables 
   float *uv_grid_real_gpu;
   float *uv_grid_imag_gpu;
   float *uv_grid_counter_gpu; 
   float *uv_grid_counter_cpu; // temporarily counter is on CPU to be able to save and calculat Sum -> later paralllel reduction on GPU 
   
   // antenna flags (0-ok, 1-flagged) and weights (1-ok, 0-remove)
   int *antenna_flags_gpu;
   int *antenna_flags_cpu;
   float *antenna_weights_gpu;
   float *antenna_weights_cpu;
   
   int m_AllocatedXYSize;    // size of allocated correlation matrix or UVW 
   int m_AllocatedImageSize; // size of memory allocated for images 
   
   // Test buffers:
   float *test_data_real_gpu; // intitialised for size of corr. matrix 
   float *test_data_imag_gpu; // intitialised for size of corr. matrix 

   // Memory management functions for GPU version - for just a single image (N parameter is not used here on purpose)
   //  N blocks version is implemented in the derived class CPacerImagerMultiHip
   virtual void AllocGPUMemory( int corr_size, 
                                int image_size
                              );
        
   // allocate GPU memory specifically for visibilities in xcorr structure - this is only for the case where xcorr structure has it allocated on CPU (for whatever reason)
   void AllocGPUMemoryForXCorr( Visibilities* p_xcorr );                              

   // Clean GPU Memory 
   virtual void CleanGPUMemory(); 
   
   // initialisation of specific ararys:
   void InitCableLengts( int n_ant );
   
   // update antenna flags:
   void UpdateAntennaFlags( int n_ant );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 
   // GOAL   : function calling both gridding_fast and imaging. This function is virtual and will be overwritten in HIP version to do GPU specific code
   // 
   // INPUT  : 
   //          fits_vis_real, fits_vis_imag : visibilities (REAL and IMAG 2D arrays as FITS class) 
   //          fits_vis_u, fits_vis_v, fits_vis_w : UVW (real values baselines in units of wavelength - see TMS)
   //          delta_u, delta_v : size of the UV cell 
   //          frequency_mhz : frequency in MHz
   //
   // OUTPUT : 
   //          - uv_grid_real, uv_grid_imag : visibilities on UV grid (real and imag arrays)
   //          - uv_grid_counter : visibility counter and 
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual void gridding_imaging( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  double delta_u, double delta_v,
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL,
                  bool do_gridding=true,                  
                  bool do_dirty_image=true,
                  const char* in_fits_file_uv_re = "", // gridded visibilities can be provided externally
                  const char* in_fits_file_uv_im = "", // gridded visibilities can be provided externally
                  bool bSaveIntermediate=false , bool bSaveImaginary=true 
                );

   virtual void gridding_imaging( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  double delta_u, double delta_v,
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL,
                  bool do_gridding=true,                  
                  bool do_dirty_image=true,
                  const char* in_fits_file_uv_re = "", // gridded visibilities can be provided externally
                  const char* in_fits_file_uv_im = "", // gridded visibilities can be provided externally
                  bool bSaveIntermediate=false , bool bSaveImaginary=true 
                );

    // virtual function to NOT DO corrections in CPU but in GPU :
    virtual bool ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz );
   
    virtual bool ApplyCableCorrections( Visibilities& xcorr, double frequency_mhz );


public :
   //-----------------------------------------------------------------------------------------------------------------------------
   // IsGPU() - returns true if CPacerImagerHip object and false here:
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual inline bool IsGPU(){ return true; }

   CPacerImagerHip();
   ~CPacerImagerHip();
};


#endif