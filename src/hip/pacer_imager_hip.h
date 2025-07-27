#ifndef _PACER_IMAGER_GPU_H__
#define _PACER_IMAGER_GPU_H__

#include <typeinfo>
#include <gpu_fft.hpp>
#include <astroio.hpp>

#include "../pacer_imager.h"

class CPacerImagerHip : public CPacerImager {
protected :
   MemoryBuffer<double> frequencies_gpu;
   MemoryBuffer<float> fnorm;
   MemoryBuffer<float> u_gpu;
   MemoryBuffer<float> v_gpu;
   MemoryBuffer<float> w_gpu;
   MemoryBuffer<float> cable_length_gpu;
   
   // CUDA / HIP FFT plan:
   gpufftHandle m_FFTPlan; // initialised on the first usage , WARNING : on Setonix was void* , int is too small and causes CORE DUMP CRASH !!!

   // Additional GPU Input variables 
   // GPU Input variables 
   VISIBILITY_TYPE* vis_gpu;
   float *cable_lengths_gpu;
   float *cable_lengths_cpu;
  
   // antenna flags (0-ok, 1-flagged) and weights (1-ok, 0-remove)
   int *antenna_flags_gpu;
   int *antenna_flags_cpu;
   float *antenna_weights_gpu;
   float *antenna_weights_cpu;
   
   // int m_AllocatedXYSize;    // size of allocated correlation matrix or UVW 
   // int m_AllocatedImageSize; // size of memory allocated for images 
   

   // Clean GPU Memory 
   void CleanGPUMemory(); 
   
   
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
   virtual Images gridding_imaging( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  double delta_u, double delta_v,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL
                );

    // virtual function to NOT DO corrections in CPU but in GPU :
    virtual void ApplyGeometricCorrections( Visibilities& xcorr, MemoryBuffer<float>& w_cpu, MemoryBuffer<double>& frequencies);
   
    virtual void ApplyCableCorrections(Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies);


public :


   CPacerImagerHip();
   ~CPacerImagerHip();
};


#endif