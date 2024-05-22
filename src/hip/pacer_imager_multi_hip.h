#ifndef _PACER_IMAGER_MULTI_H__
#define _PACER_IMAGER_MULTI_H__

#include <typeinfo>

#include "pacer_imager_hip.h"

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>


#include <vector>

// TODO :
//  - implement correlator input type, which will use different gridding kernels adjusted to the type of the correlator used


// Multi-frequency imager - here we assume that the blocks are for separate frequency channels.
// Hence, if UVW=CONST -> we can calculate N-fine-channel counters and UVW planes and keep them like this forever 
class CPacerImagerMultiFreqHip : public CPacerImagerHip
{
protected :
   gpufftComplex* m_in_cross_correlations_gpu; // input visibilities (Cross-correlations) pointer can be set externally to use memory buffer from a correlator or 
                                              // can be allocated here (may be in L-files ordering)
   float*        m_in_auto_correlations_gpu;  // auto-correlations currently not used, but can added later 
   
   int*          m_mapping_array_gpu; // equivalent of mapping matrix in array_config.h(cpp) : m_CorrelatorMatrix
   int*          m_mapping_array_cpu;
   InputMapping* m_inputs_gpu;             // equivalent of vector std::vector<InputMapping> m_Inputs in array config 

   std::vector<double> m_UVCounterSum; // sums of UV counters (normalisation factors) for all the channels 
   int m_nCrossCorrBlockSize; // not in bytes, but in count of gpufftComplex 
   int m_nAutoCorrBlockSize;  // not in bytes, but in count of floats

   // Clean GPU Memory 
   virtual void CleanGPUMemory(); 
   
   // check if memory has been allocated :
   bool CheckMemoryBuffers();
   
   // memset to zero UV grid etc :
   bool ZeroMemoryBuffers();
   
   // initialisation and cleaning of cuda/hip Streams :   
   void InitStreams();
   void CleanStreams();
   
   // Debug functions for saving intermediate control files :
   void SaveDebugFiles( std::vector< std::vector<int> >& antenna_mapping_array, std::vector<InputMapping>& inputs_cpu, int size );
   bool SaveGriddedVisibilities();
   
   // member variables: 
   int m_nStreams;    // number of streams used. Streams are like queues where kernel executions are submitted. Kernels are executed in the order of submission (call) and 
                      // if there are multiple streams they can be executed in parallel. Default 15  
   bool m_bGridBlocks; // if m_nUseBlocks>1 : use paralellisation in BLOCKS (grid of CUDA blocks of threads) rather than streams, Default = false (no blocks used, but 15 Streams are)
   
   int m_N; // i.e. currently frequency channels (fine channels) 
            // number of blocks in gridding and cuFFTPlanMany (not gridblocks!!!) in case of this class most likely frequency channels to be processed in the same time
   
   
   // Streams :
   gpuStream_t* m_GpuStreams;
   
public :
   CPacerImagerMultiFreqHip();
   ~CPacerImagerMultiFreqHip();
   
   // Initialisation of UVW and counter on GPU :
   bool InitialiseUVW();

   bool InitialiseAntennaMapping( std::vector< std::vector<int> >& antenna_mapping_array, std::vector<InputMapping>& inputs_cpu );

   // Allocate required memory 
   virtual void AllocGPUMemory( int corr_size, 
                                int image_size
                              );

   // allocate buffer for input visibilities, which can be used in gridding :
   void AllocInVisBuffer( int size_cross_corr , int size_auto_corr ); // size - number of gpufftComplex elements (cross-correlations)
   
   // Copy external CPU data to internal GPU buffers ( m_in_auto_correlations_gpu and m_in_cross_correlations_gpu )
   // WARNING : copied only exactly as many elements as allocated by AllocInVisBuffer and set in member variables m_nCrossCorrBlockSize and m_nAutoCorrBlockSize :
   void CopyInputVisibilities( float* pBaselinesAutoBufferGPU, gpufftComplex* pBaselinesCrossBuffer );

   
   inline void SetNStreams( int nStreams ){ m_nStreams = nStreams; }
   inline int  GetNStreams(){ return m_nStreams; }
   
   void SetGridBlocks( bool bGridBlocks );
   inline bool GetGridBlocks(){ return m_bGridBlocks; }
   
   void SetNBlocks( int N );
   inline int  GetNBlocks(){ return m_N; }

   inline void SetNFreqChannels( int N ){ m_N = N; }
   inline int  GetNFreqChannels(){ return m_N; }

   // accessing output data products (images etc) :
   int CopyImagesGpu2Cpu();
   bool GetShiftedImage( int freq_channel, CBgFits& real_image_out, CBgFits& imag_image_out );
   
   // Saving FITS files on request :
   bool SaveSkyImage( int ch, int timestamp, double unixtime, int image_size, const char* szOutDir, int iSaveImaginaryFITS=0 );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 
   // GOAL   : function calling both gridding_fast and imaging. It's a function using visibilities in the form out correlator output 
   //          (ordering can be decided by a to-be added parameter CorrelatorType
   // 
   // INPUT  : 
   //          void* pBaselinesAutoBuffer, void* pBaselinesCrossBuffer : visibilities (REAL and IMAG 2D arrays as FITS class) in CorrelatorType=Randall : float* and gpufftComplex* 
   //          n_fine_channels - number of fine channels in the input visibilities, bw_mhz - bandwidth of input data,  
   //          delta_u, delta_v : size of the UV cell 
   //          frequency_mhz : frequency in MHz
   //
   // OUTPUT : the paradigm behind this function is that it does not save any files and can be a part of a bigger optimised pipeline storing all intermediate data products in GPU memory 
   //          - uv_grid_real, uv_grid_imag : visibilities on UV grid (real and imag arrays)
   //          - uv_grid_counter : visibility counter and 
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual void gridding_imaging_multi_freq( int n_fine_channels, double bw_mhz, // description on input buffers (output of the correlator)
                                             // CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, - use m_U, m_V, m_W 
                                             const char* weighting="", // weighting : U for uniform (others not implemented)
                                             bool do_gridding=true,    // if gridding is performed 
                                             bool do_dirty_image=true  // if cu/hip FFT is called 
                                           );

   
};

#endif
