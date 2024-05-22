#ifndef _GRIDDING_MULTI_IMAGE_CUDA_H__
#define _GRIDDING_MULTI_IMAGE_CUDA_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>


#include <array_config_common.h>

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//
//
// Functions for multi-block version working on L-file like input data:
//
//
// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_lfile_vis( int xySize,                   // size of correlation matrix - should really be half of the corr-matrix without diagonal
                                            int channel,                  // frequency channel
                                            double wavelength_cuda,       // wavelength [m]
                                            int pol1, char cPol1, int pol2, char cPol2, // which correlation products to grid / image , cPol1 and 2 are to avoid ifs
                                            float *u_cuda, float *v_cuda, // UVW
                                            double delta_u_cuda, double delta_v_cuda, 
//                                            double min_uv_cuda, // MinUV - will try to avoid if-s here TODO : remove other baselines by 0-weighting visibilities in these cells (matrix multiplication -> Tensor Cores ?)
                                            gpufftComplex* in_visibilities_corr, int in_vis_corr_size, // CROSS-CORRELATIONS 
                                            float* in_visibilities_auto, int in_vis_auto_size,        // AUTO-CORRELATIONS 
                                            gpufftComplex* out_visibilities_gridded, int image_size_cuda, // full image size (usually n_pixels_cuda*n_pixels_cuda)
                                            int n_pixels_cuda, // side of image (number of pixels is n_pixels_cuda*n_pixels_cuda)
                                            int n_ant, // number of antennas 
                                            int n_channels,
                                            InputMapping* inputs, int n_inputs, 
                                            int* mapping_array
                                          );

// BLOCKS version of the above kernel gridding_imaging_lfile_vis - channels are implemented as separate blocks -> wave length has to be relculated inside 
// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_lfile_vis_blocks( int xySize,                   // size of correlation matrix - should really be half of the corr-matrix without diagonal
                                                   float first_channel_center_freq_mhz, // center frequency of the 1st fine channel 
                                                   float channel_bw_mhz,                // bandwidth of fine channels 
                                            int pol1, char cPol1, int pol2, char cPol2, // which correlation products to grid / image, cPol1 and 2 are to avoid ifs                                                   
                                            float *u_cuda, float *v_cuda, // UVW
                                            double delta_u_cuda, double delta_v_cuda, 
//                                            double min_uv_cuda, // MinUV - will try to avoid if-s here TODO : remove other baselines by 0-weighting visibilities in these cells (matrix multiplication -> Tensor Cores ?)
                                            gpufftComplex* in_visibilities_corr, int in_vis_corr_size, // CROSS-CORRELATIONS
                                            float* in_visibilities_auto, int in_vis_auto_size,        // AUTO-CORRELATIONS
                                            gpufftComplex* out_visibilities_gridded_param, int image_size_cuda, // full image size (usually n_pixels_cuda*n_pixels_cuda)
                                            int n_pixels_cuda, // side of image (number of pixels is n_pixels_cuda*n_pixels_cuda)
                                            int n_ant, // number of antennas 
                                            int n_channels,
                                            InputMapping* inputs, int n_inputs, 
                                            int* mapping_array
                                          );



// Cuda kernal: gridding and cuFFT 
// same as calculate_counter but does not require visibilities (parameters vis_real_cuda and vis_imag_cuda removed) :
__global__ void calculate_counter_novis( int xySize, // size of the correlation matrix
                                   float *u_cuda, float *v_cuda, 
                                   double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                   int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                   float *uv_grid_counter_cuda, double min_uv_cuda );


#endif 
