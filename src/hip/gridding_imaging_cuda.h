#ifndef _GRIDDING_IMAGING_CUDA_H__
#define _GRIDDING_IMAGING_CUDA_H__

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>


#include <array_config_common.h>

// OLD kernel is still required to calculate UV_counter at least once for CONSTANT UV version, but can be removed in the future :
__global__ void gridding_imaging_cuda(int xySize, // size of the correlation matrix
                                      int n_ant,  // number of antennas 
                                      float *u_cuda, float *v_cuda, 
                                      int* antenna_flags, float* antenna_weights,
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      float *uv_grid_counter_cuda, float *uv_grid_real_cuda, float *uv_grid_imag_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda);

// same as above but using xcorr structure instead :
__global__ void gridding_imaging_cuda_xcorr( int xySize, // size of the correlation matrix
                                      int n_ant,
                                      float *u_cuda, float *v_cuda, 
                                      int* antenna_flags, float* antenna_weights,
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      VISIBILITY_TYPE *vis_cuda,  
                                      float *uv_grid_counter_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda);

//----------------------------------- Phase corrections : cable, geometric, apply calibration ----------------------------------- 
// kernel for applying geometric corrections:
__global__ void apply_geometric_corrections( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *w_cuda, double frequency_hz, double speed_of_light );

__global__ void apply_cable_corrections( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *cable_lengths_cuda, double frequency_hz, double speed_of_light );



// New optmised kernel :
//   - removing is_ODD ifs
//   - removing UV_real, UV_imag and UV_counter
__global__ void gridding_imaging_cuda_optimised(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, float *uv_grid_counter_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda);

// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_cuda_optimised_nocounter(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda);

// Counter calculation - just using threads (no BLOCKS) :
__global__ void calculate_counter( int xySize, // size of the correlation matrix
                                   float *u_cuda, float *v_cuda, 
                                   double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                   int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                   float *vis_real_cuda, float *vis_imag_cuda, 
                                   float *uv_grid_counter_cuda, double min_uv_cuda );


//---------------------------------------------------- Gridding kernels using BLOCKS with and without counter calculations  ----------------------------------------------------
// Gridding kernels using BLOCKS with and without counter calculations :
// same as gridding_imaging_cuda_blocks but no UV_real and UV_image separately only as m_in_buffer_cuda :
// and possibly more optimisations to come :
__global__ void gridding_imaging_cuda_blocks_optimised(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      float *uv_grid_counter_cuda_param, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda_param);

// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_cuda_blocks_optimised_nocounter(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda_param);
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

__device__ int calculate_pos(float u,
                             float v,
                             double delta_u_cuda,
                             double delta_v_cuda,
                             double wavelength_cuda,
                             double min_uv_cuda,
                             int n_pixels_cuda,
                             int center_x_cuda,
                             int center_y_cuda,
                             int is_odd_x_cuda,
                             int is_odd_y_cuda,
                             int uv_sign_cuda);

// test functions :
__global__ void vis2corrmatrix( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *vis_corrmatrix_re_cuda, float* vis_corrmatrix_im_cuda );

#endif 
