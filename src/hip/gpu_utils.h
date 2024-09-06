#ifndef _GPU_UTILS_H__
#define _GPU_UTILS_H__

#include <gpu_fft.hpp>
#include <memory_buffer.hpp>

__global__ void mult_by_const( float *data, int size, double mult_value );

__global__ void mult_by_const( gpufftComplex *data, int size, double mult_value );

__global__ void mult_arrays( float* data, float* data2, int size );


void vector_sum_gpu( float* data_gpu, int image_size, int n_images, MemoryBuffer<float>& sum_gpu);

void fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, int n_images,  MemoryBuffer<float>& fnorm);

#endif
