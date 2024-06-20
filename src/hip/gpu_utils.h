#ifndef _GPU_UTILS_H__
#define _GPU_UTILS_H__

#include <gpu_fft.hpp>

__global__ void mult_by_const( float *data, int size, double mult_value );

__global__ void mult_by_const( gpufftComplex *data, int size, double mult_value );

__global__ void mult_arrays( float* data, float* data2, int size );

// two versions of array summing :
float sum_gpu_atomicadd( float* data_gpu, int size );
float sum_gpu_parallel_reduce( float* data_gpu, int size );

#endif
