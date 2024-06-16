#ifndef _GPU_UTILS_H__
#define _GPU_UTILS_H__

#include <gpu_fft.hpp>

__global__ void mult_by_const( float *data, int size, double mult_value );

__global__ void mult_by_const( gpufftComplex *data, int size, double mult_value );

__global__ void mult_arrays( float* data, float* data2, int size );

#endif
