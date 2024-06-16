#include "gpu_utils.h"

__global__ void mult_by_const( gpufftComplex *data, int size, double mult_value )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data[i].x = data[i].x*mult_value;
    data[i].y = data[i].y*mult_value;
}


__global__ void mult_by_const( float *data, int size, double mult_value )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data[i] = data[i]*mult_value;
}


__global__ void mult_arrays( float* data, float* data2, int size )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data[i] = data[i]*data2[i];
}

__global__ void mult_arrays( float* data, float* data2, float* data_out, int size )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data_out[i] = data[i]*data2[i];
}


__global__ void div_arrays( float* data, float* data2, int size )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data[i] = data[i] / data2[i];
}

__global__ void div_arrays( float* data, float* data2, float* data_out, int size )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    data_out[i] = data[i] / data2[i];
}

