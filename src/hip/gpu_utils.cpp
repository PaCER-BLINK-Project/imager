#include "gpu_utils.h"
#include "gpu_macros.hpp"
#include "pacer_imager_hip_defines.h"
#include <exception>
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


// SUM of array:
#define NTHREADS 1024
#define NWARPS (NTHREADS / warpSize)


__global__ void vector_sum(float *values, unsigned int nitems, float* result){
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    __shared__ float partial_sums[NWARPS];
    unsigned int warpId = threadIdx.x / warpSize;
    unsigned int laneId = threadIdx.x % warpSize; 
    unsigned int gridSize = gridDim.x * blockDim.x;
    unsigned int nloops = (nitems + gridSize  - 1) / gridSize;
    float myvalue = 0;
    for(unsigned int l = 0; l < nloops; l++, idx+=gridSize){
        if(idx < nitems) myvalue += values[idx]; 
    }
 
    for(unsigned int i = warpSize/2; i >= 1; i /= 2){
      float up = __shfl_down(myvalue, i); 
        if(laneId < i){
            myvalue += up; 
        }
    }
    if(laneId == 0 && warpId > 0) partial_sums[warpId] = myvalue;

    __syncthreads();
    // step 2
    if(warpId == 0){
        if(laneId > 0 && laneId < NWARPS) myvalue = partial_sums[laneId];
        for(unsigned int i = NWARPS/2; i >= 1; i >>= 1){
            float up = __shfl_down(myvalue, i, NWARPS); 
            if(laneId < i){
               myvalue += up; 
            }
        }
        if(laneId == 0) atomicAdd(result, myvalue);
    }  
}

float sum_gpu_atomicadd( float* data_gpu, int size )
{
   int nBlocks = (size + NTHREADS -1)/NTHREADS;
   float* sum_gpu=NULL;
   float sum_cpu=0.00;
   gpuMalloc((void **)&sum_gpu, sizeof(float));
   gpuMemcpy((float*)sum_gpu, (float*)(&sum_cpu), sizeof(float), gpuMemcpyHostToDevice); // initilise

//   gpuEvent_t startEvent, stopEvent; 
//   float ms;
   std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
//   checkCuda( gpuEventRecord(startEvent,0) );
   vector_sum<<<nBlocks,NTHREADS>>>( data_gpu, size, sum_gpu );
//   printf("DEBUG : using %d * %d = %d bytes of shared memory\n",nBlocks*NTHREADS,int(sizeof(float)),int(nBlocks*NTHREADS*sizeof(float)));
   gpuDeviceSynchronize();
   std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
   printf("BENCHMARK_chrono : sum_atomic_add took %.6f [ms]\n",time_span.count()*1000.00);

   gpuMemcpy((float*)(&sum_cpu), (float*)sum_gpu, sizeof(float),gpuMemcpyDeviceToHost);

   gpuFree( sum_gpu );

   return sum_cpu;
}

__device__ inline int calc_fft_shift(int pos, int side){
   int is_odd = side % 2;
   return (pos + side/2 + is_odd) % (side);
}

//-------------------------------------------------------------------- FFT shift on complex data ---------------------------------------------------------------------------------------------------
// FFT shift in X direction:
__global__ void fft_shift_and_norm_x(gpufftComplex* data, size_t image_x_side, size_t image_y_side){
   size_t tid {blockDim.x * blockIdx.x + threadIdx.x};
   if(tid >= image_x_side / 2 * image_y_side) return;

   size_t src_row = tid / (image_x_side / 2);
   size_t src_col = tid % (image_x_side / 2);
   size_t src = src_row * image_x_side + src_col;
   size_t dst_col = calc_fft_shift(src_col, image_x_side);
   size_t dst = src_row * image_x_side + dst_col;

   gpufftComplex tmp = data[dst];
   data[dst] = data[src];
   data[src] = tmp;
}



__global__ void fft_shift_and_norm_y(gpufftComplex* data, size_t image_x_side, size_t image_y_side, float fnorm ){
   size_t tid {blockDim.x * blockIdx.x + threadIdx.x};
   if(tid >= image_x_side * (image_y_side / 2)) return;

   size_t src_row = tid / image_x_side;
   size_t src_col = tid % image_x_side;
   size_t src = src_row * image_x_side + src_col;
   size_t dst_row = calc_fft_shift(src_row, image_y_side);
   size_t dst = dst_row * image_x_side + src_col;

   gpufftComplex tmp = data[dst] * fnorm;
   data[dst] = data[src] * fnorm;
   data[src] = tmp;
}



// fft shift on complex data :
void fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, float fnorm /*=1.00*/ ){
   if( (xSize % 2) != 0 ){
      throw std::invalid_argument {"function fft_shift_gpu currently only handles even image size"};
   }

   int size = xSize*ySize;
   int n_threads_needed {(xSize / 2) * ySize};
   int n_blocks {(n_threads_needed + NTHREADS - 1) / NTHREADS};
   fft_shift_and_norm_x<<<n_blocks, NTHREADS>>>(data_gpu, xSize, ySize);
   n_threads_needed = xSize * (ySize / 2);
   n_blocks = (n_threads_needed + NTHREADS - 1) / NTHREADS;
   fft_shift_and_norm_y<<<n_blocks, NTHREADS>>>(data_gpu, xSize, ySize, fnorm);
}


