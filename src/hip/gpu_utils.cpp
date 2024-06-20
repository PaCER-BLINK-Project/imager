#include "gpu_utils.h"
#include "gpu_macros.hpp"
#include "pacer_imager_hip_defines.h"

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

// simple kernel using atomicAdd - it's actually the fastest on small arrays (like correlation matrix). Probably due to shared memory allocation overhead
__global__ void sum_atomic_add( float* data, int size, float* sum )
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= size ){
       return;
    }

    atomicAdd(sum,data[i]);
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
   sum_atomic_add<<<nBlocks,NTHREADS>>>( data_gpu, size, sum_gpu );
//   printf("DEBUG : using %d * %d = %d bytes of shared memory\n",nBlocks*NTHREADS,int(sizeof(float)),int(nBlocks*NTHREADS*sizeof(float)));
   gpuDeviceSynchronize();
   std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
   printf("BENCHMARK_chrono : sum_atomic_add took %.6f [ms]\n",time_span.count()*1000.00);

   gpuMemcpy((float*)(&sum_cpu), (float*)sum_gpu, sizeof(float),gpuMemcpyDeviceToHost);

   gpuFree( sum_gpu );

   return sum_cpu;
}

// shared memory version1
__global__ void sum_shared_mem( float* data, int size, float* sum )
{
    // copy from global GPU memory (data) to SHARED memory buffer shared_memory:
    // see https://developer.nvidia.com/blog/using-shared-memory-cuda-cc/
    extern __shared__ float shared_memory[];
    int tid = threadIdx.x;
    if (tid >= size ) return;
    shared_memory[tid] = data[tid];
     __syncthreads();
//    printf("DEBUG_VALUE = %.4f\n",shared_memory[tid]);

    // parallel reduction : ~/Desktop/GPU/documents/CUDA_programming/2024/MS/reduction-1.pdf
    // do reduction in shared mem
    for(unsigned int s=1; s < blockDim.x; s *= 2) {
       if(tid % (2*s) == 0){
          shared_memory[tid] += shared_memory[tid + s];
       }
       __syncthreads();
    }

    // write result for this block to global mem
    // if(tid == 0) g_odata[blockIdx.x] = sdata[0];
    if(tid == 0){
       sum[0] += shared_memory[0];
//       printf("DEBUG: partial sum = %.8f\n",shared_memory[0]); 
    }
}

__global__ void sum_short( float* data, int size, float* sum )
{
   int tid = threadIdx.x;
   if (tid >= size ) return;
   atomicAdd(sum,data[tid]);
}

float sum_gpu_parallel_reduce( float* data_gpu, int size )
{
    float* sum_gpu=NULL;
    float sum_cpu=0.00;
    
    // initialise sum_gpu variable;
    gpuMalloc((void **)&sum_gpu, sizeof(float));
    gpuMemcpy((float*)sum_gpu, (float*)(&sum_cpu), sizeof(float), gpuMemcpyHostToDevice); // initilise

// float ms;
// gpuEvent_t startEvent, stopEvent;
// checkCuda( gpuEventRecord(startEvent,0) );
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
    int n_parts = int(size/NTHREADS);
//    printf("Processing %d parts each %d elements (total = %d vs. size = %d -> rest = %d\n",n_parts,NTHREADS,(n_parts*NTHREADS),size,(size-(n_parts*NTHREADS)));
    for(int i=0;i<n_parts;i++){
       sum_shared_mem<<<1,NTHREADS,NTHREADS*sizeof(float)>>>( data_gpu+(i*NTHREADS), NTHREADS, sum_gpu );
       gpuDeviceSynchronize();
    }
    int rest = size - int(size/NTHREADS)*NTHREADS;
//    printf("DEBUG : unprocessed number of elements = %d\n",rest);
    if( rest >= 0 ){
       sum_short<<<1,rest>>>( data_gpu+(n_parts*NTHREADS), rest, sum_gpu );
    }
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
    printf("BENCHMARK_chrono : sum_gpu_parallel_reduce took %.6f [ms]\n",time_span.count()*1000.00);

//    checkCuda( gpuEventRecord(stopEvent,0) );
//    checkCuda( gpuEventSynchronize(stopEvent) );
//    checkCuda( gpuEventElapsedTime(&ms, startEvent, stopEvent) );
//    printf("BENCHMARK : execution of sum_atomic_add took %.12f [ms]\n", ms);

    gpuMemcpy((float*)(&sum_cpu), (float*)sum_gpu, sizeof(float),gpuMemcpyDeviceToHost);
//    printf("GPU Sum function = %.8f\n",sum_cpu);

//    checkCuda( gpuEventDestroy(startEvent) );
//    checkCuda( gpuEventDestroy(stopEvent) );
    gpuFree( sum_gpu );

    return sum_cpu;
}
