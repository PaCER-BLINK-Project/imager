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


//----------------------------------------------------- FFT SHIFT ----------------------------------------------------- 
// FFT shift in X direction:
__global__ void fft_shift_x( float* data, int xSize, int size, int y )
{
   extern __shared__ float shared_memory[];
   int t = threadIdx.x;
   int center_freq_x = xSize / 2;
   int is_odd = 0;

   // initialise shared memory 
   // for(; t < xSize; t += blockDim.x) shared_memory[t] = data[t];
   for( int x=t; x < xSize; x += blockDim.x){
      if( x < center_freq_x ){
         shared_memory[center_freq_x+x] = data[x];
      }else{
         shared_memory[x-(center_freq_x+is_odd)] = data[x];
      }
   }
   __syncthreads();

  for( int x=t; x < xSize; x += blockDim.x){
    data[x] = shared_memory[x]; // back to memory (in place)
  }
  __syncthreads();
}

// FFT shift in Y direction (may be inefficient due to long strides)
__global__ void fft_shift_y( float* data, int xSize, int ySize, int size, int x )
{
   extern __shared__ float shared_memory[];
   int t = threadIdx.x;
   int center_freq_y = ySize / 2;
   int is_odd = 0;

   // 
   for( int y=t; y < ySize; y += blockDim.x){
      int pos = y*xSize + x; // calculate index in the entire array

      if( y < center_freq_y ){
         shared_memory[center_freq_y+y] = data[pos];
      }else{
         shared_memory[y-(center_freq_y+is_odd)] = data[pos];
      }
   }
   __syncthreads();

  for( int y=t; y < ySize; y += blockDim.x){
    int pos = y*xSize + x; // calculate index in the entire array
    data[pos] = shared_memory[y]; // back to memory (in place)
  }
  __syncthreads();
}


bool fft_shift_gpu( float* data_gpu, int xSize, int ySize )
{
   if( xSize != ySize ){
      printf("ERROR : function fft_shift_gpu currently only handles square images %d != %d\n",xSize,ySize);
      return false;
   }
   if( (xSize % 2) != 0 ){
      printf("ERROR : function fft_shift_gpu currently only handles even image size -> %d not supported\n",xSize);
      return false;
   }

   int size = xSize*ySize;
   int shared_memory_size_bytes = xSize*sizeof(float);
   if( shared_memory_size_bytes >= MAX_SHARED_MEMORY_BYTES ){
      printf("ERROR : maximum size of shared memory exceeded %d > %d\n",shared_memory_size_bytes,MAX_SHARED_MEMORY_BYTES);
   }

   for(int y=0;y<ySize;y++){
      fft_shift_x<<<1,NTHREADS,shared_memory_size_bytes>>>(data_gpu + y*xSize, xSize, size, y );
   }

   for(int x=0;x<xSize;x++){
      fft_shift_y<<<1,NTHREADS,shared_memory_size_bytes>>>(data_gpu, xSize, ySize, size, x );
   }

   return true;
}


//-------------------------------------------------------------------- FFT shift on complex data ---------------------------------------------------------------------------------------------------
// FFT shift in X direction:
__global__ void fft_shift_and_norm_x( gpufftComplex* data, int xSize, int size, int y )
{
   extern __shared__ gpufftComplex shared_memory_complex[];
   int t = threadIdx.x;
   int center_freq_x = xSize / 2;
   int is_odd = 0;

   // initialise shared memory 
   // for(; t < xSize; t += blockDim.x) shared_memory[t] = data[t];
   for( int x=t; x < xSize; x += blockDim.x){
      if( x < center_freq_x ){
         shared_memory_complex[center_freq_x+x] = data[x];
      }else{
         shared_memory_complex[x-(center_freq_x+is_odd)] = data[x];
      }
   }
   __syncthreads();

  for( int x=t; x < xSize; x += blockDim.x){
    data[x] = shared_memory_complex[x]; // back to memory (in place)
  }
  __syncthreads();
}

// FFT shift in Y direction (may be inefficient due to long strides)
__global__ void fft_shift_and_norm_y( gpufftComplex* data, int xSize, int ySize, int size, int x, float fnorm )
{
   extern __shared__ gpufftComplex shared_memory_complex[];
   int t = threadIdx.x;
   int center_freq_y = ySize / 2;
   int is_odd = 0;

   // 
   for( int y=t; y < ySize; y += blockDim.x){
      int pos = y*xSize + x; // calculate index in the entire array

      if( y < center_freq_y ){
         shared_memory_complex[center_freq_y+y] = data[pos];
      }else{         
         shared_memory_complex[y-(center_freq_y+is_odd)] = data[pos];
      }
   }
   __syncthreads();

  for( int y=t; y < ySize; y += blockDim.x){
    int pos = y*xSize + x; // calculate index in the entire array
    data[pos].x = shared_memory_complex[y].x*fnorm; // back to memory (in place)
    data[pos].y = shared_memory_complex[y].y*fnorm; // back to memory (in place)
  }
  __syncthreads();


}



// fft shift on complex data :
bool fft_shift_and_norm_gpu( gpufftComplex* data_gpu, int xSize, int ySize, float fnorm /*=1.00*/ )
{
   if( (xSize % 2) != 0 ){
      printf("ERROR : function fft_shift_gpu currently only handles even image size -> %d not supported\n",xSize);
      return false;
   }

   int size = xSize*ySize;
   int shared_memory_size_bytes = xSize*sizeof(gpufftComplex);
   if( shared_memory_size_bytes >= MAX_SHARED_MEMORY_BYTES ){
      printf("ERROR : maximum size of shared memory exceeded %d > %d\n",shared_memory_size_bytes,MAX_SHARED_MEMORY_BYTES);
   }

   for(int y=0;y<ySize;y++){
      fft_shift_and_norm_x<<<1,NTHREADS,shared_memory_size_bytes>>>(data_gpu + y*xSize, xSize, size, y );
   }

   shared_memory_size_bytes = ySize*sizeof(gpufftComplex);
   if( shared_memory_size_bytes >= MAX_SHARED_MEMORY_BYTES ){
      printf("ERROR : maximum size of shared memory exceeded %d > %d\n",shared_memory_size_bytes,MAX_SHARED_MEMORY_BYTES);
   }
   for(int x=0;x<xSize;x++){
      fft_shift_and_norm_y<<<1,NTHREADS,shared_memory_size_bytes>>>(data_gpu, xSize, ySize, size, x, fnorm ); // final multiplication
   }

   return true;
}


