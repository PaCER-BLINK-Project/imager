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
#define NTHREADS 1024
// TODO: fix this kernel to work for both NVIDIA and AMD
#define NWARPS 16


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


