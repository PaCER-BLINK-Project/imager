#include "gpu_utils.h"
#include <gpu_macros.hpp>
#include "pacer_imager_hip_defines.h"
#include <exception>
#include <memory_buffer.hpp>

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


__global__ void vector_sum_kernel(float *values, unsigned int nitems, int n_images, float* result){
    __shared__ float partial_sums[NWARPS];
    unsigned int warpId = threadIdx.x / warpSize;
    unsigned int laneId = threadIdx.x % warpSize; 
    unsigned int gridSize = gridDim.x * blockDim.x;

    for(int img_id = 0; img_id < n_images; img_id++){
      unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
      float myvalue = 0;
      for(; idx < nitems; idx+=gridSize){
         myvalue += values[img_id * nitems + idx]; 
      }
   
      for(unsigned int i = warpSize/2; i >= 1; i /= 2){
         float up = __shfl_down(myvalue, i); 
         if(laneId < i){
               myvalue += up; 
         }
         #ifdef __NVCC__
         __syncwarp();
         #endif
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
               #ifdef __NVCC__
               __syncwarp();
               #endif
         }
         if(laneId == 0) atomicAdd(result + img_id, myvalue);
      }
    }  
}


void vector_sum_gpu( float* data_gpu, int image_size, int n_images, MemoryBuffer<float>& sum_gpu){
   struct gpuDeviceProp_t props;
   int gpu_id = -1;
   gpuGetDevice(&gpu_id);
   gpuGetDeviceProperties(&props, gpu_id);
   unsigned int n_blocks = props.multiProcessorCount * 2;
   sum_gpu.to_gpu();
   gpuMemset(sum_gpu.data(), 0, n_images * sizeof(float));
   vector_sum_kernel<<<n_blocks, NTHREADS>>>( data_gpu, image_size, n_images, sum_gpu.data());
}

__device__ inline int calc_fft_shift(int pos, int side){
   int is_odd = side % 2;
   return (pos + side/2 + is_odd) % (side);
}

//-------------------------------------------------------------------- FFT shift on complex data ---------------------------------------------------------------------------------------------------
// FFT shift in X direction:
__global__ void fft_shift_and_norm_x(gpufftReal* data, size_t image_x_side, size_t image_y_side, int n_images ){
   size_t tid {blockDim.x * blockIdx.x + threadIdx.x};
   if(tid >= image_x_side / 2 * image_y_side) return;

   size_t src_row = tid / (image_x_side / 2);
   size_t src_col = tid % (image_x_side / 2);
   size_t src = src_row * image_x_side + src_col;
   size_t dst_col = calc_fft_shift(src_col, image_x_side);
   size_t dst = src_row * image_x_side + dst_col;
   for(int img_id = 0; img_id < n_images; img_id++){
      gpufftReal tmp = data[img_id * image_x_side * image_y_side + dst];
      data[img_id * image_x_side * image_y_side + dst] = data[img_id * image_x_side * image_y_side + src];
      data[img_id * image_x_side * image_y_side + src] = tmp;
   }
}



__global__ void fft_shift_and_norm_y(gpufftReal* data, size_t image_x_side, size_t image_y_side, int n_images, float *fnorm ){
   size_t tid {blockDim.x * blockIdx.x + threadIdx.x};
   if(tid >= image_x_side * (image_y_side / 2)) return;

   size_t src_row = tid / image_x_side;
   size_t src_col = tid % image_x_side;
   size_t src = src_row * image_x_side + src_col;
   size_t dst_row = calc_fft_shift(src_row, image_y_side);
   size_t dst = dst_row * image_x_side + src_col;

   for(int img_id = 0; img_id < n_images; img_id++){
      gpufftReal tmp = data[img_id * image_x_side * image_y_side + dst] / fnorm[img_id];
      data[img_id * image_x_side * image_y_side + dst] = data[img_id * image_x_side * image_y_side + src] / fnorm[img_id];
      data[img_id * image_x_side * image_y_side + src] = tmp;
   }
}



// fft shift on complex data :
void fft_shift_and_norm_gpu( gpufftReal* data_gpu, int xSize, int ySize, int n_images, MemoryBuffer<float>& fnorm){
   if( (xSize % 2) != 0 ){
      throw std::invalid_argument {"function fft_shift_gpu currently only handles even image size"};
   }
   int size = xSize*ySize;
   int n_threads_needed {(xSize / 2) * ySize};
   int n_blocks {(n_threads_needed + NTHREADS - 1) / NTHREADS};
   fft_shift_and_norm_x<<<n_blocks, NTHREADS>>>(data_gpu, xSize, ySize, n_images);
   n_threads_needed = xSize * (ySize / 2);
   n_blocks = (n_threads_needed + NTHREADS - 1) / NTHREADS;
   fft_shift_and_norm_y<<<n_blocks, NTHREADS>>>(data_gpu, xSize, ySize, n_images, fnorm.data());
}


__global__ void averaging_kernel(const float *data, size_t n_pixels, size_t n_images, float *out){
   size_t i {blockIdx.x * blockDim.x + threadIdx.x};
   if(i >= n_pixels) return;
   float tmp = data[i];
   for(size_t t {i + n_pixels}; t < n_pixels * n_images; t += n_pixels) tmp += data[t];
   tmp  /= n_images;
   out[i] = tmp;
}


 Images image_averaging_gpu(const Images& images){
   size_t image_size = images.image_size();
   size_t n_images = images.integration_intervals() * images.nFrequencies;
   MemoryBuffer<float> avg_image {image_size, true};
   unsigned int n_blocks {static_cast<unsigned int>((image_size + NTHREADS - 1) / NTHREADS)};
   averaging_kernel<<<n_blocks, NTHREADS>>>(images.data(), image_size, n_images, avg_image.data());
   gpuCheckLastError();
   return {std::move(avg_image), images.obsInfo, images.obsInfo.nTimesteps, images.obsInfo.nFrequencies, images.side_size};
 }