// In order to include the gridding_gpu.cu 
/*#include <cuda.h> 
#include <cuda_runtime.h>
// So that it recognises: blockIdx
#include <device_launch_parameters.h>
#include <cufftw.h>
#include <cufft.h>*/

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>

// msfitslib :
#include <array_config_common.h>

// local defines :
#include "pacer_imager_hip_defines.h"

// TODO : can this duplication of code be fixed ?
__device__ int calculate_pos_local(float u,
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
                             int uv_sign_cuda) // will be called with +1 and -1 
{
   // Operation 1: uv_lambda()
   double u_lambda = (u)/(wavelength_cuda); 
   double v_lambda = (v)/(wavelength_cuda); 

   // Calculating distance between the two antennas 
   double uv_distance = sqrt(u_lambda*u_lambda + v_lambda*v_lambda);

   if( uv_distance > min_uv_cuda )
   {            
      // (For all the rows of the Correlation Matrix)
      // Operation 2: uv_index()
      double u_pix, v_pix;
      u_pix = round(u_lambda/delta_u_cuda); 
      v_pix = round(v_lambda/delta_v_cuda);
      int u_index = uv_sign_cuda*u_pix + (n_pixels_cuda/2); 
      int v_index = uv_sign_cuda*v_pix + (n_pixels_cuda/2);

      // Operation 3: x_grid, y_grid (With FFT-UNSHIFT)
      int x_grid = 0;
      int y_grid = 0; 

      if( u_index < center_x_cuda )
         x_grid = u_index + center_x_cuda + is_odd_x_cuda;
      else
         x_grid = u_index - center_x_cuda;

      if( v_index < center_y_cuda )
         y_grid = v_index + center_y_cuda + is_odd_y_cuda;
      else
         y_grid = v_index - center_y_cuda; 
     
      // TODO : understand why this is giving wrong image with a black stripe in the centre !
      // This may optimise this code in the future (remove if-s) if it also produces the same results
      // WARNING : does not support ODD image sizes (add ASSERT !!!)
      // TODO : maybe copysign(int,int) is required - doesn't copysign use if-s too ?
      // int x_grid = round(u_index + copysignf( 1.0, float(center_x_cuda-u_index))*center_x_cuda);
      // int y_grid = round(v_index + copysignf( 1.0, float(center_y_cuda-v_index))*center_y_cuda);
      

      // Operation 4: Assignment of (re,im)vis to uv_grid
      // Position for assignment 
      return (n_pixels_cuda*y_grid) + x_grid; 
   }

   // same as else :   
   return -1;
}


//
//
// Functions for multi-block version working on L-file like input data:
//
//

// Cuda kernal: gridding and cuFFT 
// same as calculate_counter but does not require visibilities (parameters vis_real_cuda and vis_imag_cuda removed) :
__global__ void calculate_counter_novis( int xySize, // size of the correlation matrix
                                   float *u_cuda, float *v_cuda, 
                                   double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                   int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                   float *uv_grid_counter_cuda, double min_uv_cuda )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    int pos = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
    if(pos>=0 && pos<image_size_cuda)
    {
       // Allocating in uv_grid
       atomicAdd(&uv_grid_counter_cuda[pos],1);
    }   

    int pos2 = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
    if(pos2>=0 && pos2<image_size_cuda)
    {
       atomicAdd(&uv_grid_counter_cuda[pos2],1);
    }
}

// equivalent of CLfileReader::GetCrossCorrData
__device__ gpufftComplex*  get_cross_corr_data( int input1, int input2, int n_inputs, int n_channels, int* mapping_array, gpufftComplex* in_visibilities_cross_corr )
{
   if(!mapping_array){
      return NULL;
   }  
 
   // TODO : can this code, and the code in CLfileReader::GetCurrentCorrMatrix be replaced by some memcpy :
   // TODO : check what the matrix m_ArrayConfig.m_CorrelatorMatrix contains - print it somewhere !
   int index = input2*n_inputs + input1;
   int baseline_index = mapping_array[index];
   if ( baseline_index < 0 ){
       // TODO : make sure this if is elimated !!!
       // printf("ERROR : this should never happen in this code baseline_index = %d for inputs %d and %d, n_inputs = %d, n_channels = %d (index = %d)\n",baseline_index,input1,input2,n_inputs,n_channels,index);
       // baseline_index = m_ArrayConfig.m_CorrelatorMatrix[input1][input2];
       baseline_index = mapping_array[input1*n_inputs + input2];
       if( baseline_index < 0 ){
          printf("ERROR : this should never happen in this code baseline_index = %d for inputs %d and %d, n_inputs = %d, n_channels = %d (index = %d)\n",baseline_index,input1,input2,n_inputs,n_channels,index);
       }
    }

//   printf("DEBUG : get_cross_corr_data( %d , %d , %d , %d ) -> index = %d, baseline_index = %d\n",input1,input2,n_inputs,n_channels,index,baseline_index);
   gpufftComplex* out_data = in_visibilities_cross_corr + baseline_index*n_channels;

   return out_data;
}

__device__ InputMapping* get_input_info( int antenna ,int pol, char cPol, InputMapping* inputs, int n_inputs )
{
   if( n_inputs <= 0 ){
      return NULL;
   }
 
   int inputIndex = antenna*2 + pol;
   if( inputIndex < n_inputs ){
      InputMapping& input = inputs[ inputIndex ];
      
      if( input.antenna == antenna && input.pol == cPol ){
         return &( inputs[ inputIndex ] );
      }else{
         printf("ERROR : currently only very simple mapping is handled. Requested antenna/pol = %d/%c vs. found at index=%d is %d/%c\n",antenna,cPol,inputIndex,input.antenna,input.pol);
/*         static bool warnined_printed = false;
         if( !warnined_printed ){
            printf("WARNING : using experimental version of mapping ...\n");
            warnined_printed = true;
         }

         for(int i=0;i<n_inputs;i++){
            InputMapping& input = inputs[ i ];

            if( input.antenna == antenna && input.pol == cPol ){
               return (&input);
            }
         }*/
      }
   }

   printf("ERROR: could not find information on antenna %d and polarisation = %d -> cannot continue with this error as results may be wrong !!!\n",antenna,pol);
   return NULL;
}


__device__ gpufftComplex*  get_cross_corr_data_antpol( int antenna1, int pol1, char cPol1, int antenna2, int pol2, char cPol2, 
                                                      InputMapping* inputs, int n_inputs, 
                                                      int n_channels, int* mapping_array, gpufftComplex* in_visibilities_cross_corr )
{
   if(!mapping_array){
      return NULL;
   }  
 
   InputMapping* pInput1 = get_input_info( antenna1 , pol1, cPol1, inputs, n_inputs );
   InputMapping* pInput2 = get_input_info( antenna2 , pol2, cPol2, inputs, n_inputs );

   if( pInput1 && pInput2 ){
      int input1 = pInput1->input;
      int input2 = pInput2->input;

      if( input1 == input2 ){
         printf("ERROR in code : cannot ask for the same inputs %d in GetCrossCorrData (ant1=%d:%d ant2=%d:%d)\n",input1,antenna1,pol1,antenna2,pol2);
         return NULL;
      }

      return get_cross_corr_data( input1, input2, n_inputs, n_channels, mapping_array, in_visibilities_cross_corr );
   }

   printf("ERROR : could not find on of the inputs for ant,pol pair (%d,%d) or (%d,%d)\n",antenna1,pol1,antenna2,pol2);
   return NULL;
}

// equivalent of float* CLfileReader::GetAutoCorrData( int input )
__device__ float* get_auto_corr_data( int input, int n_inputs, int n_channels, int* mapping_array, float* in_visibilities_auto_corr )
{
   if(!in_visibilities_auto_corr){
       return NULL;
   }

   float* out_data = in_visibilities_auto_corr + input*n_channels;
   return out_data;
}

__device__ float*  get_auto_corr_data_antpol( int antenna1, int pol1, char cPol1,  
                                              InputMapping* inputs, int n_inputs, 
                                              int n_channels, int* mapping_array, float* in_visibilities_auto_corr )
{
   if(!mapping_array){
      return NULL;
   }  
 
   InputMapping* pInput1 = get_input_info( antenna1 , pol1, cPol1, inputs, n_inputs );

   if( pInput1 ){
      int input1 = pInput1->input;

      return get_auto_corr_data( input1, n_inputs, n_channels, mapping_array, in_visibilities_auto_corr );
   }

   printf("ERROR in get_auto_corr_data_antpol : could not find on of the input for ant,pol = (%d,%d)\n",antenna1,pol1);
   return NULL;
}

// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__device__ void gridding_imaging_lfile_vis_base( int i, 
                                            int xySize,                   // size of correlation matrix - should really be half of the corr-matrix without diagonal
                                            int freq_channel,             // frequency channel
                                            double wavelength_cuda,       // wavelength [m]
                                            int pol1, char cPol1, int pol2, char cPol2, // which correlation products to grid / image , cPol1 and 2 are to avoid ifs
                                            float *u_cuda, float *v_cuda, // UVW
                                            double delta_u_cuda, double delta_v_cuda, 
//                                            double min_uv_cuda, // MinUV - will try to avoid if-s here TODO : remove other baselines by 0-weighting visibilities in these cells (matrix multiplication -> Tensor Cores ?)
                                            gpufftComplex* in_visibilities_corr, int in_vis_corr_size,
                                            float* in_visibilities_auto, int in_vis_auto_size,
                                            gpufftComplex* out_visibilities_gridded, 
                                            int image_size_cuda, // full image size (usually n_pixels_cuda*n_pixels_cuda)
                                            int n_pixels_cuda,
                                            int n_ant, // number of antennas 
                                            int n_channels,
                                            InputMapping* inputs, int n_inputs, 
                                            int* mapping_array
                                          )
{   
    // hardcoded in this version:
    float min_uv_cuda = -1000;
    int center_x_cuda = int(n_pixels_cuda/2);
    int center_y_cuda = int(n_pixels_cuda/2);
    int is_odd_x_cuda = (n_pixels_cuda % 2);
    int is_odd_y_cuda = is_odd_x_cuda;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
// TODO : add calculation of position in the visibilies in the same way as it would be done in a correletion matrix : see functions in lfile_reader
 

    int ant1 = (i % n_ant);
    int ant2 = (i / n_ant);

    if( ant1 > ant2 ){ // ant1 > ant2 replicated CPU code in CPacerImager::gridding_fast to only use half of the correlation matrix !
        gpufftComplex* baseline_data = get_cross_corr_data_antpol( ant1, pol1, cPol1, ant2, pol2, cPol2, inputs, n_inputs, n_channels, mapping_array, in_visibilities_corr  );

       if( baseline_data ){ // ant1 > ant2 replicated CPU code in CPacerImager::gridding_fast to only use half of the correlation matrix !
          float re = baseline_data[freq_channel].x;
          float im = baseline_data[freq_channel].y;

          // TEST :
          // if( freq_channel == 1 ){
          //    printf("DEBUG : i = %d -> ant1 = %d, ant2 = %d, channel = %d -> re/im = %.6f / %.6f\n",i,ant1,ant2,freq_channel,re,im);
          // }
    
          int pos = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
          if(pos>=0 && pos<image_size_cuda)
          {
             atomicAdd(&out_visibilities_gridded[pos].x,re);
             atomicAdd(&out_visibilities_gridded[pos].y,im);
          }   

          // DEBUG :
          // leaving these temporary commted out - just in case I need to get back to this version:
          /*int x = pos % n_pixels_cuda;
          int y = pos / n_pixels_cuda;
          if( x == 5 && y == 1 && freq_channel>=0 ){
             printf("DEBUG_gridding_imaging_lfile_vis [1] (ch=%d) : (%d,%d) added %.8f %.8f -> %.8f %.8f , (u,v) = (%.4f,%.4f) , pos = %d ants %d-%d, i=%d\n",freq_channel,x,y,re,im,out_visibilities_gridded[pos].x,out_visibilities_gridded[pos].y,u_cuda[i]/wavelength_cuda, v_cuda[i]/wavelength_cuda,pos,ant1,ant2,i);
          }*/

          int pos2 = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
          if(pos2>=0 && pos2<image_size_cuda)
          {
             // Allocating inside m_in_buffer as well 
             atomicAdd(&out_visibilities_gridded[pos2].x,re);
             atomicAdd(&out_visibilities_gridded[pos2].y,-im);
          }

          // DEBUG :
          // leaving these temporary commted out - just in case I need to get back to this version:
          /*x = pos2 % n_pixels_cuda;
          y = pos2 / n_pixels_cuda;
          if( x == 5 && y == 1 && freq_channel==0 ){
             printf("DEBUG_gridding_imaging_lfile_vis [2] (ch=%d) : (%d,%d) added2 %.8f %.8f -> %.8f %.8f , (u,v) = (%.4f,%.4f) , pos2=%d ants %d-%d , i=%d\n",freq_channel,x,y,re,-im,out_visibilities_gridded[pos2].x,out_visibilities_gridded[pos2].y,u_cuda[i]/wavelength_cuda, v_cuda[i]/wavelength_cuda,pos2,ant1,ant2,i);
          }*/
       }
    }else{
       if( ant1 == ant2 ){
          // TODO : ASSESS AND DECIDE if it should become a separate kernel to avoid to many conditional statements (if-s) ?
          float* baseline_data = get_auto_corr_data_antpol( ant1, pol1, cPol1, inputs, n_inputs, n_channels, mapping_array, in_visibilities_auto  ); 

          if( baseline_data ){ // ant1 > ant2 replicated CPU code in CPacerImager::gridding_fast to only use half of the correlation matrix !
             float re = baseline_data[freq_channel];

             int pos = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
             if(pos>=0 && pos<image_size_cuda)
             {
                atomicAdd(&out_visibilities_gridded[pos].x,re);
             }   

             int pos2 = calculate_pos_local( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
             if(pos2>=0 && pos2<image_size_cuda)
             {
                // Allocating inside m_in_buffer as well 
                atomicAdd(&out_visibilities_gridded[pos2].x,re);
             }
          }
       }
    }
}



// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_lfile_vis( int xySize,                   // size of correlation matrix - should really be half of the corr-matrix without diagonal
                                            int freq_channel,             // frequency channel
                                            double wavelength_cuda,       // wavelength [m]
                                            int pol1, char cPol1, int pol2, char cPol2, // which correlation products to grid / image , cPol1 and 2 are to avoid ifs
                                            float *u_cuda, float *v_cuda, // UVW
                                            double delta_u_cuda, double delta_v_cuda, 
//                                            double min_uv_cuda, // MinUV - will try to avoid if-s here TODO : remove other baselines by 0-weighting visibilities in these cells (matrix multiplication -> Tensor Cores ?)
                                            gpufftComplex* in_visibilities_corr, int in_vis_corr_size,
                                            float* in_visibilities_auto, int in_vis_auto_size,
                                            gpufftComplex* out_visibilities_gridded, 
                                            int image_size_cuda, // full image size (usually n_pixels_cuda*n_pixels_cuda)
                                            int n_pixels_cuda,
                                            int n_ant, // number of antennas 
                                            int n_channels,
                                            InputMapping* inputs, int n_inputs, 
                                            int* mapping_array
                                          )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x; // thread index is still pos in 256x256 array of antennas. Howver, in_visibilities_corr only have half of visibilities (lower half of corr-matrix)

    gridding_imaging_lfile_vis_base( i, xySize, freq_channel, wavelength_cuda, pol1, cPol1, pol2, cPol2, u_cuda, v_cuda, delta_u_cuda, delta_v_cuda,
                   	             in_visibilities_corr, in_vis_corr_size, in_visibilities_auto, in_vis_auto_size, out_visibilities_gridded, 
                                     image_size_cuda, n_pixels_cuda, n_ant, n_channels, inputs, n_inputs, mapping_array );
}


// BLOCKS version of the above kernel gridding_imaging_lfile_vis - channels are implemented as separate blocks -> wave length has to be relculated inside 
// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_lfile_vis_blocks( 
                                            int xySize,                   // size of correlation matrix - should really be half of the corr-matrix without diagonal
                                            float first_channel_center_freq_mhz, // center frequency of the 1st fine channel 
                                            float channel_bw_mhz,                // bandwidth of fine channels 
                                            int pol1, char cPol1, int pol2, char cPol2,
                                            float *u_cuda, float *v_cuda, // UVW
                                            double delta_u_cuda, double delta_v_cuda, 
//                                            double min_uv_cuda, // MinUV - will try to avoid if-s here TODO : remove other baselines by 0-weighting visibilities in these cells (matrix multiplication -> Tensor Cores ?)
                                            gpufftComplex* in_visibilities_corr, int in_vis_corr_size, // CROSS-CORRELATIONS
                                            float* in_visibilities_auto, int in_vis_auto_size,        // AUTO-CORRELATIONS
                                            gpufftComplex* out_visibilities_gridded_param, int image_size_cuda, // full image size (usually n_pixels_cuda*n_pixels_cuda)
                                            int n_pixels_cuda, // side of image (number of pixels is n_pixels_cuda*n_pixels_cuda)
                                            int n_ant, // number of antennas 
                                            int n_channels,
                                            InputMapping* inputs, int n_inputs, 
                                            int* mapping_array
                                          )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int freq_channel = blockIdx.y; // second block dimension means IMAGE BLOCK -> here frequency fine channel

    double freqMHz = first_channel_center_freq_mhz + freq_channel*channel_bw_mhz;
    double wavelength_cuda = SPEED_OF_LIGHT/(freqMHz*1e6); // c/freq_in_Hz
    gpufftComplex* out_visibilities_gridded = out_visibilities_gridded_param + freq_channel*image_size_cuda;


    gridding_imaging_lfile_vis_base( i, xySize, freq_channel, wavelength_cuda, pol1, cPol1, pol2, cPol2, u_cuda, v_cuda, delta_u_cuda, delta_v_cuda,
                   	             in_visibilities_corr, in_vis_corr_size, in_visibilities_auto, in_vis_auto_size, out_visibilities_gridded, 
                                     image_size_cuda, n_pixels_cuda, n_ant, n_channels, inputs, n_inputs, mapping_array );
}