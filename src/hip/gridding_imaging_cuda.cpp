/* 
Header file for both gridding + cuFFT, implemented in cuda 
- gridding_imaging_cuda

References: 
https://stackoverflow.com/questions/17489017/can-we-declare-a-variable-of-type-cufftcomplex-in-side-a-kernel
*/

#include <gpu_macros.hpp>
#include <gpu_fft.hpp>

// msfitslib :
#include <array_config_common.h>
#include "../pacer_imager_defs.h"

__device__ int calculate_pos(float u,
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


// Cuda kernal: gridding and cuFFT 
__global__ void gridding_imaging_cuda( int xySize, // size of the correlation matrix
                                      int n_ant,   // number of antennas 
                                      float *u_cuda, float *v_cuda, 
                                      int* antenna_flags, float* antenna_weights,
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      float *uv_grid_counter_cuda, float *uv_grid_real_cuda, float *uv_grid_imag_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda)
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 
    
    // int n_ant = int( sqrt(xySize) );
    int ant1 = i % n_ant;
    int ant2 = i / n_ant;

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) && antenna_flags[ant1]<=0 && antenna_flags[ant2]<=0 ) // not NaN and not flagged
    {
        if( ant1 > ant2 ){ // >= to include auto-correlations !
           // grids both ant1,ant2 and ant2,ant1 UV points - this means that
           // only half of the matrix should be done !
           // if the input correlation matrix has values != NaN in the top half
           // -> they will be gridded -> 2x too many UV points !!!
     
           int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
           if(pos>=0 && pos<image_size_cuda)
           {
              // Allocating in uv_grid                
              atomicAdd(&uv_grid_real_cuda[pos],re); // TODO : these may not be necessary in the future (only for debugging and saving, but it can be done conditionally)
              atomicAdd(&uv_grid_imag_cuda[pos],im); // TODO : these may not be necessary in the future (only for debugging and saving, but it can be done conditionally)
              atomicAdd(&uv_grid_counter_cuda[pos],1);

              // Allocating inside m_in_buffer as well 
              atomicAdd(&m_in_buffer_cuda[pos].x,re);
              atomicAdd(&m_in_buffer_cuda[pos].y,im);
           }   

           int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
           if(pos2>=0 && pos2<image_size_cuda)
           {
              atomicAdd(&uv_grid_real_cuda[pos2],re);  // TODO : these may not be necessary in the future (only for debugging and saving, but it can be done conditionally)
              atomicAdd(&uv_grid_imag_cuda[pos2],-im); // TODO : these may not be necessary in the future (only for debugging and saving, but it can be done conditionally)
              atomicAdd(&uv_grid_counter_cuda[pos2],1);

              // Allocating inside m_in_buffer as well 
              atomicAdd(&m_in_buffer_cuda[pos2].x,re);
              atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
           }        
        }           
    }   
}

// Cuda kernal: gridding and cuFFT 
__global__ void gridding_imaging_cuda_xcorr( int xySize, // size of the correlation matrix
                                      int n_ant,
                                      float *u_cuda, float *v_cuda, 
                                      int* antenna_flags, float* antenna_weights,
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      VISIBILITY_TYPE *vis_cuda,  
                                      float *uv_grid_counter_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda)
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    int max_a = (i / n_ant); // row 
    int min_a = (i % n_ant); // col

    // printf("DEBUG : i = %d -> (%d,%d) vs. image_size = %d , n_ant = %d\n",i,ant1,ant2,image_size_cuda,n_ant);

    //if( i == 0 ){
    //   printf("DEBUG : gridding_imaging_cuda_xcorr vis_cuda[0] = %.8f\n",vis_cuda[0]);
    //}

    // TODO : < change to <=     
    if( max_a <= min_a ){ // <= to exclude autos, temporarily including them !
       return;
    }
    // max_a must really be > min_a to ensure the upper triangular matrix !


    // Getting corresponding real and imag visibilities 
    // double re = vis_real_cuda[i]; 
    // double im = vis_imag_cuda[i]; 
    // std::complex<double>* vis = xcorr.at( time_step, fine_channel, i, j );
    // int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2; // *2 because there are REAL/IMAG pairs of values 

    // UPPER TRIANGULAR MATRIX : index of element - skip : (max_a*(max_a+1))/2 + min_a elements, but here *2 (re/im)  *4 (4 correlation products)    
    // ant1=max_a and ant2=min_a :
    int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    // int vis_index = ( (ant1 * (ant1 + 1)) / 2 + ant2 )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    VISIBILITY_TYPE* vis = vis_cuda + vis_index;
    double re = vis[0];
    double im = vis[1];

    // TODO: comment below lines including return:    
    // debugging :
    /*u_cuda[i] = re; // vis_cuda[2*i];
    v_cuda[i] = im; // vis_cuda[2*i+1];    
    // fill also other half of the matrix :
    if( max_a != min_a ){
       // ant1 is max :
       int it = min_a*n_ant + max_a;
       u_cuda[it] = re;
       v_cuda[it] = -im;
    }    
    return;*/

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) && antenna_flags[min_a]<=0 && antenna_flags[max_a]<=0 ) // not NaN and not flagged
    {
        int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
        if(pos>=0 && pos<image_size_cuda)
        {
           // Allocating in uv_grid                
           atomicAdd(&uv_grid_counter_cuda[pos],1);

           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos].x,re);
           atomicAdd(&m_in_buffer_cuda[pos].y,im);
        }   

        int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
        if(pos2>=0 && pos2<image_size_cuda)
        {
          atomicAdd(&uv_grid_counter_cuda[pos2],1);
           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos2].x,re);
           atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
        }        
    }   
}

__global__ void apply_cable_corrections( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *cable_lengths_cuda, double frequency_hz, double speed_of_light )
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    int max_a = (i / n_ant); // row 
    int min_a = (i % n_ant); // col
    
    if( max_a >= n_ant || min_a >= n_ant ){
       return;
    }

 //   printf("DEBUG : i = %d -> (%d,%d) vs. n_ant = %d\n",i,max_a,min_a,n_ant);

    //if( i == 0 ){
    //   printf("DEBUG : gridding_imaging_cuda_xcorr vis_cuda[0] = %.8f\n",vis_cuda[0]);
    //}
    
    if( max_a < min_a ){ // auto-correlations are also calibrated (not skipped here)
       return;
    }
    // max_a must really be > min_a to ensure the upper triangular matrix !


    // Getting corresponding real and imag visibilities 
    // double re = vis_real_cuda[i]; 
    // double im = vis_imag_cuda[i]; 
    // std::complex<double>* vis = xcorr.at( time_step, fine_channel, i, j );
    // int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2; // *2 because there are REAL/IMAG pairs of values 

    // UPPER TRIANGULAR MATRIX : index of element - skip : (max_a*(max_a+1))/2 + min_a elements, but here *2 (re/im)  *4 (4 correlation products)    
    // ant1=max_a and ant2=min_a :
    int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    // int vis_index = ( (ant1 * (ant1 + 1)) / 2 + ant2 )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    VISIBILITY_TYPE* vis = vis_cuda + vis_index;
    double re = vis[0];
    double im = -vis[1]; // - due to upper triangular matrix here (from Cristian's correlator)

    // double cableDeltaLen = 0.00;
    double cableDeltaLen = (cable_lengths_cuda[max_a] - cable_lengths_cuda[min_a]);
    // double w = -w_cuda[i]; // was + but now - for upper triangular matrix
    double angle = -2.0*M_PI*cableDeltaLen*frequency_hz / speed_of_light; // was +  but changed to - due to upper triangular matrix here (from Cristian's correlator)
    double sin_angle,cos_angle;
    sincos(angle, &sin_angle, &cos_angle);
    
    double re_prim = re*cos_angle - im*sin_angle;
    double im_prim = im*cos_angle + re*sin_angle;
    
    if( max_a == 1 && min_a == 0 ){    
       printf("KERNEL : (%d,%d) , angle = %.8f , cable_max_a = %.8f , cable_min_a = %.8f, %.4f / %.4f , freq = %.2f Hz -> %.8f / %.8f\n",max_a,min_a,angle,cable_lengths_cuda[max_a],cable_lengths_cuda[min_a],re,im,frequency_hz,re_prim,im_prim);
    }
    
    vis[0] = re_prim;
    vis[1] = -im_prim; // was + but now - for upper triangular matrix

}

__global__ void apply_geometric_corrections( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *w_cuda, double frequency_hz, double speed_of_light )
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    int max_a = (i / n_ant); // row 
    int min_a = (i % n_ant); // col

 //   printf("DEBUG : i = %d -> (%d,%d) vs. n_ant = %d\n",i,max_a,min_a,n_ant);

    //if( i == 0 ){
    //   printf("DEBUG : gridding_imaging_cuda_xcorr vis_cuda[0] = %.8f\n",vis_cuda[0]);
    //}
    
    if( max_a < min_a ){ // auto-correlations are also calibrated (not skipped here)
       return;
    }
    // max_a must really be > min_a to ensure the upper triangular matrix !


    // Getting corresponding real and imag visibilities 
    // double re = vis_real_cuda[i]; 
    // double im = vis_imag_cuda[i]; 
    // std::complex<double>* vis = xcorr.at( time_step, fine_channel, i, j );
    // int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2; // *2 because there are REAL/IMAG pairs of values 

    // UPPER TRIANGULAR MATRIX : index of element - skip : (max_a*(max_a+1))/2 + min_a elements, but here *2 (re/im)  *4 (4 correlation products)    
    // ant1=max_a and ant2=min_a :
    int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    // int vis_index = ( (ant1 * (ant1 + 1)) / 2 + ant2 )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    VISIBILITY_TYPE* vis = vis_cuda + vis_index;
    double re = vis[0];
    double im = -vis[1]; // - due to upper triangular matrix here (from Cristian's correlator)

    double w = -w_cuda[i]; // was + but now - for upper triangular matrix
    double angle = -2.0*M_PI*w*frequency_hz / speed_of_light; // was +  but changed to - due to upper triangular matrix here (from Cristian's correlator)
    double sin_angle,cos_angle;
    sincos(angle, &sin_angle, &cos_angle);
    
    double re_prim = re*cos_angle - im*sin_angle;
    double im_prim = im*cos_angle + re*sin_angle;
    
    if( max_a == 1 && min_a == 0 ){    
       printf("KERNEL : (%d,%d) , angle = %.8f , w = %.8f , %.4f / %.4f , freq = %.2f Hz -> %.8f / %.8f\n",max_a,min_a,angle,w,re,im,frequency_hz,re_prim,im_prim);
    }
    
    vis[0] = re_prim;
    vis[1] = -im_prim; // was + but now - for upper triangular matrix
}


__global__ void vis2corrmatrix( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *vis_corrmatrix_re_cuda, float* vis_corrmatrix_im_cuda )
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    int max_a = (i / n_ant); // row 
    int min_a = (i % n_ant); // col

    if( max_a < min_a ){ // <= to exclude autos, temporarily including them !
       return;
    }

    int vis_index = ( (max_a * (max_a + 1)) / 2 + min_a )*2*4; // *2 because there are REAL/IMAG pairs of values, *4 correlation products !!!
    VISIBILITY_TYPE* vis = vis_cuda + vis_index;
    double re = vis[0];
    double im = vis[1];

    // TODO: comment below lines including return:    
    // debugging :
    vis_corrmatrix_re_cuda[i] = re; // vis_cuda[2*i];
    vis_corrmatrix_im_cuda[i] = im; // vis_cuda[2*i+1];    

    // fill also other half of the matrix :
    if( max_a != min_a ){
       // ant1 is max :
       int it = min_a*n_ant + max_a;
       if( it >= xySize ){
          return;
       }
       vis_corrmatrix_re_cuda[it] = re;
       vis_corrmatrix_im_cuda[it] = -im;
    }    
}

__global__ void gridding_imaging_cuda_optimised(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      float *uv_grid_counter_cuda, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda)
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) )
    {
        int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
        if(pos>=0 && pos<image_size_cuda)
        {
           // Allocating in uv_grid                
           atomicAdd(&uv_grid_counter_cuda[pos],1);

           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos].x,re);
           atomicAdd(&m_in_buffer_cuda[pos].y,im);
        }   

        int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
        if(pos2>=0 && pos2<image_size_cuda)
        {
           // Allocating in uv_grid
           atomicAdd(&uv_grid_counter_cuda[pos2],1);

           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos2].x,re);
           atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
        }
        
    }   
}

// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_cuda_optimised_nocounter(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda)
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) )
    {
        int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
        if(pos>=0 && pos<image_size_cuda)
        {
           atomicAdd(&m_in_buffer_cuda[pos].x,re);
           atomicAdd(&m_in_buffer_cuda[pos].y,im);
        }   

        int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
        if(pos2>=0 && pos2<image_size_cuda)
        {
           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos2].x,re);
           atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
        }
        
    }   
}



// Cuda kernal: gridding and cuFFT 
__global__ void calculate_counter( int xySize, // size of the correlation matrix
                                   float *u_cuda, float *v_cuda, 
                                   double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                   int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                   float *vis_real_cuda, float *vis_imag_cuda, 
                                   float *uv_grid_counter_cuda, double min_uv_cuda )
{   
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) )
    {
         int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
         if(pos>=0 && pos<image_size_cuda)
         {
            // Allocating in uv_grid
            atomicAdd(&uv_grid_counter_cuda[pos],1);
         }   

         int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
         if(pos2>=0 && pos2<image_size_cuda)
         {
            atomicAdd(&uv_grid_counter_cuda[pos2],1);
         }
    }   
}


// gridding kernels using BLOCKS : with and without counter calculation:
__global__ void gridding_imaging_cuda_blocks_optimised(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      float *uv_grid_counter_cuda_param, double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda_param)
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int block = blockIdx.y; // second block dimension means IMAGE BLOCK now 

    gpufftComplex *m_in_buffer_cuda = m_in_buffer_cuda_param + block*image_size_cuda;
    float* uv_grid_counter_cuda = uv_grid_counter_cuda_param  + (block*image_size_cuda);

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) )
    {
        int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
        if(pos>=0 && pos<image_size_cuda)
        {
           // Allocating in uv_grid                
           atomicAdd(&uv_grid_counter_cuda[pos],1);

           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos].x,re);
           atomicAdd(&m_in_buffer_cuda[pos].y,im);
        }   

        int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
        if(pos2>=0 && pos2<image_size_cuda)
        {
           // Allocating in uv_grid
           atomicAdd(&uv_grid_counter_cuda[pos2],1);

           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos2].x,re);
           atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
        }        
    }   
}

// same as above gridding_imaging_cuda_blocks_optimised but counter is not calculated, which is for the case of CONSTANT UVW :
__global__ void gridding_imaging_cuda_blocks_optimised_nocounter(int xySize, // size of the correlation matrix
                                      float *u_cuda, float *v_cuda, 
                                      double wavelength_cuda, int image_size_cuda, double delta_u_cuda, double delta_v_cuda, 
                                      int n_pixels_cuda, int center_x_cuda, int center_y_cuda, int is_odd_x_cuda, int is_odd_y_cuda,
                                      float *vis_real_cuda, float *vis_imag_cuda, 
                                      double min_uv_cuda, 
                                      gpufftComplex *m_in_buffer_cuda_param)
{
    // Calculating the required id 
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int block = blockIdx.y; // second block dimension means IMAGE BLOCK now 

    gpufftComplex *m_in_buffer_cuda = m_in_buffer_cuda_param + block*image_size_cuda;

    if ( i >= xySize ){
       // this thread should not be used as it will correspond to the cell
       // outside the size of the correlation matrix
       return;
    }

    // Getting corresponding real and imag visibilities 
    double re = vis_real_cuda[i]; 
    double im = vis_imag_cuda[i]; 

    // Checking for NaN values 
    if( !isnan(re) && !isnan(im) )
    {
        int pos = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, +1 );
        if(pos>=0 && pos<image_size_cuda)
        {
           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos].x,re);
           atomicAdd(&m_in_buffer_cuda[pos].y,im);
        }   

        int pos2 = calculate_pos( u_cuda[i], v_cuda[i], delta_u_cuda, delta_v_cuda, wavelength_cuda, min_uv_cuda, n_pixels_cuda, center_x_cuda, center_y_cuda, is_odd_x_cuda, is_odd_y_cuda, -1 );
        if(pos2>=0 && pos2<image_size_cuda)
        {
           // Allocating inside m_in_buffer as well 
           atomicAdd(&m_in_buffer_cuda[pos2].x,re);
           atomicAdd(&m_in_buffer_cuda[pos2].y,-im);
        }        
    }   
}

