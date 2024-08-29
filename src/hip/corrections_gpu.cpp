#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <gpu_macros.hpp>
#include "pacer_imager_hip_defines.h"


__global__ void apply_cable_corrections( int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, double *cable_lengths_cuda, double frequency_hz, double speed_of_light )
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
    double im = vis[1]; // - due to upper triangular matrix here (from Cristian's correlator)

    // double cableDeltaLen = 0.00;
    double cableDeltaLen = (cable_lengths_cuda[min_a] - cable_lengths_cuda[max_a]);
    // double w = -w_cuda[i]; // was + but now - for upper triangular matrix
    double angle = -2.0*M_PI*cableDeltaLen*frequency_hz / speed_of_light; // was +  but changed to - due to upper triangular matrix here (from Cristian's correlator)
    double sin_angle,cos_angle;
    sincos(angle, &sin_angle, &cos_angle);
    
    double re_prim = re*cos_angle - im*sin_angle;
    double im_prim = im*cos_angle + re*sin_angle;

    vis[0] = re_prim;
    vis[1] = im_prim; // was + but now - for upper triangular matrix

}

__global__ void apply_geometric_corrections(int xySize, int n_ant, VISIBILITY_TYPE *vis_cuda, float *w_cuda, double frequency_hz, double speed_of_light)
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
    double angle = 2.0*M_PI*w*frequency_hz / speed_of_light; // was +  but changed to - due to upper triangular matrix here (from Cristian's correlator)
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


void apply_cable_lengths_corrections_gpu(Visibilities &xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   cable_lengths.to_gpu();
   frequencies.to_gpu();

   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;

   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         double frequency_hz = frequencies[fine_channel];
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_cable_corrections<<<nBlocks, NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, cable_lengths.data(), frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
      }
   }
}


void apply_geometric_corrections_gpu(Visibilities &xcorr, float *w_gpu, MemoryBuffer<double>& frequencies){
   if(!xcorr.on_gpu()) xcorr.to_gpu();
   frequencies.to_gpu();
   int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
   int nBlocks = (xySize + NTHREADS -1)/NTHREADS;
   for(int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
      for(int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++){
         double frequency_hz = frequencies[fine_channel];
         // TODO: if ( frequency == freq_channel || freq_channel < 0 ){
         VISIBILITY_TYPE* vis_local_gpu = (VISIBILITY_TYPE*)xcorr.at(time_step,fine_channel,0,0);
         apply_geometric_corrections<<<nBlocks,NTHREADS>>>(xySize, xcorr.obsInfo.nAntennas, vis_local_gpu, w_gpu, frequency_hz, SPEED_OF_LIGHT);
         gpuCheckLastError();
      }
   }
}