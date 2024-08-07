/* 
Start: 20/03/23 
End  : 21/03/23 

Description: fftw_plan_many() code works! 

Sources: 
for fftw_plan_many_dft: https://www.fftw.org/fftw3_doc/Advanced-Complex-DFTs.html

*/

#include <stdio.h> 
#include <stdlib.h> 
#include <iostream>
using namespace std;

#include <math.h>
#include <string.h> 
#include <time.h>
#include <vector>

// For handling .fits files 
#include "bg_fits.h" 
#include <bg_globals.h>

// In order to concatenate two strings 
#include <bits/stdc++.h> 
// #include <cstring> 

// In order to use high resolution clock 
#include <ctime>
#include <ratio>
#include <chrono>
#include <fftw3.h>

#include "../src/pacer_imager.h"

// In order to use #pragma acc parallel loop directive
// #include <openacc.h>

bool gWriteFits=false; 
string u = "u.fits"; 
string v = "v.fits"; 
string w = "w.fits"; 
string vis_real = "vis_real.fits"; 
string vis_imag = "vis_imag.fits"; 
int N = 1; // default number of blocks =  1
int n_pixels = 180; // default image pixels = 180 
int n_channels = 1; // default number of channels of the code = 1 
int nStreams = 1; // number of CUDA streams (like queues for kernel executions)
bool gDebuggerCalculateControlSum=false;

// constant UVW (like for EDA2) -> does not require recaculation of UVW grid for every single timestamp:
bool gConstantUVW=true; // for EDA2

#define VEL_LIGHT  299792458.0
#define NTHREADS 1024

  // fft_shift(): Taken from Marcin's code 
  void fft_shift( CBgFits& dirty_image, CBgFits& out_image )
  {
    int xSize = dirty_image.GetXSize();
    int ySize = dirty_image.GetYSize();

    CBgFits tmp_image( xSize, ySize );
    
    int center_freq_x = int( xSize/2 );
    int center_freq_y = int( ySize/2 );
    
    int is_odd = 0;
    if ( (xSize%2) == 1 && (ySize%2) == 1 )
    {
        is_odd = 1;
    }

    for(int y=0;y<ySize;y++)
    { 
        float* tmp_data = tmp_image.get_line(y);
        float* image_data = dirty_image.get_line(y);
        
        for(int x=0;x<=center_freq_x;x++)
        { 
          tmp_data[center_freq_x+x] = image_data[x];
        }
        for(int x=(center_freq_x+is_odd);x<xSize;x++)
        {
          tmp_data[x-(center_freq_x+is_odd)] = image_data[x];
        }      
    }

    for(int x=0;x<xSize;x++)
    { 
        for(int y=0;y<=center_freq_y;y++)
        { 
          out_image.setXY(x,center_freq_y+y,tmp_image.getXY(x,y));
        }
        for(int y=(center_freq_y+is_odd);y<ySize;y++)
        {
          out_image.setXY( x , y-(center_freq_y+is_odd),tmp_image.getXY(x,y));
        }      
    }
}

void usage()
{
    printf("fftw_blocks OPTIONS\n");
    printf("-u u.fits : input FITS file with U values [default %s]\n",u.c_str());
    printf("-v v.fits : input FITS file with V values [default %s]\n",v.c_str());
    printf("-w w.fits : input FITS file with W values [default %s]\n",w.c_str());
    printf("-r VIS_REAL : inputs FITS file with REAL part of visibilities [default %s]\n",vis_real.c_str());
    printf("-i VIS_IMAG : inputs FITS file with IMAG part of visibilities [default %s]\n",vis_imag.c_str());
    printf("-F WRITE_OUTPUT_FITS_FILES : write output FITS files [default %d]\n",gWriteFits);
    printf("-n N_BLOCKS : number of blocks to test [default %d]\n",N);
    printf("-f N_CHANNELS : number of frequency channels (or iterations) [default %d]\n",n_channels);
    printf("-p SIZE : size of image (on side) -> full number of pixels is SIZE x SIZE [default SIZE = %d]\n",n_pixels);
    printf("-c : enable calculation of control sum on gridded visibilities to check if they are not zero [default %d]\n",gDebuggerCalculateControlSum);
    printf("-s NUMBER_OF_STREAMS : number of CUDA streams used in gridding [default %d]\n",nStreams);
    printf("-M : changing UVW like for the MWA [default %d , i.e. constant like for all-sky images like for EDA2/SKA-Low statations]\n",gConstantUVW);

    exit(0);
}

// In order to pass command line arguments with - option
void parse_cmdline(int argc, char * argv[]) 
{
  char optstring[] = "cn:p:f:r:i:u:v:w:F:s:M";
  int opt;

  while ((opt = getopt(argc, argv, optstring)) != -1) 
  {
     switch (opt) 
     {
        case 'M':
           gConstantUVW = false;
           break;

        case 'c':
           gDebuggerCalculateControlSum = true;
           break;

        case 'u':
           if( optarg )
           {
              u = optarg;
           }
           break;
           
        case 'v':
           if( optarg )
           {
              v = optarg;
           }
           break;
           
        case 'w':
           if( optarg )
           {
              w = optarg;
           }
           break;

        case 'n':
           if( optarg )
           {
              N = atoi(optarg);
           }
           break;
        
        case 'f':
           if( optarg )
           {
              n_channels = atoi(optarg);
           }
           break;
        
        case 'p':
           if( optarg )
           {
              n_pixels = atoi(optarg);
           }
           break;

        case 'r':
           if( optarg )
           {
              vis_real = optarg;
           }
           break;
        
        case 'i':
           if( optarg )
           {
              vis_imag = optarg;
           }
           break;

        case 'F':
           if( optarg )
           {
              gWriteFits = (atol(optarg)>0);
           }
           break;
           
        case 's':
           if( optarg )
           {
              nStreams = atol(optarg);
           }
           break;
           
        default:  
           fprintf(stderr,"Unknown option %c\n",opt);
           exit(0); 
     }
  }

/*  if( (n_pixels%2) != 0 ){
     printf("ERROR : only even image sizes are allowed in this version (to optimise kernel), change value of -p parameter to EVEN !\n");
     exit(-1);
  }*/
}      

/*
argc: Number of command line arguments 
argv: List of command line arguments 
*/
int main(int argc, char* argv[])
{      
    printf("\n"); 
    printf("\n START CODE!"); 

    using namespace std::chrono;
    if( argc>=2 && strncmp(argv[1],"-h",2)==0 )
    {
       usage();
    }  
  
    // Values specific to this program 
    double frequency_MHz = 159.375; // default value in MHz 
    double frequency_Hz = frequency_MHz*1e6;
    double wavelength = VEL_LIGHT/frequency_Hz;
    printf("\n OK frequency_Hz: %.4f, wavelength: %.4f",frequency_Hz, wavelength); 

    // Printing number of command line arguments 
    cout << "\n You have entered " << argc << " arguments:" << "\n";

    // Printing the list of command line arguments 
    for (int i = 0; i < argc; ++i)
        cout << i+1 << ":" << argv[i] << "\n";

    parse_cmdline(argc,argv);
   
    // Reading .fits files: 
    printf("\n Reading in .fits files.."); 
    // Input .fits 
    CBgFits u_fits;
    CBgFits v_fits;
    CBgFits w_fits;
    CBgFits vis_real_fits; 
    CBgFits vis_imag_fits; 

    u_fits.ReadFits( u.c_str() , 0, 1, 1 ); 
    v_fits.ReadFits( v.c_str() , 0, 1, 1 );
    w_fits.ReadFits( w.c_str() , 0, 1, 1 );
    vis_real_fits.ReadFits( vis_real.c_str() , 0, 1, 1 );
    vis_imag_fits.ReadFits( vis_imag.c_str() , 0, 1, 1 );

    // Input size: u, v and w
    int u_xSize = u_fits.GetXSize();
    int u_ySize = u_fits.GetYSize();
    int xySize = (u_xSize*u_ySize); // 256x256 for EDA2 
    printf("\n OK xySize (u,v,w size) = %d", xySize);

    // Input size: vis_real 
    int vis_real_xSize =  vis_real_fits.GetXSize(); 
    int vis_real_ySize =  vis_real_fits.GetYSize();
    int vis_real_size = (vis_real_xSize*vis_real_ySize); // 256x256 for EDA2 
    printf("\n OK vis_real_size = %d", vis_real_size); 
 
    // Input size: vis_imag 
    int vis_imag_xSize =  vis_imag_fits.GetXSize(); 
    int vis_imag_ySize =  vis_imag_fits.GetYSize();  
    int vis_imag_size = (vis_real_xSize*vis_real_ySize); // 256X156 for EDA2 
    printf("\n OK vis_imag_size = %d", vis_imag_size); 
 
    // Image dimensions 
    int width = n_pixels; 
    int height = n_pixels;
    int image_size = (width*height);
    int block_size =  (N*image_size); 
    printf("\n OK PARAMETERS : "); 
    printf("\n OK Number of blocks: %d", N); 
    printf("\n OK Number of channels: %d", n_channels); 
    printf("\n OK n_pixels: %d", n_pixels); 
    printf("\n OK Image Width: %d", width); 
    printf("\n OK Image Height: %d", height); 
    printf("\n OK Overall Image Size: %d", image_size); 
    printf("\n OK Overall Block size: %d", block_size); 
    printf("\n OK Number of CUDA streams: %d", nStreams);
    printf("\n OK Constant UVW      : %d", gConstantUVW);

    // CONSTANTS once :
    // delta_u, delta_v calculations 
    double FOV_degrees = 180.00;
    double FoV_radians = FOV_degrees*M_PI/180; 
    double delta_u = 1.00/(FoV_radians);
    double delta_v = 1.00/(FoV_radians);
    printf("\n OK M_PI: %f",M_PI); 
    printf("\n OK FOV_degrees: %f", FOV_degrees); 
    printf("\n OK FOV_radians: %f", FoV_radians); 
    printf("\n OK delta_u, delta_v (C++) : %f %f", delta_u, delta_v); 
    int center_x = int(n_pixels/2);
    int center_y = int(n_pixels/2);
    double min_uv = -1000;

    CBgFits uv_grid_real_fits(width, height); 
    CBgFits uv_grid_imag_fits(width, height); 
    CBgFits uv_grid_counter_fits(width, height); 
    uv_grid_real_fits.SetValue( 0.00 );
    uv_grid_imag_fits.SetValue( 0.00 );
    uv_grid_counter_fits.SetValue( 0.00 );

    // GAYATRI CHANGE Step 1: CPU Input/Output Variables 
    fftw_complex* m_in_buffer;
    fftw_complex* m_out_buffer;

    // GAYATRI CHANGE Step 2: Allocating m_in_buffer buffer/m_out_buffer 
    m_in_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image_size*N); // same as gridded visibilities 
    m_out_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image_size*N);

    // initialise input visibilities :
    fftw_complex* input_visibilities = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * xySize * N);
    fftw_complex* gridded_visibilities = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image_size * N);
    fftw_complex* gridded_visibilities_for_counter = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image_size * N);
    fftw_complex* sky_image = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * image_size * N);
    
    memset( input_visibilities, '\0', sizeof(fftw_complex) * xySize * N );
    memset( gridded_visibilities, '\0', sizeof(fftw_complex) * image_size * N );
    memset( sky_image, '\0', sizeof(fftw_complex) * image_size * N );
    
    int xSize = u_fits.GetXSize();
    int ySize = u_fits.GetYSize();

    for(int b=0; b<N; b++)
    {
       fftw_complex* ptr_b = input_visibilities + xySize * b;
       
       for(int ant1=0;ant1<xSize;ant1++)
       {
          for(int ant2=0;ant2<xSize;ant2++)
          {
             int pos = ant2*xSize + ant1;
             ptr_b[pos][0] = vis_real_fits.getXY(ant1,ant2);
             ptr_b[pos][1] = vis_imag_fits.getXY(ant1,ant2);
             
          }
       }
    }

    // TODO : sort out linking :
    if( nStreams > 1 )
    {
       fftwf_init_threads();
       fftwf_plan_with_nthreads( nStreams ); // Streams -> FFT threads here 
    }

    // FFT Plan_many :
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    // Define input and output arrays
    int n[2]; // array sizes
    n[0] = width; 
    n[1] = height; 
    int rank = 2; // number of dimensions

    // cufftPlanMany(&plan, 2, n, NULL, 1, image_size, NULL, 1, image_size, CUFFT_C2C, N);
    // cufftPlanMany(cufftHandle *plan, int rank, int *n, int *inembed, int istride, int idist, int *onembed, int ostride, int odist, cufftType type, int batch);

    // Create FFTW plan
    fftw_plan plan = fftw_plan_many_dft(
                  rank, 
                  n, 
                  N,
                  gridded_visibilities,
                  NULL, 
                  1, 
                  image_size, 
                  sky_image,
                  NULL, 
                  1,
                  image_size, 
                  FFTW_FORWARD, FFTW_ESTIMATE);
    
    high_resolution_clock::time_point t1a = high_resolution_clock::now();
    duration<double> time_span1 = duration_cast<duration<double>>(t1a - t1);
    printf("\n CLOCK fftw_plan_many_dft() took: %.6f seconds. PARAMETERS ( N_PIXELS , N_BLOCKS , N_STREAMS , N_CHANNELS ) = ( %d , %d , %d , %d ) \n",time_span1.count(),n_pixels,N,nStreams,n_channels);
    
    // statistics:    
    double gridding_sum=0.00,gridding_sum2=0.00;
    int gridding_count=0;
    double fftexec_sum=0.00,fftexec_sum2=0.00;
    int fftexec_count=0;


    CPacerImager pacer_imager_cpu;
    // Iterating over number of frequency channels 
    for(int c=0; c<n_channels; c++)
    {
      // re-initialise memory in UV grid before every iteration (other wise the previous gridding is there !)

      // START: visibilities into N blocks // yo dude

      // START TIME: gridding()
      high_resolution_clock::time_point t1 = high_resolution_clock::now();
      // Iterating through every block        
      for(int b=0; b<N; b++)
      {
         fftw_complex* ptr_input = input_visibilities + xySize * b;
         fftw_complex* ptr_output = gridded_visibilities + image_size * b;
         fftw_complex* ptr_sky = sky_image + image_size * b;
         
         // Call to CPU gridding function 
         pacer_imager_cpu.gridding_fast( ptr_input, u_fits, v_fits, w_fits, ptr_output, uv_grid_counter_fits, delta_u, delta_v, frequency_MHz, n_pixels, min_uv, "" ); 
         
         // fftw_plan pFwd = fftw_plan_dft_2d( width, height, (fftw_complex*)ptr_output, (fftw_complex*)ptr_sky, FFTW_FORWARD, FFTW_ESTIMATE); // was FFTW_FORWARD or FFTW_BACKWARD ???
         // fftw_execute(pFwd);
         // fftw_destroy_plan(pFwd);
      }
      high_resolution_clock::time_point t2 = high_resolution_clock::now();
      duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
      printf("\n CLOCK Gridding took: %.6f seconds. PARAMETERS ( N_PIXELS , N_BLOCKS , N_STREAMS , N_CHANNELS ) = ( %d , %d , %d , %d ) \n",time_span.count(),n_pixels,N,nStreams,n_channels);
      
      gridding_sum += time_span.count()*1e3;
      gridding_sum2 += (time_span.count()*1e3)*(time_span.count()*1e3);
      gridding_count++;
      
      // calculate counter :
      uv_grid_counter_fits.SetValue(0.00);
      pacer_imager_cpu.gridding_fast( input_visibilities, u_fits, v_fits, w_fits, gridded_visibilities_for_counter, uv_grid_counter_fits, delta_u, delta_v, frequency_MHz, n_pixels, min_uv, "" );
       
      // START: fftw_plan_many_dft()
      high_resolution_clock::time_point t1A = high_resolution_clock::now();

      // Execute FFTW plan
      fftw_execute_dft(plan, gridded_visibilities, sky_image);

      printf("\n OK fftw_execute_dft() executed!"); 

      // END: fftwPlanMany() 
      high_resolution_clock::time_point t2A = high_resolution_clock::now();
      duration<double> time_span2 = duration_cast<duration<double>>(t2A - t1A);
      printf("\n CLOCK fftw_execute_dft took: %.6f seconds. PARAMETERS ( N_PIXELS , N_BLOCKS , N_STREAMS , N_CHANNELS ) = ( %d , %d , %d , %d ) \n",time_span2.count(),n_pixels,N,nStreams,n_channels);

      fftexec_sum += time_span2.count()*1e3;
      fftexec_sum2 += (time_span2.count()*1e3)*(time_span2.count()*1e3);
      fftexec_count++;



      double fnorm = 1.00/uv_grid_counter_fits.Sum(); 
      printf("\n fnorm = %f", fnorm);       
        
      // For storing the final outputs 
      CBgFits out_image_real(width, height); 
      CBgFits out_image_imag(width, height); 
      CBgFits out_image_real_shifted(width, height); 
      CBgFits out_image_imag_shifted(width, height); 

      char filename_real[1024]; 
      char filename_imag[1024]; 

      for(int b=0;b<N;b++)
      {
        float* out_data_real = out_image_real.get_data();
        float* out_data_imag = out_image_imag.get_data();

        fftw_complex* ptr_output = sky_image + image_size * b;

        for(int i=0;i<image_size;i++)
        {
          out_data_real[i] =  ptr_output[i][0]*fnorm; 
          out_data_imag[i] =  ptr_output[i][1]*fnorm; 
        }
        
        fft_shift(out_image_real,out_image_real_shifted);
        fft_shift(out_image_imag,out_image_imag_shifted);

        if( gWriteFits )
        {
           sprintf(filename_real,"re_%d%d.fits",c,b);
           printf("\n filename_real saved: %s", filename_real); 
           out_image_real_shifted.WriteFits(filename_real);

           sprintf(filename_imag,"im_%d%d.fits",c,b);
           printf("\n filename_imag saved: %s", filename_imag); 
           out_image_imag_shifted.WriteFits(filename_imag);
        }
      }
    

    }
    
    // report :    
    double gridding_mean = gridding_sum/gridding_count;
    printf("\nGridding took %.6f +/- %.6f [ms]\n",gridding_mean,sqrt(gridding_sum2/gridding_count - gridding_mean*gridding_mean));
    double fftexec_mean = fftexec_sum/fftexec_count;
    printf("\nFFT Exec took %.6f +/- %.6f [ms]\n",fftexec_mean,sqrt(fftexec_sum2/fftexec_count - fftexec_mean*fftexec_mean));

    // Free memory 
    fftw_free( sky_image );
    fftw_free( input_visibilities );
    fftw_free( gridded_visibilities );
    fftw_free( gridded_visibilities_for_counter );
    fftw_free(m_in_buffer); 
    fftw_free(m_out_buffer); 
    
    printf("\n END CODE"); 
    printf("\n"); 

    // End of main() function 
    return 0;
}

   