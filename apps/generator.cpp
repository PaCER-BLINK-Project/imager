// program averages few FITS images of the same sizes 

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include <bg_globals.h>
#include "bg_fits.h"
#include <mystring.h>

#include <vector>
using namespace std;

// FFTW :
#include <fftw3.h>

string in_basename="uv_grid";
string gPostfix="";
string out_fits="1276619416_20200619164456_dirty_image.fits";

int n_pixels = 128; // was 512

# define MAX_VIS 1e12

// parameters :
double gAmplutude = 1.00;
double gOmega = 1.00; // number of full periods in the n_pixel x n_pixel image maximum which makes sense is n_pixel/2 for even or (n_pixel-1)/2 for odd
double gDC = 0.00;

void usage()
{
   printf("pacer_generate_vis OUT_VISIBILITY_FITS_BASENAME [default %s]\\n\n\n",in_basename.c_str());
   
   printf("-A amplitude of cos/sin [default %.4f]\n",gAmplutude);
   printf("-D DC_TERM [default %.8f]\n",gDC);
   printf("-O or -w OMEGA [default %.4f]\n",gOmega);
   
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "hn:A:D:O:w:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);   
      switch (opt) {
         case 'h':
            // antenna1 = atol( optarg );
            usage();
            break;


         case 'A' :
            gAmplutude = atof( optarg );
            break; 
            
         case 'D' :
            gDC = atof( optarg );
            break; 

         case 'O' :
            gOmega = atof( optarg );
            break; 

         case 'w' :
            gOmega = atof( optarg );
            break; 
            
         case 'n':
            n_pixels = atol( optarg );
            break;

         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
}

void print_parameters()
{
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");
    printf("Base name for FITS  = %s\n",in_basename.c_str());
    printf("Amplitude = %.8f\n",gAmplutude);
    printf("DC term   = %.8f\n",gDC);
    printf("############################################################################################\n");
}

long int gSeed=458249083;

double getrandom()
{
    if( gSeed <= 0 ){
       gSeed = 477737 + get_dttm();
       //       srand(gSeed);
    }
    double r = double(rand())/double(RAND_MAX+1.0);

    return r;
}



void generate_uv_grid( CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double amplitude=1.00, double omega=1.00 )
{
   int xSize = uv_grid_real.GetXSize();
   int ySize = uv_grid_real.GetYSize();

   int x_center = xSize/2;
   int y_center = ySize/2;

   for(int iy=0;iy<ySize;iy++){
      for(int ix=0;ix<xSize;ix++){
         double x = 2*M_PI*( ix - x_center )/xSize;
         double y = 2*M_PI*( iy - y_center )/ySize;
         double noise = getrandom();
//         noise = 0.00;
         
         double re_value = amplitude*cos( omega*x ) + amplitude*cos( omega*y ) + gDC ; // + noise*0.1;
         double im_value = amplitude*sin( omega*x ) + amplitude*cos( omega*y ) + gDC; // + noise*0.1;
         
//         double re_value = noise;
//         double im_value = noise;
         
         uv_grid_real.setXY( ix,iy, re_value );
         uv_grid_imag.setXY( ix,iy, im_value );      
      }
   }   
   
   char fits_file_re[1024],fits_file_im[1024],fits_file_counter[1024];
   
   sprintf(fits_file_re,"%s_re.fits",in_basename.c_str());
   sprintf(fits_file_im,"%s_im.fits",in_basename.c_str());
   sprintf(fits_file_counter,"%s_counter.fits",in_basename.c_str());
   
   uv_grid_real.WriteFits( fits_file_re );
   uv_grid_imag.WriteFits( fits_file_im );
   uv_grid_counter.WriteFits( fits_file_counter );
}

int main(int argc,char* argv[])
{
//  if( argc >= 2 ){
//     in_basename = argv[1];
//  }
  parse_cmdline( argc , argv );
  print_parameters();

  CBgFits uv_grid_counter( n_pixels, n_pixels ),uv_grid_real( n_pixels, n_pixels ) , uv_grid_imag( n_pixels, n_pixels );  
  uv_grid_counter.SetValue(1.00);
  
  generate_uv_grid( uv_grid_real, uv_grid_imag, uv_grid_counter, gAmplutude, gOmega );  
}

