//------------------------------------------------------------------------------------------------------------------------------------------------
// 1st version of PACER BLINK fast imager : produces a dirty images for now (should work for both MWA and SKA-Low stations)
//------------------------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>

#include <bg_globals.h>
#include "bg_fits.h"
#include <mystring.h>
#include <myparser.h>
#include <mystrtable.h>

#include <vector>
using namespace std;

// FFTW :
#include <fftw3.h>

// defines :
#include "../src/pacer_imager_defs.h"

// Pacer imager class :
#include "../src/pacer_imager.h"

#ifdef _PACER_PROFILER_ON_
#include <mydate.h>
#endif

string in_basename="1276619416_20200619164456";
string gPostfix="";
string out_fits="1276619416_20200619164456_dirty_image.fits";
// string out_rms_fits="out_rms.fits";

// REAL and IMAG external UV GRID files to make dirty image:
string g_in_fits_file_uv_re;
string g_in_fits_file_uv_im;

// options to enable particular functions:
bool gDoGridding=true;
bool gDo_dirty_image=true;

int gImageSize = 128; // was 512

// weighting :
string gUniformWeighting="N"; // uniform, means that sum of visibilities in a given UV cell is divided by a number of points added

double gFrequencyMHz = 159.375; // default of the station in channel 204 is 204*(400.00/512.00) = 159.375 MHz 
double gFOV_degrees  = 180.00;
double gMinUV        = -1000;
bool   gZenithImage  = false;
int    gProfilesNImages = 1;

// input files :
string gAntennaPositionsFile;
string gFlaggedAntennasListString;
vector<int> gFlaggedAntennasList;

// output files:
string gOutputDirectory="./";

// statistics:
string gStatPositionRadius;

void usage()
{
   printf("pacer_dirty_image VISIBILITY_FITS_BASENAME\n\n\n");
   
   printf("Exected input visibility FITS files are VISIBILITY_FITS_BASENAME_vis_realPOSTFIX.fits and VISIBILITY_FITS_BASENAME_vis_imagPOSTFIX.fits , where POSTFIX is specified by option -p\n");
   
   printf("\t-O OUTDIR  : full path to output directory where FITS files are saved [default %s]\n",gOutputDirectory.c_str());
   printf("\t-p POSTFIX : default is not postfix it's added to basename specified with VISIBILITY_FITS_BASENAME\n");
   printf("\t-g 0 / 1 for GRIDDING ON / OFF [default GRIDDING = %d]\n",gDoGridding);
   printf("\t-r FITS_FILE_RE : specify external UV grid real fits file [requires gridding=0]\n");
   printf("\t-i FITS_FILE_IM : specify external UV grid imag fits file [requires gridding=0]\n");
   printf("\t-f FREQ_MHz     : frequency in MHz [default %.3f MHz]\n",gFrequencyMHz);
   printf("\t-F FoV[deg]     : field of view in degrees [default %.4f degree]\n",gFOV_degrees);
   printf("\t-w WEIGHTING    : change weighting schema N - natural, U - uniform [default %s]\n",gUniformWeighting.c_str());
   printf("\t-m MIN_UV_DISTANCE : minimum UV distance in wavelengths for a baseline to be included [default %.4f]\n",gMinUV);
   printf("\t-a antenna_positions.txt : text file with antenna positions in a format : AntName X[m] Y[m] Z[m]\n");
   printf("\t-A flagged antennas list string (coma separated, e.g. 1,2,3,4) [default empty]\n");  
   printf("\t-n IMAGE_SIZE : single value N for image size N x N pixels [default %d]\n",gImageSize);
   printf("\t-Z : image phase centered at zenith, re-calculation of UVW is not required\n");
   printf("\t-P N_TIMES : create the same image N times for profiling or other tests [default %d]\n",gProfilesNImages);
   printf("\t-v VERBOSITY/DEBUG level [default %d]\n",CPacerImager::m_ImagerDebugLevel);
   printf("\t-V FILE_SAVE_LEVEL : to control number of FITS files saved [default %d]\n",CPacerImager::m_SaveFilesLevel);
   printf("\t-s : print statistics of the final sky image (mean,median,rms,rms_iqr etc)\n");
   printf("\t-S Xc,Yc,R : specify circular window of radius R around (Xc,Yc) position where additional statistics (as -s) are calculated\n");   
   printf("\t-c CAL_SOL_FILE : specify file with calibration solutions [default None], currently supported formats : .txt\n");
   printf("\t-I : include auto-correlations [not included by default]\n");
   printf("\t-M META_DATA_FILE : name of meta data file\n");
   printf("\t-o OPTION  : additional options specified by name, such as:\n");
   printf("\t             miriad : uses the exact settings of MIRIAD natural weighting for validation\n");

   exit(0);
}

void parse_cmdline(int argc, char * argv[], CPacerImager& imager) {
   char optstring[] = "hp:n:g:r:i:f:F:w:m:a:ZP:v:V:A:sS:o:O:c:IM:";
   int opt;
        
   while ((opt = getopt(argc, argv, optstring)) != -1) {
//      printf("opt = %c (%s)\n",opt,optarg);   
      switch (opt) {
         case 'h':
            // antenna1 = atol( optarg );
            usage();
            break;

         case 'a':
            if( optarg && strlen(optarg) ){
               gAntennaPositionsFile = optarg;
            }
            break;

         case 'A':
            if( optarg && strlen(optarg) ){
               gFlaggedAntennasListString = optarg;
            }
            break;

         case 'c':
            imager.m_CalibrationSolutions.m_filename = optarg;
            break;

         case 'f':
            gFrequencyMHz = atof( optarg );
            break;

         case 'F':
            gFOV_degrees = atof( optarg );
            break;

         case 'I':
            imager.m_bIncludeAutos = true;
            break;

         case 'M':
            imager.m_ImagerParameters.m_MetaDataFile = optarg;
            imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
            gZenithImage = false;
            break;

         case 'p':
            gPostfix = optarg;
            break;

         case 'P':
            gProfilesNImages = atol( optarg );
            break;

         case 'n':
            gImageSize = atol( optarg );
            break;

         case 'g':
            gDoGridding = atol( optarg );
            break;
            
         case 'i':
            g_in_fits_file_uv_im = optarg;
            break;

         case 'm':
            gMinUV = atof( optarg );
            break;

         case 'o':
            if( optarg ){
               if( strcmp(optarg,"miriad")==0 || strcmp(optarg,"MIRIAD")==0 ){
                  CPacerImager::m_bCompareToMiriad = true;
               }
            }
            break;

         case 'O' :
            if( optarg ){
               gOutputDirectory = optarg;
            }
            break;
         case 'r' :
            g_in_fits_file_uv_re = optarg;
            break; 
         
         case 's' :
            CPacerImager::m_bPrintImageStatistics = true;
            break;
            
         case 'S' :
            if( optarg && strlen(optarg) ){
               gStatPositionRadius = optarg;
            }
            break;

         case 'v' :
            if( optarg ){
               CPacerImager::SetDebugLevel( atol( optarg ) );
            }
            break; 

         case 'V' :
            if( optarg ){
               CPacerImager::SetFileLevel( atol( optarg ) );
            }
            break; 

         case 'w' :
            if( optarg ){
               gUniformWeighting = optarg;
            }
            break; 

         case 'Z' :
            gZenithImage = true;
            imager.m_ImagerParameters.m_bConstantUVW = true; // constant UVW (not re-calculated as a function of time due to pointing)
            break; 

         default:   
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
   
   if( strlen(gFlaggedAntennasListString.c_str()) > 0 ){
      MyParser pars=gFlaggedAntennasListString.c_str();
      CMyStrTable items;
      pars.GetItemsNew( items, "," );
      for(int i=0;i<items.size();i++){
         gFlaggedAntennasList.push_back( atol( items[i].c_str() ) );
      }
   }
}

void print_parameters( CPacerImager& imager )
{
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");
    printf("Base name for FITS  = %s\n",in_basename.c_str());
    printf("out_fits     = %s\n",out_fits.c_str());
    printf("Output dir   = %s\n",gOutputDirectory.c_str());
    printf("Postfix      = %s\n",gPostfix.c_str());
    printf("Image size   = (%d x %d)\n",gImageSize,gImageSize);
    printf("Do gridding  = %d\n",gDoGridding);
    if( !gDoGridding ){
       printf("External REAL/IMAG FITS FILES : %s / %s\n",g_in_fits_file_uv_re.c_str(),g_in_fits_file_uv_im.c_str());
    }
    printf("Frequency           = %.4f [MHz]\n",gFrequencyMHz);
    printf("Field of view (FoV) = %.4f [deg]\n",gFOV_degrees);
    printf("Weighting           = %s\n",gUniformWeighting.c_str());
    printf("UV range            = %.4f - Infinity\n",gMinUV);   
    printf("Antenna positions file = %s\n",gAntennaPositionsFile.c_str());
    printf("Meta data file      = %s\n",imager.m_ImagerParameters.m_MetaDataFile.c_str());
    printf("Zenith image        = %d\n",gZenithImage);
    printf("Profile N times the same image = %d\n",gProfilesNImages);
    printf("Verbosity level     = %d\n",CPacerImager::m_ImagerDebugLevel);
    printf("File save level     = %d\n",CPacerImager::m_SaveFilesLevel);
    if( gFlaggedAntennasList.size() > 0 ){
       mystring szFlaggedAnts;
       for(int i=0;i<gFlaggedAntennasList.size();i++){
          char szTmp[16];
          sprintf(szTmp,"%d",gFlaggedAntennasList[i]);
          
          szFlaggedAnts += szTmp;
          szFlaggedAnts += ",";
       }
       printf("Flagged antennas    = %s\n",szFlaggedAnts.c_str());
    }
    printf("Show statistics = %d\n",CPacerImager::m_bPrintImageStatistics);
    if( strlen(gStatPositionRadius.c_str()) > 0 ){
       printf("Show statistics around position = %s\n",gStatPositionRadius.c_str());
    }
    printf("MIRIAD compatibility = %d\n",CPacerImager::m_bCompareToMiriad);
    printf("Calibration solutions file = %s\n",imager.m_CalibrationSolutions.m_filename.c_str());
    printf("Include auto-correlations  = %d\n",imager.m_bIncludeAutos);
    printf("############################################################################################\n");
}

void print_header()
{
#ifdef __BUILD_DATE__
   printf("pacer_imager version 1.00 (C++/object), compiled on %s\n",__BUILD_DATE__);
#else
   // this is becuse __TIMESTAMP__ seems to be stuck on the 19th March and does not change ???!!!
   printf("pacer_imager version 1.00 (C++/object), compiled on %s\n",__TIMESTAMP__);
#endif   
}

int main(int argc,char* argv[])
{
   print_header();

  // parsing paratemeters with input and output file names:
  if( argc >= 2 ){
     in_basename = argv[1];
  }
  out_fits="out.fits";
  if( argc >= 3 ){
     out_fits = argv[2];
  }  

  CPacerImager imager;
  
  // parse and print parameters :
  parse_cmdline( argc , argv , imager );
  print_parameters( imager );
  
  imager.m_ImagerParameters.SetGlobalParameters( gAntennaPositionsFile.c_str(), (gZenithImage) ); // Constant UVW when zenith image (-Z)
  imager.m_ImagerParameters.m_szOutputDirectory = gOutputDirectory.c_str();
  
  imager.Initialise();
  
  // set antenna flags (if there are):
  if( gFlaggedAntennasList.size() > 0 ){
     imager.SetFlaggedAntennas( gFlaggedAntennasList );
  }

  CBgFits uv_grid_real_param,uv_grid_imag_param;
  
  if( uv_grid_real_param.ReadFits( "layer_re.fits" ) ){
     printf("ERROR : could not read FITS layer_re.fits\n");
     exit(-1);
  }
  
  if( uv_grid_imag_param.ReadFits( "layer_im.fits" ) ){
     printf("ERROR : could not read FITS layer_im.fits\n");
     exit(-1);
  }
  
  printf("TEST values [0,0] = %.8f and [1,0] = %.8f\n", uv_grid_real_param.getXY(0,0),  uv_grid_real_param.getXY(1,0) );
     
  CBgFits uv_grid_counter( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() );
  uv_grid_counter.SetValue(1);
  
  imager.dirty_image( uv_grid_real_param, uv_grid_imag_param, uv_grid_counter, true, "test", true, false );
}


