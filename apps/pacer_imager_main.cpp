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
#include <mydate.h>
#include <myfile.h>

#include <vector>
using namespace std;

// defines :
#include "../src/pacer_imager_defs.h"
#include <images.hpp>

// Pacer imager class :
#ifdef IMAGER_HIP
#include "../src/hip/pacer_imager_hip.h"
// #include "../src/pacer_imager.h"
#else
// FFTW :
#include <fftw3.h>
#include "../src/pacer_imager.h"
#include "../src/gridding.hpp"
#endif

#ifdef _PACER_PROFILER_ON_
#include <mydate.h>
#endif

bool read_corr_matrix( const char* basename, CBgFits& fits_vis_real, CBgFits& fits_vis_imag, const char* szPostfix, int pol=0 );

Visibilities* ConvertFits2XCorr( CBgFits& vis_re, CBgFits& vis_im, ObservationInfo& obsInfo, CBgFits* vis_re_yy=nullptr, CBgFits* vis_im_yy=nullptr );

const char* PolarizationToString( Polarization pol )
{
   switch( pol ){
      case Polarization::XX :
         return "XX";
      case Polarization::YY :
         return "YY";
      case Polarization::I :
         return "I";
         
      default :
         break;   
   }
   
   return "XX";
}

class CImagerParameters
{
public :
   // static parameters, which will be the same for all CPacerImager objects in the current program execution
   // for example antenna positions are the same for all frequency channels (if they are imaged by separate CPacerImager objects as is currently implemented)
   static std::string m_AntennaPositionsFile; // antenna positions are same for all freq. channels -> static, similarly m_bConstantUVW
   static std::string m_MetaDataFile; // metafits file (.txt or .metafits)
   static bool   m_bAutoFixMetaData;  // automatically recalculate RA,DEC,TIME using just standard METAFITS file (no need to create special METAFITS using fix_metafits_time_radec_all.py )
   static bool   m_bConstantUVW; // default false and only true to zenith phase centered images (all-sky from EDA2)
   static bool   m_bAntennaPositionsXYZ; // default false, true when antenna poisitions already in XYZ coordinates system (WG54 or whatever it is)
   static bool   m_bCalcEarthXYZ; // for the MWA we require conversion from local (x,y,z) to Earth's (X,Y,Z) coordinates , for EDA2 we can keep using local (x,y,z) 
   std::string m_szOutputDirectory;   // can save fits files to different directories for different freq. channels
   double m_fUnixTime;
   static bool   m_bApplyGeomCorr;
   static bool   m_bApplyCableCorr;
   
   // frequency and bandwidth 
   double m_fCenterFrequencyMHz;
   double m_fBandwidthMHz;
   
   // which polarisation to image:
   Polarization pol_to_image { Polarization::XX };
   
   // UV range :
   double m_fMinUV;
   
   // output image size :
   int m_ImageSize;
   double m_integrationTime;
   
   // Image Field-of-View :
   double m_ImageFOV_degrees;
   double m_PixsizeInRadians; // pixel size in radians 
   
   std::string m_szWeighting; // U for uniform, N for natural etc 
   
   // W-range 
   // bool m_bW_range_specified;
   double m_MinW; // -INF
   double m_MaxW; // +INF
   
   // gridding and convolution kernels :
   int m_nConvolvingKernelSize; // default <0 -> no convolition kernel, >0 -> kernel size in pixels 
   

   CImagerParameters();
   CImagerParameters( const CImagerParameters& right );
   static void SetGlobalParameters( const char* szAntennaPositionsFile, bool bConstantUVW=false );
   CImagerParameters& operator=( const CImagerParameters& right );
   
   inline void SetConstantUVW(){ m_bConstantUVW = true; }
   inline void SetUnixTime( double fUnixTime ){ m_fUnixTime = fUnixTime; printf("DEBUG : unixtime set to %.6f\n",m_fUnixTime); }
   
//   void SetAntPositionsFile( const char* szAntennaPositionsFile );
   
};

void Initialise( CPacerImager& imager );

string in_basename="1276619416_20200619164456";
string gPostfix="";
bool   gUnixTimeParamProvided=false; // this is to indicate that unix time is explicitly provided and there is not need to get it from input filename or header
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
bool   gCalcEarthXYZ = false; // for the MWA we require conversion from local (x,y,z) to Earth's (X,Y,Z) coordinates , for EDA2 we can keep using local (x,y,z)
                              // TODO : in the future revise if this parameter is required or we can tell based on the data and other parameter (-Z option or -X ???)

// input files :
string gAntennaPositionsFile="antenna_locations.txt";
string gFlaggedAntennasListString;
vector<int> gFlaggedAntennasList;

// output files:
string gOutputDirectory="./";
string gBaseOutFitsName;

// statistics:
string gStatPositionRadius;

CImagerParameters gImagerParameters;

double parse_input_filename( const char* szInBaseName )
{
   char szDateTime[64],szRest[128];
   int channel;
   double unixtime = get_dttm();
   
   // chan_204_20230205T034618_simul_XX
   int ret = sscanf(in_basename.c_str(),"chan_%d_%s_%s",&channel,szDateTime,szRest);
   if( ret >= 2 ){
      char szDate[16],szTime[16],szDateTimeUTC[64];;
      memset(szDate,'\0',16);
      memset(szTime,'\0',16);
      strncpy(szDate,szDateTime,8);
      strncpy(szTime,szDateTime+8+1,6);
      sprintf(szDateTimeUTC,"%s_%s",szDate,szTime);
            
      unixtime = get_gmtime_from_string( szDateTimeUTC );      
      PRINTF_DEBUG("DEBUG : successfully parsed input file name %s into %d/%s -> %s UTC -> unixtime = %.8f\n",in_basename.c_str(),channel,szDateTime,szDateTimeUTC,unixtime);
   }else{
      PRINTF_WARNING("WARNING : could not parse file name %s using format chan_%%s_%%s_%%s -> only %d strings parsed out\n",in_basename.c_str(),ret);
   }
   
   return unixtime;
}

void usage()
{
   double default_pixsize_deg = gImagerParameters.m_PixsizeInRadians*(180.00/M_PI);

   printf("pacer_dirty_image VISIBILITY_FITS_BASENAME\n\n\n");
   
   printf("Exected input visibility FITS files are VISIBILITY_FITS_BASENAME_vis_realPOSTFIX.fits and VISIBILITY_FITS_BASENAME_vis_imagPOSTFIX.fits , where POSTFIX is specified by option -p\n");
   
   printf("\t-O OUTDIR  : full path to output directory where FITS files are saved [default %s]\n",gOutputDirectory.c_str());
   printf("\t-J outfile base : base name for the output FITS files [default %s]\n",gBaseOutFitsName.c_str());
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
   printf("\t-N N_TIMES : create the same image N times for profiling or other tests [default %d]\n",gProfilesNImages);
   printf("\t-P POLARISATION_TO_IMAGE : XX, YY or I [default %s]\n",PolarizationToString(gImagerParameters.pol_to_image));
   printf("\t-v VERBOSITY/DEBUG level [default %d]\n",CPacerImager::m_ImagerDebugLevel);
   printf("\t-V FILE_SAVE_LEVEL : to control number of FITS files saved [default %d]\n",CPacerImager::m_SaveFilesLevel);
   printf("\t-s : print statistics of the final sky image (mean,median,rms,rms_iqr etc)\n");
   printf("\t-S Xc,Yc,R : specify circular window of radius R around (Xc,Yc) position where additional statistics (as -s) are calculated\n");   
   printf("\t-c CAL_SOL_FILE : specify file with calibration solutions [default None], currently supported formats : .txt\n");
   printf("\t-I : include auto-correlations [not included by default]\n");
   printf("\t-M META_DATA_FILE : name of meta data file\n");
   printf("\t-E : recalculate antenna locations (E,N,Z) in the config/metafits files to Earth's coordinates (X,Y,Z) [default disabled]\n");
   printf("\t-X : test option to use antenna positions in XYZ (WG54) format rather than local coordinates (default disabled, i.e. local XYZ), for testing only\n");      
   printf("\t-U unixtime : unix time of the start of the data [default none -> use current time]\n");
   printf("\t-L : apply cable correction based on cable lengths in the .metafits file [default disabled]\n");
   printf("\t-K CONVOLITION_KERNEL_SIZE : size of convolution kernel in pixels [default %d], <=0 -> no convolution kernel (no anti-aliasing)\n",gImagerParameters.m_nConvolvingKernelSize);
   printf("\t-o OPTION  : additional options specified by name, such as:\n");
   printf("\t             miriad : uses the exact settings of MIRIAD natural weighting for validation\n");
   printf("\t             min_w  : minimum limint on W\n");
   printf("\t             max_w  : maximum limint on W\n");
   printf("\t             pixsize : pixel size in degrees [default = %.3f [deg]]\n",default_pixsize_deg);

   exit(0);
}

void parse_cmdline(int argc, char * argv[] ) {
   char optstring[] = "hp:n:g:r:i:f:F:w:m:a:ZP:v:V:A:sS:o:O:c:IM:XEU:J:LK:N:";
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

         case 'E':
            gCalcEarthXYZ = true;
            break;

         case 'c':
//            gCalibrationSolutions.m_filename = optarg;
            break;

         case 'f':
            gFrequencyMHz = atof( optarg );
            break;

         case 'F':
            gFOV_degrees = atof( optarg );
            break;

         case 'I':
//            imager.m_bIncludeAutos = true;
            break;
            
         case 'J':
            if( optarg && strlen(optarg) ){
               gBaseOutFitsName = optarg;
            }
            break;   

         case 'K':
            if( optarg && strlen(optarg) ){
               gImagerParameters.m_nConvolvingKernelSize = atol(optarg);
            }
            break;

         case 'L':
            gImagerParameters.m_bApplyCableCorr = true;
            break;
            
         case 'M':
            gImagerParameters.m_MetaDataFile = optarg;
            gImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
            gZenithImage = false;
            break;

         case 'p':
            gPostfix = optarg;
            break;

         case 'P': // parse polarisation
            if( optarg ){
               switch( ::toupper(optarg[0]) ){
                  case 'X' :
                     gImagerParameters.pol_to_image = Polarization::XX;
                     break;
                  case 'Y' :
                     gImagerParameters.pol_to_image = Polarization::YY;
                     break;
                  case 'I' :
                     gImagerParameters.pol_to_image = Polarization::I;
                     break;
                     
                  default :
                     gImagerParameters.pol_to_image = Polarization::XX;   
                     
               }
            }
            break;
            

         case 'N':
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
               // printf("DEBUG0 : value = %s\n",optarg);
/*               if( strcmp(optarg,"miriad")==0 || strcmp(optarg,"MIRIAD")==0 ){
                  CPacerImager::m_bCompareToMiriad = true;
               }*/
               if( strncmp(optarg,"min_w",5)==0 || strncmp(optarg,"MIN_W",5)==0 ){
                  // printf("DEBUG : value = %s\n",optarg+6);
                  gImagerParameters.m_MinW = atof(optarg+6);
               }
               if( strncmp(optarg,"max_w",5)==0 || strncmp(optarg,"MAX_W",5)==0 ){
                  gImagerParameters.m_MaxW = atof(optarg+6);
               }
               if( strncmp(optarg,"pixsize",7)==0 || strncmp(optarg,"PIXSIZE",7)==0 ){
                  gImagerParameters.m_PixsizeInRadians = atof(optarg+8)*(M_PI/180.00);
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

         case 'X' :
            CImagerParameters::m_bAntennaPositionsXYZ = true;
            break;
            
         case 'U' :
            if( optarg && strlen(optarg) ){
               gImagerParameters.m_fUnixTime = atof( optarg );
               gUnixTimeParamProvided = true;
            }
            break;
 
         case 'Z' :
            gZenithImage = true;
            gImagerParameters.m_bConstantUVW = true; // constant UVW (not re-calculated as a function of time due to pointing)
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
   
   if( !gUnixTimeParamProvided ){
      gImagerParameters.m_fUnixTime = parse_input_filename( in_basename.c_str() );
   }
}

void print_parameters()
{
    double pixsize_deg = gImagerParameters.m_PixsizeInRadians*(180.00/M_PI);
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");
    printf("Base name for FITS  = %s\n",in_basename.c_str());
    printf("Output dir   = %s\n",gOutputDirectory.c_str());
    printf("Output FITS base name = %s\n",gBaseOutFitsName.c_str());
    printf("Postfix      = %s\n",gPostfix.c_str());
    printf("Image size   = (%d x %d)\n",gImageSize,gImageSize);
    printf("Do gridding  = %d\n",gDoGridding);    
    if( !gDoGridding ){
       printf("External REAL/IMAG FITS FILES : %s / %s\n",g_in_fits_file_uv_re.c_str(),g_in_fits_file_uv_im.c_str());
    }
    printf("Frequency           = %.4f [MHz]\n",gFrequencyMHz);
    printf("Polarisation        = %s\n",PolarizationToString(gImagerParameters.pol_to_image));
    printf("Field of view (FoV) = %.4f [deg]\n",gFOV_degrees);
    printf("Pixel angular size  = %.8f [deg]\n",pixsize_deg);
    printf("Weighting           = %s\n",gUniformWeighting.c_str());
    printf("UV range            = %.4f - Infinity\n",gMinUV);   
    printf("W range             = %e - %e\n",gImagerParameters.m_MinW,gImagerParameters.m_MaxW);
    printf("Antenna positions file = %s\n",gAntennaPositionsFile.c_str());
    printf("Meta data file      = %s\n",gImagerParameters.m_MetaDataFile.c_str());
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
//    printf("MIRIAD compatibility = %d\n",CPacerImager::m_bCompareToMiriad);
//    printf("Calibration solutions file = %s\n",gCalibrationSolutions.m_filename.c_str());
//    printf("Include auto-correlations  = %d\n",imager.m_bIncludeAutos);
    printf("Antenna positions XYZ      = %d\n",gImagerParameters.m_bAntennaPositionsXYZ);
    printf("Data start time            = %.8f\n",gImagerParameters.m_fUnixTime);   
    printf("Cable correction           = %d\n",gImagerParameters.m_bApplyCableCorr);
    printf("Size of convolution kernel = %d (<=0 anti-aliasing convolution disabled)\n",gImagerParameters.m_nConvolvingKernelSize);
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

// parse and print parameters :
  parse_cmdline( argc , argv );
  print_parameters();


#ifdef IMAGER_HIP
  CPacerImagerHip imager;   
#else
  printf("DEBUG0 : pol_to_image = %d\n", int(gImagerParameters.pol_to_image) );
  CPacerImager imager( gImagerParameters.m_MetaDataFile.c_str(), gFlaggedAntennasList, false, gImagerParameters.pol_to_image ) ;
#endif  
  
  gImagerParameters.SetGlobalParameters( gAntennaPositionsFile.c_str(), (gZenithImage) ); // Constant UVW when zenith image (-Z)
  gImagerParameters.m_szOutputDirectory = gOutputDirectory.c_str();
  gImagerParameters.m_bCalcEarthXYZ = gCalcEarthXYZ;
  gImagerParameters.m_ImageSize = gImageSize;
  gImagerParameters.m_fCenterFrequencyMHz = gFrequencyMHz;
  gImagerParameters.m_ImageFOV_degrees = gFOV_degrees;
  gImagerParameters.m_fMinUV = gMinUV;
  gImagerParameters.m_szWeighting = gUniformWeighting.c_str();  
  
  Initialise( imager );
  imager.apply_geom_correction = false;
  imager.apply_cable_correction = false;
  
  // set antenna flags (if there are):
//  if( gFlaggedAntennasList.size() > 0 ){
//     imager.SetFlaggedAntennas( gFlaggedAntennasList );
//  }

  // this loop is just for testing purposes , if -P 100 specified it means the same image will be generated 100 times so 
  // that testing some initialisations, UVW relculation or profiling can be done
  
  CBgFits re_xx,im_xx;
  read_corr_matrix( in_basename.c_str(), re_xx, im_xx, gPostfix.c_str(), 0 );
  
  CBgFits re_yy, im_yy;
  // try to read YY correlation matrix if possible
  bool bYY=read_corr_matrix( in_basename.c_str(), re_yy, im_yy, gPostfix.c_str(), 3 );
    
  ObservationInfo obsInfo;
  obsInfo = VCS_OBSERVATION_INFO;
/*  obsInfo.nAntennas = re.GetXSize();
  obsInfo.nFrequencies = 1;
  obsInfo.nPolarizations = 1;
  obsInfo.nTimesteps = 1;
  obsInfo.timeResolution = 1;
  obsInfo.frequencyResolution = 0.040;*/
  // TODO : fill the rest 
  
  
  Visibilities* xcorr = nullptr;
  if( bYY ){
     xcorr = ConvertFits2XCorr( re_xx, im_xx, obsInfo, &re_yy, &im_yy );
  }else{
     xcorr = ConvertFits2XCorr( re_xx, im_xx, obsInfo );
  }
//  int vis_size = re.GetXSize()*re.GetYSize();
//  MemoryBuffer<std::complex<float>> vis_buffer( vis_size );
//  Visibilities xcorr( std::move(vis_buffer), obsInfo, 1, 1 );
// how to use this constructor to initialise with some specific correlation matrix ???
//  Visibilities(MemoryBuffer<std::complex<float>>&& data, const ObservationInfo& obsInfo, unsigned int nIntegrationSteps, unsigned int nAveragedChannels) 
// Intialise with data in re and im 

  std::string outfits{ "test.fits" };
  PACER_PROFILER_START  
  for(int i=0;i<gProfilesNImages;i++){  
     printf("\n\nProcessing image : %d\n",i);
     
     // TODO : 
     // convert input data to xcorr structure - something I've had code for ready to go !!!
     
     imager.m_nConvolvingKernelSize = gImagerParameters.m_nConvolvingKernelSize;
     imager.m_fFrequencyMHz = gFrequencyMHz;
     Images image = imager.run_imager( (*xcorr), gImageSize, -1000, "" );
     image.to_fits_files( gOutputDirectory.c_str() );
  
/*       imager.run_imager( in_basename.c_str(),  gPostfix.c_str(), 
                          gDoGridding, gDo_dirty_image, 
                          g_in_fits_file_uv_re.c_str(), g_in_fits_file_uv_im.c_str(),
                          gBaseOutFitsName.c_str() );                      */
  }                   
  
  delete xcorr;
  
  PACER_PROFILER_END("Full execution of imager (including I/O) took ")  
}


// imager parameters :
// CImagerParameters CPacerImager::m_ImagerParameters;
std::string CImagerParameters::m_AntennaPositionsFile; // antenna positions are same for all freq. channels -> static, similarly m_bConstantUVW
std::string CImagerParameters::m_MetaDataFile;         // meta data file (.txt or .metafits)
bool   CImagerParameters::m_bAutoFixMetaData=false;
bool   CImagerParameters::m_bConstantUVW=false; // default false and only true to zenith phase centered images (a
bool   CImagerParameters::m_bAntennaPositionsXYZ=false; 
bool   CImagerParameters::m_bCalcEarthXYZ=false; // for the MWA we require conversion from local (x,y,z) to Earth's (X,Y,Z) coordinates , for EDA2 we can keep using local (x,y,z)
bool   CImagerParameters::m_bApplyGeomCorr=false;
bool   CImagerParameters::m_bApplyCableCorr=false;


CImagerParameters::CImagerParameters()
  : m_fUnixTime(0), m_fCenterFrequencyMHz(0), m_fBandwidthMHz(0), m_ImageSize(0), m_fMinUV(-1000), m_ImageFOV_degrees(180.00), m_integrationTime(1.00), m_nConvolvingKernelSize(-1)
//  : m_bConstantUVW(false)
{
   m_PixsizeInRadians = ((230.00/300.00)/(2.00*35.00)); // Lambda/B_max at freq = 230 MHz -> Lambda = (230/300) m ~= 1.304347 m
   m_szOutputDirectory = "./";   
   m_MinW = -std::numeric_limits<double>::max();
   m_MaxW = std::numeric_limits<double>::max();
}

CImagerParameters::CImagerParameters( const CImagerParameters& right )
{
   (*this) = right;
}

CImagerParameters& CImagerParameters::operator=( const CImagerParameters& right )
{
//   m_AntennaPositionsFile = right.m_AntennaPositionsFile;
//   m_bConstantUVW = right.m_bConstantUVW;
    m_szOutputDirectory = right.m_szOutputDirectory;
    m_fUnixTime = right.m_fUnixTime;
    m_fCenterFrequencyMHz = right.m_fCenterFrequencyMHz;
    m_fBandwidthMHz = right.m_fBandwidthMHz;
    m_ImageSize = right.m_ImageSize;
    m_fMinUV = right.m_fMinUV;
    m_ImageFOV_degrees = right.m_ImageFOV_degrees;
    m_MinW = right.m_MinW;
    m_MaxW = right.m_MaxW;
    m_integrationTime = right.m_integrationTime;
    m_nConvolvingKernelSize = right.m_nConvolvingKernelSize;
    
    return (*this);
}


// WARNING : this function may cause troubles for multithreaded situation - 
// TODO : protect with a mutex if necessary 
void CImagerParameters::SetGlobalParameters( const char* szAntennaPositionsFile, bool bConstantUVW )
{
   if( szAntennaPositionsFile && strlen( szAntennaPositionsFile ) ){
      m_AntennaPositionsFile = szAntennaPositionsFile;
   }
   
   m_bConstantUVW = bConstantUVW;
}


bool gInitialised = false;

void Initialise( CPacerImager& imager )
{
  
  if( !gInitialised )
  {
     gInitialised = true;
     
    // test file with antenna positions can be used to overwrite whatever was in .metafits
    if( strlen( gImagerParameters.m_AntennaPositionsFile.c_str() ) && MyFile::DoesFileExist(  gImagerParameters.m_AntennaPositionsFile.c_str() ) )
    {
       bool bConvertToXYZ = false;
       if( !gImagerParameters.m_bConstantUVW )
       { // if non-constant UVW -> non zenith phase centered all-sky image
          bConvertToXYZ = true;
       }
       if( gImagerParameters.m_bAntennaPositionsXYZ )
       { // text file already has XYZ in WG54 system - for example CASA dump 
          printf("INFO : antenna positions already in XYZ coordinate system (WG54) no need to convert\n");
          bConvertToXYZ = false;
       }
              
       printf("DEBUG0 :  m_bAntennaPositionsXYZ = %d , bConvertToXYZ = %d\n",gImagerParameters.m_bAntennaPositionsXYZ,bConvertToXYZ);
       // read antenna positions and do whatever else is necessary (update flags etc)
       // ReadAntennaPositions( bConvertToXYZ );
       int m_nAntennas = imager.m_MetaData.m_AntennaPositions.ReadAntennaPositions( gImagerParameters.m_AntennaPositionsFile.c_str(), bConvertToXYZ  );
       printf("DEBUG0 : read %d antenna positions , bConvertToXYZ = %d\n",m_nAntennas,bConvertToXYZ);
       imager.UpdateFlags(); // if antenna positions are read now for the first time, the flags need to be updated

       
       if( /*true ||*/ strlen( gImagerParameters.m_MetaDataFile.c_str() ) == 0 )
       { // only calculate UVW here when Metadata is not required
          // initial recalculation of UVW at zenith (no metadata provided -> zenith):       
          // WARNING : bool CPacerImager::CalculateUVW() - could be used, but it also calls this function itself which may cause infinite recursise call
          // m_Baselines = m_MetaData.m_AntennaPositions.CalculateUVW( m_U, m_V, m_W, (CPacerImager::m_SaveFilesLevel>=SAVE_FILES_DEBUG), gImagerParameters.m_szOutputDirectory.c_str(), m_bIncludeAutos );
          // UpdateParameters();
          imager.CalculateUVW(); // bForce=true to force call of m_MetaData.m_AntennaPositions.CalculateUVW and , bInitialise=false to avoid call to this (CPacerImager::Initialise) function in a recursive way !
          PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines (include Autos = %d)\n",imager.m_Baselines,imager.m_bIncludeAutos);
       }else
       {
          printf("INFO : non-zenith pointing meta data is required to calculate UVW\n");
       }
    }else
    {
       PRINTF_WARNING("WARNING : antenna position file %s not specified or does not exist\n",gImagerParameters.m_AntennaPositionsFile.c_str());
    }
    
    // read all information from metadata 
    if( strlen( gImagerParameters.m_MetaDataFile.c_str() ) && MyFile::DoesFileExist( gImagerParameters.m_MetaDataFile.c_str() ) ){
      PRINTF_INFO("INFO : reading meta data from file %s\n",gImagerParameters.m_MetaDataFile.c_str());
      
      double obsid = -1;
      if( CImagerParameters::m_bAutoFixMetaData ){
         obsid = CObsMetadata::ux2gps( gImagerParameters.m_fUnixTime );
      }
//      if( !imager.m_MetaData.ReadMetaData( gImagerParameters.m_MetaDataFile.c_str(), obsid, gImagerParameters.m_integrationTime ) ){
      if( !imager.m_MetaData.ReadMetaData( gImagerParameters.m_MetaDataFile.c_str() ) ){
         PRINTF_ERROR("ERROR : could not read meta data from file %s\n",gImagerParameters.m_MetaDataFile.c_str() );
      }
      if( obsid > 0 ){
          imager.m_MetaData.fix_metafits( obsid, gImagerParameters.m_integrationTime ); 
      }
    }       
    
    
  }
}

Visibilities* ConvertFits2XCorr( CBgFits& vis_re, CBgFits& vis_im, ObservationInfo& obsInfo, CBgFits* vis_re_yy, CBgFits* vis_im_yy )
{
//   printf("ERROR : function CPacerImager::ConvertFits2XCorr not implemented yet !\n");   
//   exit(-1); 
   
   // TODO :
   // xcorr - needs to be allocated so that xcorr.data pointer is not NULL !
   
   // TODO :
   // use : std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, i, j );
   // but other way around 
   int n_ant = vis_re.GetXSize(); // xcorr.obsInfo.nAntennas;
   int n_corrs = 4; // 4 correlation products : XX XY YX YY 
   int n_baselines = n_ant*(n_ant+1)/2;
   
   int time_step = 0;
   int fine_channel = 0;
   
   obsInfo.nAntennas = n_ant;
   obsInfo.nFrequencies = 1;
   obsInfo.nPolarizations = 2;
   obsInfo.nTimesteps = 1;
   obsInfo.timeResolution = 1.0;
   obsInfo.frequencyResolution = 0.040;
   obsInfo.coarseChannelBandwidth = 1.28;
   obsInfo.coarseChannel = 109;
   obsInfo.startTime = 1534699078;
   obsInfo.coarse_channel_index = 0;
   obsInfo.metadata_file = "1218734296.metafits";
      
   const size_t matrixSize = n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations;
   const size_t nIntervals  = (obsInfo.nTimesteps); // TODO  + voltages.nIntegrationSteps - 1) / voltages.nIntegrationSteps;
   unsigned int nChannelsToAvg = 1; // TODO : verify
   const size_t nOutFrequencies = obsInfo.nFrequencies / nChannelsToAvg;
   const size_t nValuesInTimeInterval = matrixSize * nOutFrequencies;
   const size_t outSize = nValuesInTimeInterval * nIntervals;
   unsigned int avgCh = fine_channel / nChannelsToAvg;
   
//   std::complex<VISIBILITY_TYPE>* data = new std::complex<VISIBILITY_TYPE>[matrixSize]; // was n_ant*n_ant*obsInfo.nPolarizations
//   Visibilities* xcorr = new Visibilities( data, obsInfo, 1, 1 );
   int vis_size = vis_re.GetXSize()*vis_re.GetYSize();
   MemoryBuffer<std::complex<float>> vis_buffer( vis_size*4 ); // 4 correlations (2pol x 2pol)
   Visibilities* xcorr = new Visibilities( std::move(vis_buffer), obsInfo, 1, 1 );

   xcorr->nIntegrationSteps = 1;
   xcorr->nAveragedChannels = 1;
   xcorr->nFrequencies = 1;
   printf("DEBUG : CPacerImager::ConvertFits2XCorr - allocations ok, visibility array size = %d\n",int(matrixSize));fflush(stdout);

   for(int i=0;i<n_ant;i++){ // loop over ant1          
     for(int j=0;j<=i;j++){ // loop over ant2 
        std::complex<float>* vis_xx = xcorr->at( time_step, fine_channel, i, j );
        std::complex<float>* vis_yy = xcorr->at( time_step, fine_channel, i, j ) + 3;
        
        if( vis_xx ){        
           double re = vis_re.getXY( j , i );
           double im = vis_im.getXY( j , i );

           std::complex<float> vis_value( re, im );
           (*vis_xx) = vis_value;
           
           if( vis_re_yy && vis_im_yy ){
              double re_yy = vis_re_yy->getXY( j , i );
              double im_yy = vis_im_yy->getXY( j , i );
              std::complex<float> vis_value_yy( re_yy, im_yy );
              (*vis_yy) = vis_value_yy;
           }else{
              (*vis_yy) = std::complex<float>(0.0,0.0); // if no YY correlation matrix provided set YY visibilities to 0
              // Can set 2x vis_xx for testing 
              // (*vis_yy) = vis_value*std::complex<float>(2.00,0); // just to test if I get image YY = 2x image in XX 
           }
        }else{
           printf("ERROR in code : no visibility at (%d,%d,%d,%d)\n",time_step, fine_channel, i, j );
           exit(-1);
        }
  
// not required as only one side of the matrix is stored ?       
//        std::complex<VISIBILITY_TYPE>* vis_conj = xcorr->at( time_step, fine_channel, j , i );
//        (*vis_conj) = conj( vis_value );
     }
   }
   printf("DEBUG : CPacerImager::ConvertFits2XCorr end\n");

   return xcorr;
}

bool read_corr_matrix( const char* basename, CBgFits& fits_vis_real, CBgFits& fits_vis_imag, const char* szPostfix, int pol )
{
  // ensures initalisation of object structures 
//  Initialise();

  // creating FITS file names for REAL, IMAG and U,V,W input FITS files :
  string fits_file_real = basename;
  fits_file_real += "_vis_real";
  if( strlen( szPostfix ) ){
     fits_file_real += szPostfix;
  }
  char szFitsFile[256];
  sprintf(szFitsFile,"%s_pol%d.fits",fits_file_real.c_str(),pol);
  fits_file_real = szFitsFile;

  string fits_file_imag = basename;
  fits_file_imag += "_vis_imag";
  if( strlen( szPostfix ) ){
     fits_file_imag += szPostfix;
  }
  sprintf(szFitsFile,"%s_pol%d.fits",fits_file_imag.c_str(),pol);
  fits_file_imag = szFitsFile;


  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_INFO_LEVEL){
     printf("Expecting the following files to exist:\n");
     printf("\t%s\n",fits_file_real.c_str()); 
     printf("\t%s\n",fits_file_imag.c_str()); 
  }
  
  // REAL(VIS)
  PRINTF_INFO("Reading fits file %s ...\n",fits_file_real.c_str());
  if( fits_vis_real.ReadFits( fits_file_real.c_str(), 0, 1, 1 ) ){
     // try just with _real.fits (not _vis_real.fits) :
     printf("INFO : could not read visibility FITS file %s\n",fits_file_real.c_str());
     fits_file_real = basename;
     fits_file_real += "_real.fits";
     if( fits_vis_real.ReadFits( fits_file_real.c_str(), 0, 1, 1 ) ){
        printf("ERROR : could not read visibility FITS file %s\n",fits_file_real.c_str());
        return false;
     }
  }else{
     PRINTF_INFO("OK : fits file %s read ok\n",fits_file_real.c_str());
  }

  // IMAG(VIS)
  PRINTF_INFO("Reading fits file %s ...\n",fits_file_imag.c_str());
  if( fits_vis_imag.ReadFits( fits_file_imag.c_str(), 0, 1, 1 ) ){
     // try just with _imag.fits (not _vis_imag.fits) :
     printf("INFO : could not read visibility FITS file %s\n",fits_file_imag.c_str());      
     fits_file_imag = basename;
     fits_file_imag += "_imag.fits";
     if( fits_vis_imag.ReadFits( fits_file_imag.c_str(), 0, 1, 1 ) ){  
        printf("ERROR : could not read visibility FITS file %s\n",fits_file_imag.c_str());
        return false;
     }
  }else{
     PRINTF_INFO("OK : fits file %s read ok\n",fits_file_imag.c_str());
  }

//  bool ret = ReadOrCalcUVW( basename , szPostfix );
  return true;  
}
