//------------------------------------------------------------------------------------------------------------------------------------------------
//Multi-image imager - only for GPU (no CPU version)
//------------------------------------------------------------------------------------------------------------------------------------------------
// TODO :
// - image both polarisations !!!!
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <complex.h>

#include <bg_globals.h>
#include <bg_fits.h>
#include <mystring.h>
#include <myparser.h>
#include <mystrtable.h>

#include <vector>
using namespace std;

// L-file reader : https://github.com/marcinsokolowski/lfile2corrmatrix 
#include <lfile_reader.h>

// defines :
#include "../src/realtime_image_analyser.h"
#include "../src/pacer_imager_defs.h"
#include "../src/pacer_imager_parameters.h"
#include "../src/pacer_common.h"

#ifdef IMAGER_HIP
// Pacer imager class :
#include "../src/hip/pacer_imager_multi_hip.h"
#else
#include "../src/pacer_imager.h"
#endif

#ifdef _PACER_PROFILER_ON_
#include <mydate.h>
#endif

#ifdef IMAGER_HIP
#include <gpu_fft.hpp>
#endif

string in_basename="1276619416_20200619164456";
string gPostfix="";
string out_fits="1276619416_20200619164456_dirty_image.fits";
string gOutDir="fits_images";
// string out_rms_fits="out_rms.fits";

// REAL and IMAG external UV GRID files to make dirty image:
string g_in_fits_file_uv_re;
string g_in_fits_file_uv_im;

// options to enable particular functions:
bool gDoGridding=true;
bool gDo_dirty_image=true;

int gImageSize = 180; // was 512
bool gGPUImaging=false; // also has to be compiled with GPU/HIP (USE_HIP=ON enabled)

// weighting :
string gUniformWeighting="N"; // uniform, means that sum of visibilities in a given UV cell is divided by a number of points added

double gFrequencyMHz = 159.375; // default of the station in channel 204 is 204*(400.00/512.00) = 159.375 MHz 
double delta_freq_mhz = (400.00/512.00)*(32.00/27.00);
double gFOV_degrees  = 180.00;
double gMinUV        = -1000;
bool   gZenithImage  = true; // yes for EDA2 
int    gProfilesNImages = 1;
bool   gCalcEarthXYZ = false; // for the MWA we require conversion from local (x,y,z) to Earth's (X,Y,Z) coordinates , for EDA2 we can keep using local (x,y,z)
                              // TODO : in the future revise if this parameter is required or we can tell based on the data and other parameter (-Z option or -X ???)

double gIntegrationTimeInSec=0.10001664; // default 100ms or should be 100.01664 for 32channels x 2894 averages x 1.08 usec = 100016.64 usec = 100.01664 msec 

// input files :
string gAntennaPositionsFile="antenna_locations.txt";
string gFlaggedAntennasListString;
vector<int> gFlaggedAntennasList;

// input visibilities (L-files):
int n_ch=32;
int n_antennas=256; // mostly for CPU version as forGPU version the number is obtained from parameter file 
int n_timesteps=-1; // unknown number of timestamps
string gMappingFile = "instr_config.txt";
int gTimestampsToProcess=-1; // number of timestamps to process (set to >0 to process smaller number than all of them)
bool gDoImaging=true;
int  gImageChannel=0;
bool gWriteCorrMatrixAll=false;
bool gIncludeAutos=false;
string gMetaDataFile;
bool gConstantUVW=false;
double gUnixTime=0.00;

// calibration :
string gCalibrationSolutionsFilename;

// output files:
int   gSaveFitsFilesFreq=0; // specifies frequency of saving output FITS files, 0 - nothing saved, 1 - all saved, >1 -> every N-th image saved
int   gSaveImaginaryFITS=0; // specifies frequency of saving output FITS files of imaginary part, 0 - nothing saved, 1 - all saved, >1 -> every N-th image saved
string gOutputDirectory="./";

// statistics:
string gStatPositionRadius;

// options of CUDA gridding etc :
bool gUseGridBlocks=false;

void usage()
{
   printf("pacer_imager_multi LFILE_BASENAME\n\n\n");
   
   printf("Exected input visibility FITS files are VISIBILITY_FITS_BASENAME_vis_realPOSTFIX.fits and VISIBILITY_FITS_BASENAME_vis_imagPOSTFIX.fits , where POSTFIX is specified by option -p\n");

   printf("\t-C N_CH : number of frequency channels [default %d]\n",n_ch);
   printf("\t-j IMAGE_CHANNEL : image a specific frequency channel [default %d], <0 - means image all channels\n",gImageChannel);
   printf("\t-b : use grid blocks to separate into images of channels [default %d]\n",gUseGridBlocks);
   printf("\t-B : save resulting sky images to FITS files every N-th image [default %d], 0 - nothing saved, 1 - all saved, >1 -> every N-th image saved\n",gSaveFitsFilesFreq);      
   printf("\t-G : execute GPU version (default %d , 0 - is CPU)\n",gGPUImaging);   
   printf("\t-O OUTDIR  : full path to output directory where FITS files are saved [default %s]\n",gOutputDirectory.c_str());
   printf("\t-N NUMBER_OF_TIMESTAMPS : specify maximum number of timestamps to image [default %d]\n",gTimestampsToProcess);   
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
   printf("\t-n IMAGE_SIZE : single value N for image size N x N pixels [default %d]. WARNING : only even image sizes are allowed by this program\n",gImageSize);
   printf("\t-Z : image phase centered at zenith, re-calculation of UVW is not required\n");
   printf("\t-P N_TIMES : create the same image N times for profiling or other tests [default %d]\n",gProfilesNImages);
   printf("\t-v VERBOSITY/DEBUG level [default %d]\n",CPacerImager::m_ImagerDebugLevel);
   printf("\t-V FILE_SAVE_LEVEL : to control number of FITS files saved [default %d]\n",CPacerImager::m_SaveFilesLevel);
   printf("\t-s : print statistics of the final sky image (mean,median,rms,rms_iqr etc)\n");
   printf("\t-S Xc,Yc,R : specify circular window of radius R around (Xc,Yc) position where additional statistics (as -s) are calculated\n");   
   printf("\t-c CAL_SOL_FILE : specify file with calibration solutions [default None], currently supported formats : .txt\n");
   printf("\t-I : include auto-correlations [not included by default]\n");
   printf("\t-M META_DATA_FILE : name of meta data file\n");
   printf("\t-E : recalculate antenna locations (E,N,Z) in the config/metafits files to Earth's coordinates (X,Y,Z) [default disabled]\n");
   printf("\t-X : test option to use antenna positions in XYZ (WG54) format rather than local coordinates (default disabled, i.e. local XYZ), for testing only\n");      
   printf("\t-t INTTIME : integration time in seconds [default %.6f seconds = %.6f ms]\n",gIntegrationTimeInSec,gIntegrationTimeInSec*1000.00);
   printf("\t-U unixtime : unix time of the start of the data [default none -> use current time]\n");
   printf("\t-o OPTION  : additional options specified by name, such as:\n");
   printf("\t             miriad : uses the exact settings of MIRIAD natural weighting for validation\n");
   printf("\t             save_imag : save imaginary files [default %d]\n",gSaveImaginaryFITS);

   exit(0);
}

void parse_cmdline(int argc, char * argv[] ) {
   char optstring[] = "hp:n:g:r:i:f:F:w:m:a:ZP:v:V:A:sS:o:O:c:IM:XEU:t:N:GB:bC:j:";
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

         case 'b':
            gUseGridBlocks = true;
            break;

         case 'B':
            if( optarg && strlen(optarg) ){
               gSaveFitsFilesFreq = atol( optarg );
            }             
            break;

         case 'C':
            if( optarg && strlen(optarg) ){
              n_ch = atol( optarg );
            }
            break;

         case 'j':
            if( optarg && strlen(optarg) ){
              gImageChannel = atol( optarg );
            }
            break;

         case 'E':
            gCalcEarthXYZ = true;
            break;

         case 'c':
            gCalibrationSolutionsFilename = optarg;
            break;

         case 'f':
            gFrequencyMHz = atof( optarg );
            break;

         case 'F':
            gFOV_degrees = atof( optarg );
            break;

         case 'G':
            gGPUImaging = true;
            break;

         case 'I':
            gIncludeAutos = true;
            break;

         case 'N':
            gTimestampsToProcess = atol( optarg );
            break;
            
         case 'M':
            gMetaDataFile = optarg;
            gConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
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
               if( strcmp(optarg,"save_imag")==0 || strcmp(optarg,"SAVE_IMAG")==0 ){
                  gSaveImaginaryFITS = 1;
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
            
         case 't' :
            if( optarg && strlen(optarg) ){
               gIntegrationTimeInSec = atof( optarg );
            }
            break;
 
         case 'U' :
            if( optarg && strlen(optarg) ){
               gUnixTime = atof( optarg );
            }
            break;
 
         case 'Z' :
            gZenithImage = true;
            gConstantUVW = true; // constant UVW (not re-calculated as a function of time due to pointing)
            // imager.m_ImagerParameters.m_bConstantUVW = true; // constant UVW (not re-calculated as a function of time due to pointing)
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
   
   if( (gImageSize % 2) != 0 ){
      printf("ERROR : only even image sizes are allowed to make the code as optimal as possible use -n EVEN_PIXEL_SIZE\n");
      usage();            
   }
}

void print_parameters()
{
    printf("############################################################################################\n");
    printf("PARAMETERS :\n");
    printf("############################################################################################\n");    
    printf("GPU version         = %d\n",gGPUImaging);
    printf("CUDA blocks         = %d\n",gUseGridBlocks);    
    printf("Image channel       = %d\n",gImageChannel);   
    printf("Base name for FITS  = %s\n",in_basename.c_str());
    printf("Save FITS files freq. = %d\n",gSaveFitsFilesFreq);
    printf("Number of frequency channels = %d\n",n_ch);    
    printf("Number of timestamps to process = %d\n",gTimestampsToProcess);
    printf("Integration time    = %.6f [ms]\n",gIntegrationTimeInSec*1000.00);
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
    printf("Meta data file      = %s\n",gMetaDataFile.c_str());
    printf("Zenith image        = %d ( Constant UVW = %d )\n",gZenithImage,gConstantUVW);
    printf("Profile N times the same image = %d\n",gProfilesNImages);
    printf("Verbosity level     = %d\n",CPacerImager::m_ImagerDebugLevel);
    printf("File save level     = %d\n",CPacerImager::m_SaveFilesLevel);
    printf("\tSave imaginary files = %d\n",gSaveImaginaryFITS);        
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
    printf("Calibration solutions file = %s\n",gCalibrationSolutionsFilename.c_str());
    printf("Include auto-correlations  = %d\n",gIncludeAutos);
    printf("Antenna positions XYZ      = %d\n",CImagerParameters::m_bAntennaPositionsXYZ);
    printf("Data start time            = %.8f\n",gUnixTime);   
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

#ifdef IMAGER_HIP
// it's done differently in CPU version :
bool ApplyCalibration( CLfileReader& lfile_reader, CPacerImagerMultiFreqHip* pRealTimeImagerGPU, int pol );
#endif

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

#ifdef IMAGER_HIP
  CPacerImagerMultiFreqHip* pRealTimeImagerGPU = NULL;
#endif  
  CRealTimeImageAnalyser*   pRealTimeImagerCPU = NULL; // CPU imaging object it's either this one of CPacerImagerMultiFreqHip 
  
  
  // parse and print parameters :
  parse_cmdline( argc , argv );
  print_parameters();

  // this loop is just for testing purposes , if -P 100 specified it means the same image will be generated 100 times so 
  // that testing some initialisations, UVW relculation or profiling can be done
  PACER_PROFILER_START
  
  CLfileReader lfile_reader( in_basename.c_str(), n_antennas, n_ch, n_timesteps, gMappingFile.c_str(), gIntegrationTimeInSec );
  // int CLfileReader::GetCorrMatrix( int timestamp, int channel, CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, int pol )
  lfile_reader.Open();
  
  int corr_size = n_antennas*n_antennas;
  int image_size = gImageSize*gImageSize;
  
  // CPU IMAGER VERSION :
  // set parameters for imaging (even if not done):
  if( gGPUImaging ){     
     // initialisation of GPU imager:
#ifdef IMAGER_HIP
     pRealTimeImagerGPU = new CPacerImagerMultiFreqHip();

     pRealTimeImagerGPU->m_ImagerParameters.SetGlobalParameters( gAntennaPositionsFile.c_str(), (gZenithImage) ); // Constant UVW when zenith image (-Z)
     pRealTimeImagerGPU->m_ImagerParameters.m_szOutputDirectory = gOutputDirectory.c_str();
     pRealTimeImagerGPU->m_ImagerParameters.m_bCalcEarthXYZ = gCalcEarthXYZ;
     pRealTimeImagerGPU->m_CalibrationSolutions.m_filename = gCalibrationSolutionsFilename;
     pRealTimeImagerGPU->m_bIncludeAutos = gIncludeAutos;   
     pRealTimeImagerGPU->m_ImagerParameters.m_MetaDataFile = gMetaDataFile.c_str();
     pRealTimeImagerGPU->m_ImagerParameters.m_bConstantUVW = gConstantUVW; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
     pRealTimeImagerGPU->m_ImagerParameters.m_fUnixTime = gUnixTime;
     pRealTimeImagerGPU->m_ImagerParameters.m_fCenterFrequencyMHz = gFrequencyMHz;
     pRealTimeImagerGPU->m_ImagerParameters.m_fBandwidthMHz = delta_freq_mhz;
     pRealTimeImagerGPU->m_ImagerParameters.m_ImageSize = gImageSize;
     pRealTimeImagerGPU->m_ImagerParameters.m_fMinUV = gMinUV;

     printf("DEBUG : n_ch = %d ???\n",n_ch);
     pRealTimeImagerGPU->SetNFreqChannels( n_ch ); // number of blocks by default the same as number of streams   
     if( gUseGridBlocks ){
        pRealTimeImagerGPU->SetGridBlocks( n_ch );
     }

     // TODO : possibly group all these initialisations into one function in class CPacerImagerMultiFreqHip or parant class  
     // only CPU variable :
     pRealTimeImagerGPU->Initialise();
     n_antennas = pRealTimeImagerGPU->m_nAntennas;
     corr_size = n_antennas*n_antennas;
     
     printf("INFO : allocating GPU memory ...\n");fflush(stdout);
     pRealTimeImagerGPU->AllocGPUMemory( corr_size, image_size );
     pRealTimeImagerGPU->AllocInVisBuffer( lfile_reader.n_nBaselines*lfile_reader.m_nChannels, lfile_reader.m_nBaselinesAuto*lfile_reader.m_nChannels );
     printf("INFO : GPU memory allocated OK\n");
     
     if( !pRealTimeImagerGPU->InitialiseAntennaMapping( lfile_reader.m_ArrayConfig.m_CorrelatorMatrix, lfile_reader.m_ArrayConfig.m_Inputs ) ){
        printf("ERROR : initialisation of antenna mappings failed\n");
        exit(-1);
     }
     
     // set antenna flags (if there are):
     if( gFlaggedAntennasList.size() > 0 ){
        pRealTimeImagerGPU->SetFlaggedAntennas( gFlaggedAntennasList );
     }
     
     if( gConstantUVW ){
        // initialise counter only for Constant UVW setup (EDA2/SKA-Low all-sky imaging)
        // copies UVW from m_U, m_V -> u_gpu and v_gpu, calculates counter using GPU kernel
        pRealTimeImagerGPU->InitialiseUVW();
     }
#else
     printf("ERROR : pacer_imager_multi program was compiled in CPU-only version (add -DUSE_HIP=ON in cmake or gpu option for build.sh to enable GPU code)\n");
     exit(-1);
#endif
  }else{
     if( gDoImaging ){
         // this path has already been well tested and produces ok images (but flipped)
         pRealTimeImagerCPU = new CRealTimeImageAnalyser( n_ch, gImageSize, gFrequencyMHz, delta_freq_mhz, gIncludeAutos );
         // set parameters for imaging (even if not done):
         pRealTimeImagerCPU->SetParameters( gAntennaPositionsFile.c_str(), true, gOutDir.c_str(), gCalibrationSolutionsFilename.c_str() );
         
         // set antenna flags (if there are):
         if( gFlaggedAntennasList.size() > 0 ){
            pRealTimeImagerCPU->SetFlaggedAntennas( gFlaggedAntennasList );
         }
     }
  }
  
  CBgFits corr_matrix_real( n_antennas, n_antennas ), corr_matrix_imag( n_antennas, n_antennas );
  // same as in src/pacer_imager.cpp :
  double FoV_radians = gFOV_degrees*M_PI/180.;
  double delta_u = 1.00/(FoV_radians); // should be 2*FoV_radians - see TMS etc 
  double delta_v = 1.00/(FoV_radians); // Rick Perley page 16 : /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
  
   // fetch next timestamp get all correlation matrices and save 
   int timestamp = 0;
   while( lfile_reader.ReadNextTimestamp() && (gTimestampsToProcess<=0 || timestamp<gTimestampsToProcess) ){
      if( lfile_reader.m_pBaselinesAutoBuffer && lfile_reader.m_nLastBytesReadAuto == lfile_reader.m_nBaselineSizeAutos || 
          lfile_reader.m_pBaselinesCrossBuffer && lfile_reader.m_nLastBytesReadCross == lfile_reader.m_nBaselineSizeCross ){ 
          // these sizes lfile_reader.m_nBaselineSizeCross and lfile_reader.m_nBaselineSizeAutos are in bytes so they can be used to allocate required amount of memory, see :
          // float _Complex is the same as cufftComplex (see /home/msok/Desktop/SOFTWARE/logbook/20230330__Complex_vs_cufftComplex.odt ) so I can cast (tested) one to another (do not need to copy element by element !!!)          
          
          printf("INFO : processing timestamp = %d\n",timestamp);
                                        
          // (char*)lfile_reader.m_pBaselinesAutoBuffer, lfile_reader.m_nLastBytesReadAuto - pointer and size of buffer for Autos
          // (char*)lfile_reader.m_pBaselinesCrossBuffer, lfile_reader.m_nLastBytesReadCross - pointer and size of buffer for Cross-Corrs
          // m_nOneBaselineSizeCross = m_nChannels*sizeof(DATA_TYPE _Complex); // was *sizeof(DATA_TYPE _Complex); but this should be number of elements not bytes, so I am remvoing it to check if it still works ok (2022-0
          // m_nOneBaselineSizeAutos = m_nChannels*sizeof(float); // was *sizeof(float); but this should be number of elements not bytes, so I am remvoing it to check if it still works ok (2022-05-04)
          // m_nBaselineSizeAutos = m_nBaselinesAuto*m_nOneBaselineSizeAutos;
          // m_nBaselineSizeCross = n_nBaselines*m_nOneBaselineSizeCross;
          
 
          // get correlation matrix for each channel and process it :   
          
          if( pRealTimeImagerCPU ){
             for(int ch=0;ch<n_ch;ch++){
               // get correlation matrix for channel ch from current buffer 
               if( lfile_reader.GetCurrentCorrMatrix( ch, corr_matrix_real, corr_matrix_imag, 0 ) ){
                  bool bWriteCorrMatrix = false; // if interesting event found or trigger comes 

                  if( gDoImaging ){
                     printf("DEBUG : imaging %d channels\n",n_ch);
                     if( pRealTimeImagerCPU ){
                        if( gImageChannel<0 || gImageChannel==ch ){
                           pRealTimeImagerCPU->Image( ch, corr_matrix_real, corr_matrix_imag, lfile_reader.m_CurrentTimestampUnixTime );
                        }else{
                           printf("WARNING : imaging of channel %d is not required\n",ch);
                        }
                     }else{
                        printf("ERROR in code : pRealTimeImagerCPU pointer to imaging object is NULL -> cannot do imaging !!!\n");
                     }
                  }else{
                     printf("WARNING : imaging is not required\n");
                  }

                  // save correlation matrix if required :
                  if( gWriteCorrMatrixAll || bWriteCorrMatrix ){
                     char szFitsReal[1024],szFitsImag[1024];
  
                     sprintf(szFitsReal,"%s/%08d/%03d/corr_matrix_ch%03d_vis_real.fits",gOutDir.c_str(),timestamp,ch,ch);
                     sprintf(szFitsImag,"%s/%08d/%03d/corr_matrix_ch%03d_vis_imag.fits",gOutDir.c_str(),timestamp,ch,ch);

                     if( corr_matrix_real.WriteFits( szFitsReal ) ){
                        printf("ERROR : could not write output FITS file %s\n",szFitsReal);
                     }else{
                        printf("OK : FITS file %s written OK\n",szFitsReal);
                     }

                     if( corr_matrix_imag.WriteFits( szFitsImag ) ){
                        printf("ERROR : could not write output FITS file %s\n",szFitsImag);
                     }else{
                        printf("OK : FITS file %s written OK\n",szFitsImag);
                     }
                  }else{
                     if( CLfileReader::m_DebugLevel > 1 ){
                        printf("INFO : writting all correlation matrices is not required\n");
                     }
                  }
               }
            }
         }else{
            if( gGPUImaging ){
#ifdef IMAGER_HIP
               if( strlen(gCalibrationSolutionsFilename.c_str()) ){
                  // bool ApplyCalibration( CLfileReader& lfile_reader, CPacerImagerMultiFreqHip* pRealTimeImagerGPU, int pol )
                  ApplyCalibration( lfile_reader, pRealTimeImagerGPU, 0 );
               }

               // Copy external CPU buffer to internal GPU buffers inside class CPacerImagerMultiFreqHip :
               pRealTimeImagerGPU->CopyInputVisibilities( lfile_reader.m_pBaselinesAutoBuffer, (gpufftComplex*)lfile_reader.m_pBaselinesCrossBuffer );
                                                                         
               // perform imaging on n_ch freq. channels (blocks in cuFFT plan many) :                                                           
               pRealTimeImagerGPU->gridding_imaging_multi_freq( n_ch, delta_freq_mhz, // input visibilites and their description 
                                                                gUniformWeighting.c_str()
                                                              );
               if( (gSaveFitsFilesFreq>0 && (timestamp % gSaveFitsFilesFreq)==0 ) && CPacerImager::m_SaveFilesLevel >= SAVE_FILES_ALL ){
                  // saving FITS files is only when required and when special option is enabled 
                  if( pRealTimeImagerGPU->CopyImagesGpu2Cpu() > 0 )
                  {
                     for(int ch=0;ch<n_ch;ch++){
                         pRealTimeImagerGPU->SaveSkyImage(ch, timestamp, lfile_reader.m_CurrentTimestampUnixTime, gImageSize, gOutputDirectory.c_str(), gSaveImaginaryFITS );
                     }
                  }                                     
               }
#else
            printf("ERROR : pacer_imager_multi program was compiled in CPU-only version (add -DUSE_HIP=ON in cmake or gpu option for build.sh to enable GPU code)\n");
            exit(-1);
#endif            
            }
         }

          
          // warning : may need slightly different UVW for different channels, right ?
          // for the start we can use the same UVW=CONST
      }          
      timestamp++;
   }
  PACER_PROFILER_END("Full execution of imager (including I/O) took ")
  
  if( pRealTimeImagerCPU ){
     delete pRealTimeImagerCPU;
  }
#ifdef IMAGER_HIP
  if( pRealTimeImagerGPU ){
     delete pRealTimeImagerGPU;
  }
#endif
}

#ifdef IMAGER_HIP
bool ApplyCalibration( CLfileReader& lfile_reader, CPacerImagerMultiFreqHip* pRealTimeImagerGPU, int pol )
{
   int pol_idx = 0; // TODO : can just parameter int pol be used ???
//   if ( strcmp(szPol,"Y") == 0 ){
//      pol_idx = 3;
//   }

   CCalSols& calsol = pRealTimeImagerGPU->m_CalibrationSolutions;
   if( calsol.size() < lfile_reader.m_nAnts ){
      printf("ERROR : number of calibration solutions = %d is smaller than number of antennas -> cannot apply cal. solutions !!!\n",int(calsol.size()));
      return false;
   }
   
   
   for( int channel=0; channel<lfile_reader.m_nChannels; channel++){      
      // WARNING : order of loops changed from ant1, ant2 to ant2, ant1 to account for the transpose which is done in CLfileReader::GetCurrentCorrMatrix
      //           but not in the GPU version !!! May now be changed in gridding function !
      for(int ant2=0;ant2<lfile_reader.m_nAnts;ant2++){
         for(int ant1=(ant2);ant1<lfile_reader.m_nAnts;ant1++){
            if( ant1 != ant2 ){
               float _Complex* baseline_data = lfile_reader.GetCrossCorrData( ant1, pol, ant2, pol );

               float re = crealf( baseline_data[channel] );
               float im = cimagf( baseline_data[channel] );
               
               std::complex<double> z(re,im);

               std::complex<double> g1 = calsol[ant2].m_cal[pol_idx]; // or y ? TODO , TBD 
               std::complex<double> g2_c = std::conj( calsol[ant1].m_cal[pol_idx] ); // or x ? TODO , TBD
               std::complex<double> z_cal = (g1)*z*(g2_c);
               
               baseline_data[channel] = z_cal.real() + z_cal.imag()*_Complex_I; // https://www.gnu.org/software/libc/manual/html_node/Complex-Numbers.html 
            }else{
               float* auto_data = lfile_reader.GetAutoCorrData( ant1, pol );
               float re = auto_data[channel];
               float im = 0.00;
               
               std::complex<double> z(re,im);               
               std::complex<double> g1 = calsol[ant2].m_cal[pol_idx]; // or y ? TODO , TBD 
               std::complex<double> g2_c = std::conj( calsol[ant1].m_cal[pol_idx] ); // or x ? TODO , TBD
               std::complex<double> z_cal = (g1)*z*(g2_c);
               auto_data[channel] = z_cal.real();
            }
         }
      }
   }
   
   return true;   
}
#endif
