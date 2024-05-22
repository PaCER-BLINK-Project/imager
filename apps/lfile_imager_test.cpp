#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <string>

#include <myfile.h>
#include <bg_fits.h>
#include <lfile_reader.h>
#include "realtime_image_analyser.h"

using namespace std;


int n_ants   = 128;
int n_inputs = 256;
int n_ch  = 768;
int n_pols = 1; // 4 moved to n_inputs = 256 
int n_avg = 100; // number of averages executed by Randall's correlator 
long int n_timesteps = -1;
int baseline = 0;

// input basename :
string gBasename = "last_dump";
string gOutDir   = "fits/";
string gOutDirVis;

// configuration mapping file :
string gMappingFile = "instr_config.txt";
string gAntennaPositionsFile = "antenna_locations.txt";

// real-time imaging 
bool gDoImaging=false;
int  gImageChannel=-1; // <0 -> image all channels 
int  n_pix=180; // image size (n_pix x n_pix)
bool gIncludeAutos=false; // include auto-correlations in imaging 

string out_filename = "out.txt";
string out_amp_vs_time_file = "amp_vs_time.txt";
bool gDumpReIm = false;
double gUnixTime = 1617786293.317148;
double gFreqMHz = 204*(400.00/512.00);
double delta_freq_mhz = (400.00/512.00)*(32.00/27.00);
int gTimestep=0;

string szOutFits;
CBgFits* pOutFits = NULL;
int yFitsSize = 20970;
string gObject="TEST";
double gRA = 77.01525833;
double gDEC = 26.06106111;
bool gOneBigFITS = false;
bool gNormalise = false;

// what to write and what not :
bool gWriteCorrMatrixAll = false;
int  gTimestampsToProcess = -1;

#define DATA_TYPE float

void print_parameters()
{
   printf("#####################################\n");
   printf("PARAMETERS :\n");
   printf("#####################################\n");
   printf("basename     = %s\n",gBasename.c_str());
   printf("output dir   = %s\n",gOutDir.c_str());
   printf("output vis dir = %s\n",gOutDirVis.c_str());
   printf("Configuration files:\n");
   printf("\tMapping file = %s\n",gMappingFile.c_str());
   printf("\tAntenna positions file = %s\n",gAntennaPositionsFile.c_str());
   printf("Data :\n");
   printf("\tFrequency = %.4f [MHz]\n",gFreqMHz);
   printf("Imaging :\n");
   printf("\tImage size = %d x %d\n",n_pix,n_pix);
   printf("\tImage channel = %d\n",gImageChannel);
   printf("\tInclude autos = %d\n",gIncludeAutos);
   printf("Amount of data to process:\n");
   printf("\tN timestamps = %d\n",gTimestampsToProcess);
   printf("Dump RE/IM  = %d\n",gDumpReIm);
   printf("outfilename = %s\n",out_filename.c_str());
   printf("N_inputs = %d\n",n_inputs);
   printf("N_ants   = %d\n",n_ants);
   printf("N_pols   = %d\n",n_pols);
   printf("N_timesteps = %ld\n",n_timesteps);
   printf("N_channels  = %d\n",n_ch);
   printf("N_avg       = %d\n",n_avg);
   printf("Baseline = %d\n",baseline);
   printf("Amp vs. time file = %s\n",out_amp_vs_time_file.c_str());
   printf("Output FITS file  = %s\n",szOutFits.c_str());
   printf("\tOne big FITS = %d\n",gOneBigFITS);
   printf("\tOBJECT = %s\n",gObject.c_str());
   printf("\t(RA,DEC) = (%.4f,%.4f) [deg]\n",gRA,gDEC);
   printf("\tySize    = %d\n",yFitsSize);
   printf("\tunix time = %.6f\n",gUnixTime);
   printf("Save transposed correlation matrix = %d\n",CLfileReader::m_bTransposed);   
   printf("Verbosity level = %d\n",CLfileReader::m_DebugLevel);
   printf("Write all corr. matrices = %d\n",gWriteCorrMatrixAll);
   printf("#####################################\n");
}

void usage()
{
   printf("lfile_dumper LFILE -a ANTS -i -b BASELINE -c N_CHANNELS -o OUTFILENAME -m OUT_AMP_VS_TIME_FILE -p N_POLS -f OUT_FITS -y Y_FITS_SIZE -A N_AVG -F FREQ_MHZ\n");
   printf("-i : enable real-time imaging\n");
   printf("-C IMAGE_CHANNEL : to image only one channel [default %d], <0 -> all\n",gImageChannel);
   printf("-n N_TIMESTAMPs : maximum number of timestamps to process [default all]\n");
   printf("-r : dump RE/IM instead of default MAG/PHASE\n");
   printf("-y Y_FITS_SIZE : output Y size of FITS file [default %d]\n",yFitsSize);
   printf("-u UNIX_TIME\n");
   printf("-F FREQ_CHANNEL\n");
   printf("-N : normalise by mean power spectrum\n");
   printf("-I mapping_file [default %s]\n",gMappingFile.c_str());
   printf("-T saves transpose correlation matrix [default %d]\n",CLfileReader::m_bTransposed);
   printf("-t TIMESTEP [default %d]\n",gTimestep);
   printf("-O OUTDIR   [default %s]\n",gOutDir.c_str());
   printf("-V OUTDIR VISIBILITY  [default %s]\n",gOutDirVis.c_str());
   printf("-A antenna_positions.txt : path to file with antenna positions (4 columns)\n");
   printf("-s N_PIXELS : number of pixels for image size [default %d]\n",n_pix);
   printf("-w write correlation matrix\n");
   printf("-E : include auto-correlations [default %d]\n",gIncludeAutos);
   exit(0);
}

void parse_cmdline(int argc, char * argv[]) {
   char optstring[] = "ha:ip:t:b:c:o:m:vrf:By:A:u:F:TNI:o:V:O:s:n:C:wE";
   int opt,opt_param,i;

   while ((opt = getopt(argc, argv, optstring)) != -1) {
      switch (opt) {
         case 'a':
            n_ants = atol(optarg);
            n_inputs = n_ants*2;
            break;

         case 'A':
            if( optarg ){
               gAntennaPositionsFile = optarg;
            }
            break;

         case 'b':
            baseline = atol(optarg);
            break;

         case 'B':
            gOneBigFITS = true;
            break;

         case 'c':
            n_ch = atol(optarg);
            break;

         case 'C':
            if( optarg ){
               gImageChannel = atol(optarg);
            }
            break;

         case 'E':
            gIncludeAutos = true;
            break;

         case 'i':
            gDoImaging = true;
            break;

         case 'I':
            if( optarg ){
               gMappingFile = optarg;
            }
            break;

         case 'f':
            szOutFits = optarg;
            break;

         case 'F':
            gFreqMHz = atof( optarg );
            break;

         case 'n':
            if( optarg ){
               gTimestampsToProcess = atol(optarg);
            }
            break;

         case 'o':
            out_filename = optarg;
            break;

         case 'p':
            n_pols = atol(optarg);
            break;
            
         case 't':
            gTimestep = atol(optarg);
            break;

         case 'm':
            out_amp_vs_time_file = optarg;
            break;
            
         case 'O':
            if( optarg ){
               gOutDir = optarg;
            }
            break;
            
         case 'V':
            if( optarg ){
               gOutDirVis = optarg;
            }
            break;
            
         case 'r':
            gDumpReIm = true;
            break;

         case 's':
            if( optarg ){
               n_pix = atol( optarg );
            }
            break;

         case 'y':
            yFitsSize = atol(optarg);
            break;

         case 'u':
            gUnixTime = atof(optarg);
            break;

         case 'N':
            gNormalise = true;
            break;

         case 'T':
            CLfileReader::m_bTransposed = true;
            break;

         case 'w':
            gWriteCorrMatrixAll = true;
            break;            
             
         case 'h':
            usage();
            break;

         case 'v':
            CLfileReader::m_DebugLevel++;
            break;
            
         default:  
            fprintf(stderr,"Unknown option %c\n",opt);
            usage();
      }
   }
}


int main(int argc, char* argv[])
{
   if( argc >= 2 ){
      gBasename = argv[1];
   }
   
   parse_cmdline( argc, argv );
   print_parameters();
   fflush(stdout);
   
   // set parameters for imaging (even if not done):
   CRealTimeImageAnalyser* pRealTimeImager = NULL;
   if( gDoImaging ){
      pRealTimeImager = new CRealTimeImageAnalyser( n_ch, n_pix, gFreqMHz, delta_freq_mhz, gIncludeAutos );

      // set parameters for imaging (even if not done):
      pRealTimeImager->SetParameters( gAntennaPositionsFile.c_str(), true, gOutDir.c_str() );
   }


   CBgFits corr_matrix_real( n_ants, n_ants ), corr_matrix_imag( n_ants, n_ants );
   
   CLfileReader lfile_reader( gBasename.c_str(), n_ants, n_ch, n_timesteps, gMappingFile.c_str() );
   // int CLfileReader::GetCorrMatrix( int timestamp, int channel, CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, int pol )
   lfile_reader.Open();

   if( strlen(gOutDirVis.c_str()) <= 0 ){
      printf("INFO : writting of visibilities is not required (use option -V to specify output directory and enable it)\n");
   }
   
   // fetch next timestamp get all correlation matrices and save 
   int timestamp = 0;
   while( lfile_reader.ReadNextTimestamp() && (gTimestampsToProcess<=0 || timestamp<gTimestampsToProcess) ){
   
      // save visibilities (each timestamps into a separate L-file) :
      if( strlen(gOutDirVis.c_str()) > 0 ){
         printf("DEBUG : saving visibilities to directory %s is requested\n",gOutDirVis.c_str());
         
         if( lfile_reader.m_pBaselinesAutoBuffer && lfile_reader.m_nLastBytesReadAuto == lfile_reader.m_nBaselineSizeAutos ){
            char szVisAutoCorrFile[1024];
            sprintf(szVisAutoCorrFile,"%s/vis_autocorr_time%05d.LACSPC",gOutDirVis.c_str(),timestamp);
         
            printf("DEBUG : saving autocorrelations for timestamp %d to file %s\n",timestamp,szVisAutoCorrFile);
            
            int written = MyFile::WriteFile( szVisAutoCorrFile , (char*)lfile_reader.m_pBaselinesAutoBuffer, lfile_reader.m_nLastBytesReadAuto );
            if( written == lfile_reader.m_nLastBytesReadAuto ){
               printf("OK : written %d as expected\n",int(written));
            }else{
               printf("ERROR : written %d vs. expected %d\n",written,lfile_reader.m_nLastBytesReadAuto);
            }
         }

         if( lfile_reader.m_pBaselinesCrossBuffer && lfile_reader.m_nLastBytesReadCross == lfile_reader.m_nBaselineSizeCross ){
            char szVisCrossCorrFile[1024];
            sprintf(szVisCrossCorrFile,"%s/vis_crosscorr_time%05d.LACSPC",gOutDirVis.c_str(),timestamp);
         
            printf("DEBUG : saving autocorrelations for timestamp %d to file %s\n",timestamp,szVisCrossCorrFile);
            
            int written = MyFile::WriteFile( szVisCrossCorrFile , (char*)lfile_reader.m_pBaselinesCrossBuffer, lfile_reader.m_nLastBytesReadCross );
            if( written == lfile_reader.m_nLastBytesReadCross ){
               printf("OK : written %d as expected\n",int(written));
            }else{
               printf("ERROR : written %d vs. expected %d\n",written,lfile_reader.m_nLastBytesReadCross);
            }
         }
      }

      // get correlation matrix for each channel and process it :   
      printf("DEBUG : imaging %d channels\n",n_ch);
      for(int ch=0;ch<n_ch;ch++){
         // get correlation matrix for channel ch from current buffer 
         if( lfile_reader.GetCurrentCorrMatrix( ch, corr_matrix_real, corr_matrix_imag, 0 ) ){
            bool bWriteCorrMatrix = false; // if interesting event found or trigger comes 
            
            if( gDoImaging ){
               printf("DEBUG : imaging channel %d\n",ch);
               if( pRealTimeImager ){
                  // just test:
                  //if( ch==0 ){
                  //   corr_matrix_real.WriteFits("test_re.fits");
                  //}
                  
                  if( gImageChannel<0 || gImageChannel==ch ){
                     pRealTimeImager->Image( ch, corr_matrix_real, corr_matrix_imag, lfile_reader.m_CurrentTimestampUnixTime );
                  }else{
                     printf("WARNING : imaging of channel %d is not required\n",ch);
                  }
               }else{
                  printf("ERROR in code : pRealTimeImager pointer to imaging object is NULL -> cannot do imaging !!!\n");
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
               if( CLfileReader::m_DebugLevel > 0 ){
                  printf("WARNING : writting all correlation matrices is not required\n");
               }
            }            
         }
      }  
      
      timestamp++;          
   }
      
   printf("Processed %d timestamps out of all %d\n",timestamp,int(lfile_reader.m_nTimestamps));
   
   if( pRealTimeImager ){
      delete pRealTimeImager;
   }
}
