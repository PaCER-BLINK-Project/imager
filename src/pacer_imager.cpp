#include "pacer_imager.h"
#include "pacer_common.h"
#include <bg_fits.h>

// FFTW, math etc :
#include <fftw3.h>
#include <math.h>

// local defines :
#include "pacer_imager_defs.h"

// msfitslib library :
#include <myfile.h>
#include <libnova_interface.h>

#ifdef _PACER_PROFILER_ON_
#include <mydate.h>
#endif

// AstroIO for Visibilities class :
#include <astroio.hpp>

// TEST OPTIONS to compare with MIRIAD image
// see memos : PAWSEY/PaCER/logbook/20220305_pacer_imager_validation.odt, MIRIAD natural weighting (sup=0) etc:
// invert vis=chan_204_20211116T203000_yx.uv map=chan_204_20211116T203000_iyx.map imsize=180,180 beam=chan_204_20211116T203000_iyx.beam  sup=0 options=imaginary stokes=yx select='uvrange(0.0,100000)'
bool CPacerImager::m_bCompareToMiriad = false;

// debug level : see pacer_imager_defs.h for SAVE_FILES_NONE etc
int CPacerImager::m_ImagerDebugLevel = IMAGER_INFO_LEVEL; // IMAGER_ALL_MSG_LEVEL;

// level of saving intermediate and test files , see pacer_imager_defs.h for defines SAVE_FILES_NONE
int CPacerImager::m_SaveFilesLevel = SAVE_FILES_ALL;

// can also save control files every N-th file 
int CPacerImager::m_SaveControlImageEveryNth=-1;

// show final image statistics :
bool CPacerImager::m_bPrintImageStatistics = false; // default disabled to make imaging as fast as possible

void CPacerImager::SetDebugLevel( int debug_level )
{
   CPacerImager::m_ImagerDebugLevel = debug_level;
   gBGPrintfLevel = debug_level;
}

void CPacerImager::SetFileLevel( int filesave_level )
{
   CPacerImager::m_SaveFilesLevel = filesave_level;
}

int CPacerImager::UpdateFlags()
{
   int count_flagged = 0;
   
   if( m_FlaggedAntennas.size() > 0 ){
      if( m_MetaData.m_AntennaPositions.size() > 0 ){
         // flagging antennas in the list :
         for(int i=0;i<m_FlaggedAntennas.size();i++){
            int ant_index = m_FlaggedAntennas[i];
            
            if( ant_index >= 0 && ant_index < m_MetaData.m_AntennaPositions.size() ){
               m_MetaData.m_AntennaPositions[ant_index].flag = 1;
               count_flagged++;
            }
         }
         PRINTF_DEBUG("CPacerImager::SetFlaggedAntennas : Flagged %d antennas in the imager object\n",count_flagged);
      }else{
         PRINTF_DEBUG("CPacerImager::SetFlaggedAntennas : No antennas in object m_MetaData.m_AntennaPositions\n");
      }
   }

   return count_flagged;
}

void CPacerImager::SetFlaggedAntennas( vector<int>& flagged_antennas )
{
   m_FlaggedAntennas = flagged_antennas;
   UpdateFlags();
}

CPacerImager::CPacerImager()
: m_bInitialised(false), m_Baselines(0), m_pSkyImageReal(NULL), m_pSkyImageImag(NULL), m_pSkyImageRealTmp(NULL), m_pSkyImageImagTmp(NULL), m_bLocalAllocation(false), m_SkyImageCounter(0),
  m_in_buffer(NULL), m_out_buffer(NULL), m_in_size(-1), m_out_size(-1), m_bIncludeAutos(false),
  m_uv_grid_counter(NULL), m_uv_grid_real(NULL), m_uv_grid_imag(NULL), m_nAntennas(0),
  u_mean(0.00), u_rms(0.00), u_min(0.00), u_max(0.00), 
  v_mean(0.00), v_rms(0.00), v_min(0.00), v_max(0.00), 
  w_mean(0.00), w_rms(0.00), w_min(0.00), w_max(0.00)  
{
   m_PixscaleAtZenith = 0.70312500; // deg for ch=204 (159.375 MHz) EDA2 
}

CPacerImager::~CPacerImager()
{
   CleanLocalAllocations();
   
   if( m_in_buffer ){
      fftw_free(  (fftw_complex*)m_in_buffer );
   }
   if( m_out_buffer ){
      fftw_free(  (fftw_complex*)m_out_buffer );
   }
}


void CPacerImager::CleanLocalAllocations()
{
   if( m_bLocalAllocation )
   { // only remove when locally allocated, as it can also be passed from outside using SetOutputImagesExternal()
      if( m_pSkyImageReal )
      {
         delete m_pSkyImageReal;
         m_pSkyImageReal = NULL;
      }
      if( m_pSkyImageImag )
      {
         delete m_pSkyImageImag;
         m_pSkyImageImag = NULL;
      }
      
      m_bLocalAllocation = false;
   }
   
   // gridded visibilites are always local allocations (cannot be passed from external function)
   if( m_uv_grid_counter ){
      delete m_uv_grid_counter;
      m_uv_grid_counter = NULL;
   }
   if( m_uv_grid_real ){
      delete m_uv_grid_real;
      m_uv_grid_real = NULL;
   }
   if( m_uv_grid_imag ){
      delete m_uv_grid_imag;
      m_uv_grid_imag = NULL;
   }   
   
   // these are always only allocated locally :
   if( m_pSkyImageRealTmp ){
      delete m_pSkyImageRealTmp;
      m_pSkyImageRealTmp = NULL;
   }
   if( m_pSkyImageImagTmp ){
      delete m_pSkyImageImagTmp;
      m_pSkyImageImagTmp = NULL;
   }
}

bool CPacerImager::AllocOutPutImages( int sizeX, int sizeY )
{
   bool bRet = false;
   if( !m_pSkyImageReal )
   {
      m_pSkyImageReal = new CBgFits( sizeX, sizeY );   
      m_bLocalAllocation = true;
      bRet = true;
   }else
   {
      CheckSize( *m_pSkyImageReal, sizeX, sizeY );
   }

   if( !m_pSkyImageImag )
   {
      m_pSkyImageImag = new CBgFits( sizeX, sizeY );
      m_bLocalAllocation = true;
      bRet = true;
   }else
   {
      CheckSize( *m_pSkyImageImag, sizeX, sizeY );
   }

   // temporary buffers, always allocated locally, but flag m_bLocalAllocation is only for output images which can go out of to external calling routines:
   if( !m_pSkyImageRealTmp )
   {
      m_pSkyImageRealTmp = new CBgFits( sizeX, sizeY );   
      bRet = true;
   }else
   {
      CheckSize( *m_pSkyImageRealTmp, sizeX, sizeY );
   }

   if( !m_pSkyImageImagTmp )
   {
      m_pSkyImageImagTmp = new CBgFits( sizeX, sizeY );   
      bRet = true;
   }else
   {
      CheckSize( *m_pSkyImageImagTmp, sizeX, sizeY );
   }

   // Allocate input and output buffers for FFTW : 
   int size = sizeX*sizeY;
   if( !m_in_buffer || m_in_size != size ){
       // WARNING : cannot be in constructor where size is not yet known :       
       if( m_in_buffer && m_in_size != size && m_in_size > 0 ){
          // if image size changes we can use this code
          PRINTF_INFO("INFO : freeing m_in_buffer memory buffer");
          fftw_free( (fftw_complex*)m_in_buffer );
       }
  
       PRINTF_INFO("INFO : allocating m_in_buffer buffer of size %d * %d bytes\n",size,int(sizeof(fftw_complex)));
       m_in_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
       m_in_size = size;
   }
   if( !m_out_buffer || m_out_size != size ){
       // WARNING : cannot be in constructor where size is not yet known :
       if( m_out_buffer && m_out_size != size && m_out_size > 0 ){
          // if image size changes we can use this code
          PRINTF_INFO("INFO : freeing m_out_buffer memory buffer");
          fftw_free( m_out_buffer );
       }

       PRINTF_INFO("INFO : allocating m_out_buffer buffer of size %d * %d bytes\n",size,int(sizeof(fftw_complex)));
       m_out_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
       m_out_size = size;
   }

   
   return bRet;
}


bool CPacerImager::AllocGriddedVis( int sizeX, int sizeY )
{
   if( !m_uv_grid_counter ){
      m_uv_grid_counter = new CBgFits( sizeX, sizeY );
   }else{
      CheckSize( *m_uv_grid_counter, sizeX, sizeY );
   }

   if( !m_uv_grid_real ){
      m_uv_grid_real = new CBgFits( sizeX, sizeY );
   }else{
      CheckSize( *m_uv_grid_real, sizeX, sizeY );
   }
   
   if( !m_uv_grid_imag ){
      m_uv_grid_imag = new CBgFits( sizeX, sizeY );
   }else{
      CheckSize( *m_uv_grid_imag, sizeX, sizeY );
   }

   return true;
}

void CPacerImager::SetOutputImagesExternal( CBgFits* pSkyImageRealExt, CBgFits* pSkyImageImagExt )
{
   CleanLocalAllocations();
      
   m_pSkyImageReal = pSkyImageRealExt;
   m_pSkyImageImag = pSkyImageImagExt;
   m_bLocalAllocation = false;
}

int CPacerImager::ReadAntennaPositions( bool bConvertToXYZ )
{
   m_nAntennas = m_MetaData.m_AntennaPositions.ReadAntennaPositions( m_ImagerParameters.m_AntennaPositionsFile.c_str(), bConvertToXYZ  );
   UpdateFlags(); // if antenna positions are read now for the first time, the flags need to be updated
   PRINTF_INFO("INFO : read %d antenna positions from file %s\n",m_nAntennas,m_ImagerParameters.m_AntennaPositionsFile.c_str());

   return m_nAntennas;
}



void CPacerImager::Initialise()
{
  if( !m_bInitialised )
  {
     m_bInitialised = true;
     
    // test file with antenna positions can be used to overwrite whatever was in .metafits
    if( strlen( m_ImagerParameters.m_AntennaPositionsFile.c_str() ) && MyFile::DoesFileExist(  m_ImagerParameters.m_AntennaPositionsFile.c_str() ) )
    {
       bool bConvertToXYZ = false;
       if( !m_ImagerParameters.m_bConstantUVW )
       { // if non-constant UVW -> non zenith phase centered all-sky image
          bConvertToXYZ = true;
       }
       if( m_ImagerParameters.m_bAntennaPositionsXYZ )
       { // text file already has XYZ in WG54 system - for example CASA dump 
          printf("INFO : antenna positions already in XYZ coordinate system (WG54) no need to convert\n");
          bConvertToXYZ = false;
       }
       
       // read antenna positions and do whatever else is necessary (update flags etc)
       ReadAntennaPositions( bConvertToXYZ );
       
       if( /*true ||*/ strlen( m_ImagerParameters.m_MetaDataFile.c_str() ) == 0 )
       { // only calculate UVW here when Metadata is not required
          // initial recalculation of UVW at zenith (no metadata provided -> zenith):       
          // WARNING : bool CPacerImager::CalculateUVW() - could be used, but it also calls this function itself which may cause infinite recursise call
          // m_Baselines = m_MetaData.m_AntennaPositions.CalculateUVW( m_U, m_V, m_W, (CPacerImager::m_SaveFilesLevel>=SAVE_FILES_DEBUG), m_ImagerParameters.m_szOutputDirectory.c_str(), m_bIncludeAutos );
          // UpdateParameters();
          CalculateUVW( true, false ); // bForce=true to force call of m_MetaData.m_AntennaPositions.CalculateUVW and , bInitialise=false to avoid call to this (CPacerImager::Initialise) function in a recursive way !
          PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines (include Autos = %d)\n",m_Baselines,m_bIncludeAutos);
       }else
       {
          printf("INFO : non-zenith pointing meta data is required to calculate UVW\n");
       }
    }else
    {
       PRINTF_WARNING("WARNING : antenna position file %s not specified or does not exist\n",m_ImagerParameters.m_AntennaPositionsFile.c_str());
    }
    
    // read all information from metadata 
    if( strlen( m_ImagerParameters.m_MetaDataFile.c_str() ) && MyFile::DoesFileExist( m_ImagerParameters.m_MetaDataFile.c_str() ) ){
      PRINTF_INFO("INFO : reading meta data from file %s\n",m_ImagerParameters.m_MetaDataFile.c_str());
      if( !m_MetaData.ReadMetaData( m_ImagerParameters.m_MetaDataFile.c_str() ) ){
         PRINTF_ERROR("ERROR : could not read meta data from file %s\n",m_ImagerParameters.m_MetaDataFile.c_str() );
      }
    }       
    
    
    if( strlen(m_CalibrationSolutions.m_filename.c_str()) > 0 ){
       if( MyFile::DoesFileExist(m_CalibrationSolutions.m_filename.c_str()) ){
          // read calibration solutions (if specified) :
          PRINTF_INFO("INFO : reading calibration solutions from file |%s|\n",m_CalibrationSolutions.m_filename.c_str());
          double sign_value = -1.00; // for EDA2 
          bool bInvertAmp = false;
          if( m_MetaData.eTelescopeName == eMWA ){
             // sign_value = 1.00; // for MWA use cal. solutions with the same sign of phase as in the file 
             bInvertAmp = true; // invert amplitudes when using MWA calibration solutions from .bin file
          }
          if( m_CalibrationSolutions.read_calsolutions( sign_value, bInvertAmp ) > 0 ) // sign_value=-1, and inverse amplitude is the best ?
          {
             m_CalibrationSolutions.show();
          }
       }else{
          PRINTF_WARNING("WARNING : calibration solutions file %s does not exist, ignoring as it may just be a template for multiple channel processing\n",m_CalibrationSolutions.m_filename.c_str());
       }
    }
  }
}

void CPacerImager::UpdateCalSolFile( const char* calsol_filename )
{
   m_CalibrationSolutions.m_filename = calsol_filename;
    
   if( strlen(m_CalibrationSolutions.m_filename.c_str()) > 0 ){
      if( MyFile::DoesFileExist(m_CalibrationSolutions.m_filename.c_str()) ){
         // read calibration solutions (if specified) :
         PRINTF_INFO("INFO : reading calibration solutions from file |%s|\n",m_CalibrationSolutions.m_filename.c_str());
         double sign_value = -1.00; // for EDA2 
         bool bInvertAmp = false;
         if( m_MetaData.eTelescopeName == eMWA ){
            // sign_value = 1.00; // for MWA use cal. solutions with the same sign of phase as in the file 
            bInvertAmp = true; // invert amplitudes when using MWA calibration solutions from .bin file
         }
         if( m_CalibrationSolutions.read_calsolutions( sign_value, bInvertAmp, true ) > 0 ) // sign_value=-1, and inverse amplitude is the best. true to force reading in this function
         {
            m_CalibrationSolutions.show();
         }
      }else{
         PRINTF_WARNING("WARNING : calibration solutions file %s does not exist, ignoring as it may just be a template for multiple channel processing\n",m_CalibrationSolutions.m_filename.c_str());
      }
   }
}

bool CPacerImager::CheckSize( CBgFits& image, int sizeX, int sizeY )
{
   if( image.GetXSize() != sizeX || image.GetYSize() != sizeY ){
      image.Realloc( sizeX, sizeY );      
      PRINTF_INFO("DEBUG : change of image size to (%d,%d) was required\n",sizeX,sizeY);
      
      return true;
   }
   
   // if image size was ok and nothing was required
   return false;
}

// See https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/ 
// for explanations why it is needed
void CPacerImager::fft_shift( CBgFits& dirty_image, CBgFits& out_image )
{
   int xSize = dirty_image.GetXSize();
   int ySize = dirty_image.GetYSize();

   // TODO : create member object m_tmp_image to avoid allocation every time this function is called 
   CBgFits tmp_image( xSize, ySize );
   
   int center_freq_x = int( xSize/2 );
   int center_freq_y = int( ySize/2 );
   
   int is_odd = 0;
   if ( (xSize%2) == 1 && (ySize%2) == 1 ){
      is_odd = 1;
   }

   // TODO : optimise similar to gridder.c in RTS or check imagefromuv.c , LM_CopyFromFFT which is totally different and may have to do with image orientation, but also is faster !!!
   // X (horizontal FFT shift) :
   for(int y=0;y<ySize;y++){ 
      float* tmp_data = tmp_image.get_line(y);
      float* image_data = dirty_image.get_line(y);
      
      // TODO / WARNING : lools like here for x=center_freq_x and images size = 2N -> center_freq_x = N -> center_freq_x+x can by N+N=2N which is outside image !!!
      for(int x=0;x<center_freq_x;x++){ // 2024-08-01 changed <= to < to fix the bug (small buffer overflow)
         tmp_data[center_freq_x+x] = image_data[x];
      }
      for(int x=(center_freq_x+is_odd);x<xSize;x++){
         tmp_data[x-(center_freq_x+is_odd)] = image_data[x];
      }      
   }

   for(int x=0;x<xSize;x++){ 
      for(int y=0;y<=center_freq_y;y++){ // check <= -> <
         out_image.setXY(x,center_freq_y+y,tmp_image.getXY(x,y));
      }
      for(int y=(center_freq_y+is_odd);y<ySize;y++){
         out_image.setXY( x , y-(center_freq_y+is_odd),tmp_image.getXY(x,y));
      }      
   }
}

// UV data are with DC in the center -> have to be FFTshfted to form input to FFT function :
// See https://www.gaussianwaves.com/2015/11/interpreting-fft-results-complex-dft-frequency-bins-and-fftshift/ 
// for explanations why it is needed
void CPacerImager::fft_unshift( CBgFits& dirty_image, CBgFits& out_image )
{
   int xSize = dirty_image.GetXSize();
   int ySize = dirty_image.GetYSize();
   
   // TODO : create member object m_tmp_image to avoid allocation every time this function is called
   CBgFits tmp_image( xSize, ySize );
   
   int center_freq_x = int( xSize/2 );
   int center_freq_y = int( ySize/2 );
   
   int is_odd = 0;
   if ( (xSize%2) == 1 && (ySize%2) == 1 ){
      is_odd = 1;
   }
   
   // TODO : optimise similar to gridder.c in RTS or check imagefromuv.c , LM_CopyFromFFT which is totally different and may have to do with image orientation, but also is faster !!!
   // X (horizontal FFT shift) :
   for(int y=0;y<ySize;y++){ 
      float* tmp_data = tmp_image.get_line(y);
      float* image_data = dirty_image.get_line(y);
      
      for(int x=0;x<center_freq_x;x++){ // check <= -> <
         tmp_data[center_freq_x+x+is_odd] = image_data[x];
      }
      for(int x=center_freq_x;x<xSize;x++){
         tmp_data[x-center_freq_x] = image_data[x];
      }      
   }

   for(int x=0;x<xSize;x++){ 
      for(int y=0;y<center_freq_y;y++){ // check <= -> <
         out_image.setXY( x, center_freq_y+y+is_odd, tmp_image.getXY(x,y));
      }
      for(int y=center_freq_y;y<ySize;y++){
         out_image.setXY(x,y-center_freq_y,tmp_image.getXY(x,y));
      }      
   }
}

bool CPacerImager::SaveSkyImage( const char* outFitsName, CBgFits* pFits, double unixtime /*=0.00*/ )
{
   if( !pFits ){
      PRINTF_ERROR("ERROR in code SaveSkyImage, pFits pointer not set\n");
      return false;
   }

   PRINTF_INFO("INFO : saving image %s\n",outFitsName);   
   pFits->SetFileName( outFitsName );
   
   // fill FITS header :
   // TODO :
   pFits->SetKeyword("TELESCOP","EDA2");
   
   // scripts/add_fits_header.py
   // TODO - use properly calculated values :
   // double pixscale = m_PixscaleAtZenith/3; // ??? why divided by 3 seems to be best ???
//   double pixscale = m_ImagerParameters.m_ImageFOV_degrees/pFits->GetXSize();
   double pixscale = m_ImagerParameters.m_PixsizeInRadians*(180.00/M_PI);
   if( m_bCompareToMiriad && false ){ // 2023-12-17 temporary disabled as WCS is not correctly saved then !!! see 20231215_repeat_processing_of_small_part_of_20230709.odt
       PRINTF_WARNING("WARNING (CPacerImager::SaveSkyImage) : MIRIAD-like option -> changing pixscale  %.8f -> %.8f\n",pixscale,m_PixscaleAtZenith);
       pixscale = m_PixscaleAtZenith;       
   }

   PRINTF_DEBUG("DEBUG : m_PixscaleAtZenith = %.6f [deg] -> pixscale = %.6f [deg] , FoV = %.4f [deg], ImageSize = %ld x %ld\n",m_PixscaleAtZenith,pixscale,m_ImagerParameters.m_ImageFOV_degrees,pFits->GetXSize(),pFits->GetYSize());
   // azh2radec 1581220006 mwa 0 90
   // (RA,DEC) = ( 312.07545047 , -26.70331900 )
   // 20.80503003133333333333
   // double ra_deg = 312.07545047; // = 20.80503003133333333333 hours ; // was 2.13673600000E+02;
   // double dec_deg = -2.67033000000E+01;
   
   // libnova_interface.h
   // void azh2radec( double az, double alt, time_t unix_time, double geo_long_deg, double geo_lat_deg, double& out_ra, double& out_dec );
   // for all-sky images at zenith :
   double ra_deg=-1000, dec_deg=-1000;
   if ( m_ImagerParameters.m_bConstantUVW ){
      PRINTF_DEBUG("DEBUG : calculating RA,DEC from (AZ,EL) = (0,90) [deg] for unixtime = %d, geo location (%.4f,%.4f)\n",int(unixtime),m_MetaData.geo_long,m_MetaData.geo_lat);
      ::azh2radec( 0.00, 90.00, unixtime, m_MetaData.geo_long, m_MetaData.geo_lat, ra_deg, dec_deg );
   }else{
      // use RA DEC from metafits or something like this 
      ra_deg = m_MetaData.raHrs*15.00;
      dec_deg = m_MetaData.decDegs;
      PRINTF_DEBUG("DEBUG : using RA,DEC = (%.8f,%.8f) [deg] from metafits\n",ra_deg,dec_deg);
   }
   
   int crpix1 = int(pFits->GetXSize()/2) + 1;
   pFits->SetKeyword("CTYPE1","RA---SIN");
   pFits->SetKeyword("CRPIX1", crpix1 );   
   pFits->SetKeywordFloat("CDELT1", -pixscale ); // WARNING / TODO : this may be related to image flip and can differ for different data -> really need to sort this out ASAP !!!
   pFits->SetKeywordFloat("CRVAL1", ra_deg ); // RA of the centre 
   pFits->SetKeyword("CUNIT1", "deg    " );

   int crpix2 = int(pFits->GetYSize()/2) + 1;
   pFits->SetKeyword("CTYPE2","DEC--SIN");
   pFits->SetKeyword("CRPIX2", crpix2 );   
   pFits->SetKeywordFloat("CDELT2", pixscale );
   pFits->SetKeywordFloat("CRVAL2", dec_deg ); // RA of the centre 
   pFits->SetKeyword("CUNIT2", "deg    " );
   
   // add UTC - also required for scripts pix2sky.py to work ok !
   struct tm gmtime_tm;
   double fraction_sec=0.00;
   time_t unix_time_time_t = (time_t)unixtime; 
   fraction_sec = unixtime - double(unix_time_time_t);
   if(gmtime_r( &unix_time_time_t, &gmtime_tm )){
      char tempstring[64];
      // 2023-06-01T10:42:50.1
      sprintf(tempstring,"%04u%02u%02u_%02u%02u%02u.%03d",gmtime_tm.tm_year+1900,(gmtime_tm.tm_mon+1),gmtime_tm.tm_mday,
                                                     gmtime_tm.tm_hour,gmtime_tm.tm_min,gmtime_tm.tm_sec,int(fraction_sec*1000.00));
       
      pFits->SetKeyword("DATE-OBS",tempstring);                                                      
   }else{
      printf("ERROR : could not convert unixtime = %.4f to UTC using gmtime_r\n",unixtime);
   }



   
   pFits->WriteFits( outFitsName );
   
   return true;
}

// Based on example : https://github.com/AccelerateHS/accelerate-examples/blob/master/examples/fft/src-fftw/FFTW.c
void CPacerImager::dirty_image( CBgFits& uv_grid_real_param, CBgFits& uv_grid_imag_param, CBgFits& uv_grid_counter,
                                bool bSaveIntermediate /*=false*/ , const char* szBaseOutFitsName /*=NULL*/, 
                                bool bSaveImaginary /*=true*/ , bool bFFTUnShift /*=true*/ )
{
   PACER_PROFILER_START

//   CBgFits uv_grid_real( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() ), uv_grid_imag( uv_grid_imag_param.GetXSize(), uv_grid_imag_param.GetYSize() );
  
// 2022-12-29 : temporary leaving the code commented out. If using old gridding (::gridding( ... )) functions is required there needs to be fft_unshift call
//              this cannont be done "inplace" so there will be a need for a temporary buffer to make it. 
//              member variables m_uv_grid_real_tmp , m_uv_grid_imag_tmp can be added of local temporary variable created and used
//   if( bFFTUnShift ){
//      fft_unshift( uv_grid_real_param, uv_grid_real );
//      fft_unshift( uv_grid_imag_param, uv_grid_imag );
//   }else{
      // TODO : this is totally suboptimal because it copies !!!
      // [0,0] = %.8f and [1,0] = %.8f\n", uv_grid_real_param.getXY(0,0),  uv_grid_real_param.getXY(1,0) );
//      uv_grid_real = uv_grid_real_param;
//      uv_grid_imag = uv_grid_imag_param;
//      printf("DEBUG : assigning UV grid using parameter FITS file , [0,0] = %.8f and [1,0] = %.8f\n", uv_grid_real.getXY(0,0),  uv_grid_real.getXY(1,0) );
//   }
// 
// 2022-12-29 :
   // using a temporary reference to avoid copy , see above comments in case old gridding(...) function + fft_unshift has to be used
   CBgFits& uv_grid_real = uv_grid_real_param;
   CBgFits& uv_grid_imag = uv_grid_imag_param;

   if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){   
      printf("DEBUG : saving UV grid after unshift is required\n");
      char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
      sprintf(szOutPutFitsRE,"%s/uv_grid_real_unshift.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
      sprintf(szOutPutFitsIM,"%s/uv_grid_imag_unshift.fits",m_ImagerParameters.m_szOutputDirectory.c_str());   
      uv_grid_real.WriteFits( szOutPutFitsRE );
      uv_grid_imag.WriteFits( szOutPutFitsIM );
   }


   int width = uv_grid_real.GetXSize();
   int height = uv_grid_real.GetYSize();
   int size = width*height;
   
   // Allocate input and output buffers
/*   if( !m_in_buffer || m_in_size != size ){
       // WARNING : cannot be in constructor where size is not yet known :       
       if( m_in_buffer && m_in_size != size && m_in_size > 0 ){
          // if image size changes we can use this code
          PRINTF_INFO("INFO : freeing m_in_buffer memory buffer");
          fftw_free( (fftw_complex*)m_in_buffer );
       }
       
       PRINTF_INFO("INFO : allocating m_in_buffer buffer of size %d * %d bytes\n",size,int(sizeof(fftw_complex)));
       m_in_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
       m_in_size = size;
   }
   if( !m_out_buffer || m_out_size != size ){
       // WARNING : cannot be in constructor where size is not yet known :
       if( m_out_buffer && m_out_size != size && m_out_size > 0 ){
          // if image size changes we can use this code
          PRINTF_INFO("INFO : freeing m_out_buffer memory buffer");
          fftw_free( m_out_buffer );
       }

       PRINTF_INFO("INFO : allocating m_out_buffer buffer of size %d * %d bytes\n",size,int(sizeof(fftw_complex)));
       m_out_buffer = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
       m_out_size = size;
   }*/
   // 2022-12-29 later : moved here to allocate m_pSkyImageRealTmp and m_pSkyImageImagTmp 
   //              it now allocates both temporary (Tmp) buffers and output images m_pSkyImageReal and m_pSkyImageImag
   //              and also to allocate m_in_buffer and m_out_buffer which are input and output buffers respectively for FFTW 
   AllocOutPutImages( uv_grid_real.GetXSize(), uv_grid_real.GetYSize() );
   
   // 
   float* real_data = uv_grid_real.get_data();
   float* imag_data = uv_grid_imag.get_data();

   // Copy in image data as real values for the transform.
   for(int i = 0; i < size; i++) {
      double re = real_data[i];
      double im = imag_data[i];
   
      ((fftw_complex*)m_in_buffer)[i][0]	= re;
      ((fftw_complex*)m_in_buffer)[i][1]	= im;
   }

   // should be the same and it is the same !
/*   for(int y=0;y<height;y++){
      for(int x=0;x<width;x++){
         double re = uv_grid_real.getXY(x,y);
         double im = uv_grid_imag.getXY(x,y);
         
         int pos = y*width + x;
         in_buffer[pos][0]   = re;
         in_buffer[pos][1]   = im;
      }
   }*/
   
   // Transform to frequency space.
   // https://www.fftw.org/fftw3_doc/Complex-Multi_002dDimensional-DFTs.html
   // WARNING : did not work OK:  ant2-ant1 in CalculateUVW in antenna_positions.cpp and this works OK with FFT_BACKWARD :
   //           Correct orientation is with V -> -V and FFTW_FORWARD - not clear why it is like this , see ::gridding (  double v = -fits_vis_v.getXY(ant1,ant2) / wavelength_m; )   
   // ???? Is there any good reaons for this - see also gridder.c and 	imagefromuv.c , LM_CopyFromFFT in RTS, especially the latter does some totally crazy re-shuffling from FFT output to image ...            
   fftw_plan pFwd = fftw_plan_dft_2d( width, height, (fftw_complex*)m_in_buffer, (fftw_complex*)m_out_buffer, FFTW_FORWARD, FFTW_ESTIMATE); // was FFTW_FORWARD or FFTW_BACKWARD ???
//   printf("WARNING : fftw BACKWARD\n");
  
   // neither agrees with MIRIAD :
   // fftw_plan pFwd = fftw_plan_dft_2d( width, height, in_buffer, out_buffer, FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw_execute(pFwd);
   fftw_destroy_plan(pFwd);
   
   
   
//   CBgFits out_image_real( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() ), out_image_imag( uv_grid_imag_param.GetXSize(), uv_grid_imag_param.GetYSize() );
   // 2022-12-29 : moved here to allocate m_pSkyImageRealTmp and m_pSkyImageImagTmp :
   //              it now allocates both temporary (Tmp) buffers and output images m_pSkyImageReal and m_pSkyImageImag
   // AllocOutPutImages( uv_grid_real_param.GetXSize(), uv_grid_real_param.GetYSize() );
   m_pSkyImageRealTmp->SetValue( 0.00 );
   m_pSkyImageImagTmp->SetValue( 0.00 );

   // copy resulting data from out_buffer to out_image_real :
   // WARNING : this is image in l,m = cos(alpha), cos(beta) coordinates and still needs to go to SKY COORDINATES !!!
   float* out_data_real = m_pSkyImageRealTmp->get_data();
   float* out_data_imag = m_pSkyImageImagTmp->get_data();
//   double fnorm = 1.00/sqrt(size); //normalisation see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220119_testing_new_versions_dirty_image_polishing.odt
   double counter_sum = uv_grid_counter.Sum();
   double fnorm = 1.00/counter_sum; // see RTS : /home/msok/mwa_software/RTS_128t/src/newgridder.cu SumVisibilityWeights and gridKernel.c:650 
                                              // also read TMS (Thomson, Moran, Swenson) about this 
   PRINTF_DEBUG("DEBUG : size = %d (%d x %d), fnorm = %e (counter sum = %.8f)\n",size,width,height,fnorm,counter_sum);
   for(int i = 0; i < size; i++) {
//      out_data[i] = out_buffer[i][0]*out_buffer[i][0] + out_buffer[i][1]*out_buffer[i][1]; // amplitude ?
     out_data_real[i] = ((fftw_complex*)m_out_buffer)[i][0]*fnorm; // real 
     out_data_imag[i] = ((fftw_complex*)m_out_buffer)[i][1]*fnorm; // imag
   }   

   char outDirtyImageReal[1024],outDirtyImageImag[1024];   
   
   if( bSaveIntermediate ){ // I will keep this if - assuming it's always TRUE, but there is still control using , if bSaveIntermediate=false it has priority over m_SaveFilesLevel
      if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
         sprintf(outDirtyImageReal,"%s/dirty_test_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),width,height);
         sprintf(outDirtyImageImag,"%s/dirty_test_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),width,height);
   
         m_pSkyImageRealTmp->WriteFits( outDirtyImageReal );
         m_pSkyImageImagTmp->WriteFits( outDirtyImageImag );
      }
   }
   
   // 2022-04-02 : test change to use member variable for final image (have to be careful with threads and to not use this class as global variable):
   // calculate and save FFT-shifted image :
   // CBgFits out_image_real2( out_image_real.GetXSize(), out_image_real.GetYSize() ), out_image_imag2( out_image_real.GetXSize(), out_image_real.GetYSize() );
   // AllocOutPutImages( m_pSkyImageRealTmp->GetXSize(), m_pSkyImageRealTmp->GetYSize() );
   
   if( !m_pSkyImageReal || !m_pSkyImageImag ){
      printf("ERROR in code : internal image buffers not allocated -> cannot continue\n");
      return;
   }
   if( !m_pSkyImageRealTmp || !m_pSkyImageImagTmp ){
      printf("ERROR in code : internal temporary image buffers m_pSkyImageRealTmp and/or m_pSkyImageImagTmp not allocated -> cannot continue\n");
      return;
   }

   fft_shift( *m_pSkyImageRealTmp, *m_pSkyImageReal );
   fft_shift( *m_pSkyImageImagTmp, *m_pSkyImageImag );
   
   int rest = 1; // just so that by default it is !=0 -> image not saved 
   if( CPacerImager::m_SaveControlImageEveryNth > 0 ){
      rest = (m_SkyImageCounter % CPacerImager::m_SaveControlImageEveryNth);
      if( rest == 0 ){
          PRINTF_INFO("INFO : saving %d-th control sky image\n",m_SkyImageCounter);
      }
   }

   if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_FINAL || rest==0 ){   
      if( szBaseOutFitsName && strlen(szBaseOutFitsName) ){
         sprintf(outDirtyImageReal,"%s/%s_real.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
      }else{
         // sprintf(outDirtyImageReal,"dirty_test_real_fftshift_%dx%d.fits",width,height);
         // const char* get_filename(  time_t ut_time , char* out_buffer, int usec=0, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix="", const char* formater="%.2u%.2u%.2uT%.2u%.2u%.2u" );
         get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageReal, m_ImagerParameters.m_szOutputDirectory.c_str(), "dirty_image_", "_real" ); // uxtime=0 -> it will be taken as current system time
      }
      SaveSkyImage( outDirtyImageReal , m_pSkyImageReal, m_ImagerParameters.m_fUnixTime );
      PRINTF_DEBUG("Saved read file to %s\n",outDirtyImageReal);
   
      if( bSaveImaginary ){
         if( szBaseOutFitsName && strlen(szBaseOutFitsName) ){
            sprintf(outDirtyImageImag,"%s/%s_imag.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
         }else{
            // sprintf(outDirtyImageImag,"dirty_test_imag_fftshift_%dx%d.fits",width,height);
            get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageImag, m_ImagerParameters.m_szOutputDirectory.c_str(), "dirty_image_", "_imag" );
         }

         m_pSkyImageImag->SetFileName( outDirtyImageImag );      
         m_pSkyImageImag->WriteFits( outDirtyImageImag );
         PRINTF_DEBUG("Saved imaginary file to %s\n",outDirtyImageImag);
      }
   }

   // free memory :
   // fftw_free( in_buffer ); 
   // fftw_free( out_buffer );
   
   if( CPacerImager::m_bPrintImageStatistics ){
      double mean, rms, minval, maxval, median, iqr, rmsiqr;
      int cnt;
      int radius = int( sqrt( m_pSkyImageReal->GetXSize()*m_pSkyImageReal->GetXSize() + m_pSkyImageReal->GetYSize()*m_pSkyImageReal->GetYSize() ) ) + 10;
      // m_SkyImageReal.GetStat( mean, rms, minval, maxval );
      m_pSkyImageReal->GetStatRadiusAll( mean, rms, minval, maxval, median, iqr, rmsiqr, cnt, radius, true );
      printf("STAT : full image %s statistics in radius = %d around the center using %d pixels : mean = %.6f , rms = %.6f, minval = %.6f, maxval = %.6f, median = %.6f, rms_iqr = %.6f\n",outDirtyImageReal,radius,cnt,mean, rms, minval, maxval, median, rmsiqr );
      
      
      // TODO : this will be parameterised to specified requested windown in the image to get RMS value from:
      double mean_window, rms_window, minval_window, maxval_window, median_window, iqr_window, rmsiqr_window;
      radius = 10; // TODO : make it use parameter and also position in the image 
      m_pSkyImageReal->GetStatRadiusAll( mean_window, rms_window, minval_window, maxval_window, median_window, iqr_window, rmsiqr_window, cnt, radius, true );
      printf("STAT : statistics of %s in radius = %d around the center using %d pixels : mean = %.6f , rms = %.6f, minval = %.6f, maxval = %.6f, median = %.6f, rms_iqr = %.6f\n",outDirtyImageReal,radius,cnt,mean_window, rms_window, minval_window, maxval_window, median_window, rmsiqr_window );
   }
   
   // TODO : re-grid to SKY COORDINATES !!!
   // convert cos(alpha) to alpha - see notes !!!
   // how to do it ???

   PACER_PROFILER_END("dirty imaging took")     
}

bool CPacerImager::CalculateUVW( bool bForce /*=false*/, bool bInitialise /*=true*/ )
{
   if( bInitialise ){ // is not required in one case of call from ::Initialise itself -> to avoid recursive call 
      Initialise();
   }
 
   bool bRecalculationRequired = false;
   
   if( m_Baselines <=0 || m_U.GetXSize() <= 0 || m_V.GetXSize() <= 0 || m_W.GetXSize() <= 0 || !m_ImagerParameters.m_bConstantUVW || bForce ){
      bRecalculationRequired = true;      
   }
   
   if( bRecalculationRequired ){
      PRINTF_DEBUG("DEBUG : recalculation of UVW is required\n");
      
      m_Baselines = m_MetaData.m_AntennaPositions.CalculateUVW( m_U, m_V, m_W, (CPacerImager::m_SaveFilesLevel>=SAVE_FILES_DEBUG), m_ImagerParameters.m_szOutputDirectory.c_str(), m_bIncludeAutos );
      PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines\n",m_Baselines); 
      
      UpdateParameters(); // update parameters after change in UVW 
   }
   
   return (m_Baselines>0);
}

bool CPacerImager::UpdateParameters()
{
   PRINTF_DEBUG("DEBUG : updating parameters (pixscale etc.) based on new UVW\n");
   // U :   
   m_U.GetStat( u_mean, u_rms, u_min, u_max );
 
   // V : 
   m_V.GetStat( v_mean, v_rms, v_min, v_max );
 
   // W : 
   m_W.GetStat( w_mean, w_rms, w_min, w_max );
 
   // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
   // was -35 / +35 
   u_min = -u_max;
   //  u_max = +35;  
   v_min = -v_max;
   //  v_max = +35;

   // recalculate PIXSCALE at zenith :
   double frequency_hz = m_ImagerParameters.m_fCenterFrequencyMHz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;      
   double u_max_lambda = u_max/wavelength_m;
   double pixscale_radians = 1.00/(2.00*u_max_lambda);
   m_PixscaleAtZenith = pixscale_radians*(180.00/M_PI);
   PRINTF_INFO("INFO : UVW updated UV range (%.6f,%.6f) - (%.6f,%.6f), pixscale at zenith calculated as %.6f [deg] ( = %.3f [arcmin] ) for frequency_hz = %.8f Hz\n",u_min/wavelength_m,v_min/wavelength_m,u_max/wavelength_m,v_max/wavelength_m,m_PixscaleAtZenith,m_PixscaleAtZenith*60.00,frequency_hz);

   return true;
}

bool CPacerImager::ReadOrCalcUVW( const char* basename, const char* szPostfix )
{
  string fits_file_u = basename;
  fits_file_u += "_u";
  if( strlen( szPostfix ) ){
     fits_file_u += szPostfix;
  }
  fits_file_u += ".fits";

  string fits_file_v = basename;
  fits_file_v += "_v";
  if( strlen( szPostfix ) ){
     fits_file_v += szPostfix;
  }
  fits_file_v += ".fits";

  string fits_file_w = basename;
  fits_file_w += "_w";
  if( strlen( szPostfix ) ){
     fits_file_w += szPostfix;
  }
  fits_file_w += ".fits";

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){  
     printf("DEBUG : Expecting the following files to exist:\n");
     printf("\t%s\n",fits_file_u.c_str()); 
     printf("\t%s\n",fits_file_v.c_str()); 
     printf("\t%s\n",fits_file_w.c_str()); 
  }


  int n_ant = m_MetaData.m_AntennaPositions.size();
  bool bCalculateMetaFits = false;
  if( n_ant>0 && (m_MetaData.HasMetaFits() || m_ImagerParameters.m_bConstantUVW ) ){
//  if( n_ant > 0 ){
     printf("DEBUG : can calculate UVW n_ant = %d , has_metadata = %d , constant UVW = %d\n",n_ant,m_MetaData.HasMetaFits(),m_ImagerParameters.m_bConstantUVW);
     bCalculateMetaFits = true;
  }else{
     printf("DEBUG : cannot calculate UVW n_ant = %d , has_metadata = %d , constant UVW = %d\n",n_ant,m_MetaData.HasMetaFits(),m_ImagerParameters.m_bConstantUVW);
  }

  if( bCalculateMetaFits ){
     if( !CalculateUVW() ){
        printf("ERROR : could not calculate UVW coordinates\n");
        return false;
     }
  }else{
     PRINTF_WARNING("WARNING : antenna position file %s not specified or does not exist -> will try using UVW FITS files : %s,%s,%s\n",m_ImagerParameters.m_AntennaPositionsFile.c_str(),fits_file_u.c_str(),fits_file_v.c_str(),fits_file_w.c_str());

     // U :
     PRINTF_INFO("Reading fits file %s ...\n",fits_file_u.c_str());
     if( m_U.ReadFits( fits_file_u.c_str(), 0, 1, 1 ) ){
        printf("ERROR : could not read U FITS file %s\n",fits_file_u.c_str());
        return false;
     }else{
        PRINTF_INFO("OK : fits file %s read ok\n",fits_file_u.c_str());
     }
  
     // V : 
     PRINTF_INFO("Reading fits file %s ...\n",fits_file_v.c_str());
     if( m_V.ReadFits( fits_file_v.c_str(), 0, 1, 1 ) ){
        printf("ERROR : could not read V FITS file %s\n",fits_file_v.c_str());
        return false;
     }else{
        PRINTF_INFO("OK : fits file %s read ok\n",fits_file_v.c_str());
     }
  
     // W : 
     PRINTF_INFO("Reading fits file %s ...\n",fits_file_w.c_str());
     if( m_W.ReadFits( fits_file_w.c_str(), 0, 1, 1 ) ){
        printf("ERROR : could not read W FITS file %s\n",fits_file_w.c_str());
        return false;
     }else{
        PRINTF_INFO("OK : fits file %s read ok\n",fits_file_w.c_str());
     }
  }

  return true;
}

bool CPacerImager::read_corr_matrix( const char* basename, CBgFits& fits_vis_real, CBgFits& fits_vis_imag, 
                                     const char* szPostfix )
{
  // ensures initalisation of object structures 
  Initialise();
  
  // creating FITS file names for REAL, IMAG and U,V,W input FITS files :
  string fits_file_real = basename;
  fits_file_real += "_vis_real";
  if( strlen( szPostfix ) ){
     fits_file_real += szPostfix;
  }
  fits_file_real += ".fits";
 
  string fits_file_imag = basename;
  fits_file_imag += "_vis_imag";
  if( strlen( szPostfix ) ){
     fits_file_imag += szPostfix;
  }
  fits_file_imag += ".fits";
 

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
  
  // Test horizontal flip :
/*  fits_vis_real.HorFlip();
  fits_vis_imag.HorFlip();
  fits_vis_real.WriteFits("re_hor_flip.fits");
  fits_vis_imag.WriteFits("im_hor_flip.fits");*/
  
  // TEST Conjugate :
  // Conjugate is equivalent to changing FFTW BACKWARD TO FORWARD 
/*  for(int y=0;y<fits_vis_imag.GetYSize();y++){
     for(int x=0;x<fits_vis_imag.GetXSize();x++){
        if( x!=y ){
           double val = fits_vis_imag.getXY(x,y);
           fits_vis_imag.setXY(x,y,-val);
        }
     }
  }*/

  // reads or calculates UVW coordinates
  bool ret = ReadOrCalcUVW( basename , szPostfix );
  return ret;  
}

void CPacerImager::gridding( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
               CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
               double frequency_mhz, 
               int n_pixels,
               double min_uv /*=-1000*/,
               const char* weighting /*="" weighting : U for uniform (others not implemented) */
            )
{
  // debug :
  // fits_vis_real.WriteFits("fits/test_re.fits");
  // fits_vis_imag.WriteFits("fits/test_im.fits");
  // fits_vis_u.WriteFits("fits/test_u.fits");
  PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f\n",min_uv);

  PACER_PROFILER_START
  bool bStatisticsCalculated = false;
  if( m_bPrintImageStatistics ){ // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
     fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );
  
     // V : 
     fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );
  
     // W : 
     fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );
     

     // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
     // was -35 / +35 
     u_min = -u_max;
     //  u_max = +35;  
     v_min = -v_max;
     //  v_max = +35;
     
     bStatisticsCalculated = true;
  }


// calculate using CASA formula from image_tile_auto.py :
// synthesized_beam=(lambda_m/max_baseline)*(180.00/math.pi)*60.00 # in arcmin
// lower=synthesized_beam/5.00
// higher=synthesized_beam/3.00
// cellside_float=(lower+higher)*0.5
// NEW :
//   double alpha = 4.00/15.00; // from (1/5+1/3)*0.5 = 4/15
//   double wrong_factor = 1.00; // 4.00 factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_u = wrong_factor*alpha*(u_max-u_min)/(n_pixels); // factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_v = wrong_factor*alpha*(v_max-v_min)/(n_pixels);
//   int freq_channel = 204;
//   double frequency_mhz = freq_channel*(400.00/512.00);
//   if( gFrequencyMHz > 0 ){
//      frequency_mhz = gFrequencyMHz;
//   }
   double frequency_hz = frequency_mhz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;
   // UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c
   // delta_u = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   // delta_v = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
// OLD :
//  double delta_u = (u_max - u_min) / n_pixels;  
//  double delta_v = (v_max - v_min) / n_pixels;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){  
     // see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220909_image_sizes_verification.odt
     // TODO : try this in the future or so
     // it should rathr be : 1.00/(2 U_max) = Lambda/(2 MAX_baseline) :
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg] [u_max = %.8f , delta_u = %.8f ] (using OLD)\n",pixscale_zenith_deg,pixscale_deg_version2,u_max,delta_u);
// 2024-03-16 - after pixscale param added     m_PixscaleAtZenith = pixscale_deg_version2; // pixscale_zenith_deg;
     if( m_bCompareToMiriad && false  ){ // 2023-12-17 temporary disabled as WCS is not correctly saved then !!! see 20231215_repeat_processing_of_small_part_of_20230709.odt
        PRINTF_WARNING("WARNING : MIRIAD-like option -> pixscale set to %.8f\n",pixscale_zenith_deg);
        m_PixscaleAtZenith = pixscale_zenith_deg;
     }
     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  // is it ok to chose the UV plane center based on this:  
//  double u_center = (u_min + u_max)/2.00;  
//  double v_center = (v_min + v_max)/2.00;  

  // Limits of UVW :
  // double GetStat( double& mean, double& rms, double& minval, double& maxval, 


  // simple gridding :
  uv_grid_real.SetValue( 0.00 );
  uv_grid_imag.SetValue( 0.00 );
  uv_grid_counter.SetValue( 0.00 );
      
  int added=0, high_value=0;
  for( int ant1 = 0 ; ant1 < fits_vis_real.GetXSize(); ant1++ ){
     for( int ant2 = 0 ; ant2 < fits_vis_real.GetXSize(); ant2++ ){
        if( ant1 > ant2 || (m_bIncludeAutos && ant1==ant2) ){ // was ant1 > ant2 
           if( m_MetaData.m_AntennaPositions.size() > 0 ){
              if( m_MetaData.m_AntennaPositions[ant1].flag > 0 || m_MetaData.m_AntennaPositions[ant2].flag > 0 ){
                 // skip flagged antennas 
                 // printf("Fast flagging used\n");
                 continue;
              }
           }else{
              if( m_FlaggedAntennas.size() > 0 ){
                 // WARNING : this is much less efficient so better to have antenna positions and check flags there
                 if( find_value( m_FlaggedAntennas, ant1 ) >=0 || find_value( m_FlaggedAntennas, ant2 ) >=0 ){
                    continue;
                 }
              }
           }                     
        
           double re = fits_vis_real.getXY(ant1,ant2);
           double im = fits_vis_imag.getXY(ant1,ant2);
           
           if( ant1==ant2 ){
              PRINTF_DEBUG("Auto-correlation debug values %.4f / %.4f\n",re,im);
           }
           
           if( !isnan(re) && !isnan(im) ){
              if ( fabs(re) < MAX_VIS && fabs(im) < MAX_VIS ){
                 // TODO convert [m] -> wavelength 
                 double u = fits_vis_u.getXY(ant1,ant2) / wavelength_m;
                 
                 // 2022-09-24 : - removed for a test on MWA data
                 // 2022-09-09 - for now sticking to - sign here to have back comatible test EDA2 data
                 //  but looks like this (-) should not be here at least does not match UV coverage from WSCEAN (it is flipped then see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                 double v = fits_vis_v.getXY(ant1,ant2) / wavelength_m; // the - sign here fixes the Y flip, but I am not sure why needed ??? check RTS : imagefromuv.c , LM_CopyFromFFT
                                                                        // where some interesting re-shuffling happens too !
                                                                        // - was for EDA2 data, but + is ok for MWA data - it may be consistent once I start using TMS equation 4.1 consistently for 
                                                                        // both EDA2 and MWA  
                 double w = fits_vis_w.getXY(ant1,ant2) / wavelength_m;
                 double uv_distance = sqrt(u*u + v*v);
                 
/* this is in WSCLEAN, but here seems to have no effect ...
                 if (w < 0.0 ) { // && !_isComplex
                    u = -u;
                    v = -v;
                    w = -w;
                    im = -im;
                 }
*/ 
                 
                 if( ant1==ant2 ){
                    PRINTF_DEBUG("Auto-correlation debug2 values %.4f / %.4f , uv_distance = %.8f vs. min_uv = %.8f (u = %.8f , v = %.8f , wavelength_m = %.8f [m])\n",re,im,uv_distance,min_uv,u,v,wavelength_m);
                 }
              
                 if( uv_distance > min_uv ){ // check if larger than minimum UV distance 
//                 int u_index = round( (u - u_min)/delta_u );
//                 int v_index = round( (v - v_min)/delta_v );
                    double u_pix = round( u/delta_u );
                    double v_pix = round( v/delta_v );
                    int u_index = u_pix + n_pixels/2; // was u - u_center
                    int v_index = v_pix + n_pixels/2; // was v - v_center
                    
              
                 // Using CELL averaging method or setXY ?
                    int ret = uv_grid_real.addXY( u_index, v_index, re );
                    if( ret < 0 ){
                       printf("WARNING : (u_index,v_index) = (%d,%d) outside the grid size %d x %d\n",u_index, v_index, n_pixels,n_pixels);
                    }
                    uv_grid_imag.addXY( u_index, v_index, im );
                    uv_grid_counter.addXY( u_index, v_index ,  1.00 );
                    
                    if( m_ImagerDebugLevel > 101 ){
                         printf("DEBUG ADDING : PACER_gridding at (x,y) =  %.1f %.1f += %.8f + j %.8f (u,v,w) = ( %.4f , %.4f , %.4f )\n",u_pix,v_pix,re,im,u,v,w);
//                       printf("DEBUG gridding : (u,v) = (%.8f,%.8f) -> index = (%d,%d) or (%.2f,%.2f), value = %.8f + j %.8f\n",u,v,u_index,v_index,u_pix,v_pix,re,im);
                    }
              
                 // add conjugates :
//                 u_index = round( (-u - u_min)/delta_u );
//                 v_index = round( (-v - v_min)/delta_v );
                    u_index = -u_pix + n_pixels/2; // was round( (-u - u_center)/delta_u ) + ...
                    v_index = -v_pix + n_pixels/2; // was round( (-v - v_center)/delta_v ) + ...
                    
                    ret = uv_grid_real.addXY( u_index, v_index, re );
                    if( ret < 0 ){
                       printf("WARNING : (u_index,v_index) = (%d,%d) outside the grid size %d x %d\n",u_index, v_index, n_pixels,n_pixels);
                    }
                    uv_grid_imag.addXY( u_index, v_index, -im );
                    uv_grid_counter.addXY( u_index, v_index ,  1.00 );
                    
                    
                    if( ant1==ant2 ){
                       PRINTF_DEBUG("Auto-correlation added values %.4f / %.4f\n",re,im);
                    }
                                                                           
                    added++;
                 }
              }else{
                 PRINTF_DEBUG("DEBUG : visibility value %e +j%e higher than limit %e -> skipped\n",re,im,MAX_VIS);
                 high_value++;
              }
           }
        }
     }
  }  
  
  // This division is in fact UNIFORM weighting !!!! Not CELL-avareging 
  // normalisation to make it indeed CELL-averaging :
  if( strcmp(weighting, "U" ) == 0 ){
     uv_grid_real.Divide( uv_grid_counter );
     uv_grid_imag.Divide( uv_grid_counter );
  }  
  PACER_PROFILER_END("gridding (no I/O) took")
  PRINTF_DEBUG("DEBUG : added %d UV points to the grid, %d too high values skipped\n",added,high_value);
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
     char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
     sprintf(uv_grid_re_name,"%s/uv_grid_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_im_name,"%s/uv_grid_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_counter_name,"%s/uv_grid_counter_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
    
     if( uv_grid_real.WriteFits( uv_grid_re_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_re_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_re_name);
     }

     if( uv_grid_imag.WriteFits( uv_grid_im_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_im_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_im_name);
     }
  
     if( uv_grid_counter.WriteFits( uv_grid_counter_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_counter_name);
     }
  }

}

// TODO : use ConvertFits2XCorr to convert fits_vis_real,fits_vis_imag into Visibitilies-xcorr and use gridding_fast( Visibitilies& xcorr ... ) function to avoid code duplication
void CPacerImager::gridding_fast( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
               CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
               double frequency_mhz, 
               int n_pixels,
               double min_uv /*=-1000*/,
               const char* weighting /*="" weighting : U for uniform (others not implemented) */
            )
{
  // debug :
  // fits_vis_real.WriteFits("fits/test_re.fits");
  // fits_vis_imag.WriteFits("fits/test_im.fits");
  // fits_vis_u.WriteFits("fits/test_u.fits");
  PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f, IncludeAutos = %d , flagged antennas = %d , metadata.ant.size = %d\n",min_uv,m_bIncludeAutos,int(m_FlaggedAntennas.size()),int(m_MetaData.m_AntennaPositions.size()));

  PACER_PROFILER_START
  bool bStatisticsCalculated = false;

  if( m_bPrintImageStatistics ){ // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
     fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );
  
     // V : 
     fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );
  
     // W : 
     fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );

     // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
     // was -35 / +35 
     u_min = -u_max;
     //  u_max = +35;  
     v_min = -v_max;
     //  v_max = +35;
     
     bStatisticsCalculated = true;
  }   


// calculate using CASA formula from image_tile_auto.py :
// synthesized_beam=(lambda_m/max_baseline)*(180.00/math.pi)*60.00 # in arcmin
// lower=synthesized_beam/5.00
// higher=synthesized_beam/3.00
// cellside_float=(lower+higher)*0.5
// NEW :
//   double alpha = 4.00/15.00; // from (1/5+1/3)*0.5 = 4/15
//   double wrong_factor = 1.00; // 4.00 factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_u = wrong_factor*alpha*(u_max-u_min)/(n_pixels); // factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_v = wrong_factor*alpha*(v_max-v_min)/(n_pixels);
//   int freq_channel = 204;
//   double frequency_mhz = freq_channel*(400.00/512.00);
//   if( gFrequencyMHz > 0 ){
//      frequency_mhz = gFrequencyMHz;
//   }
   double frequency_hz = frequency_mhz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;
   // UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c
   // delta_u = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   // delta_v = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
// OLD :
//  double delta_u = (u_max - u_min) / n_pixels;  
//  double delta_v = (v_max - v_min) / n_pixels;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){  
     // see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220909_image_sizes_verification.odt
     // TODO : try this in the future or so
     // it should rathr be : 1.00/(2 U_max) = Lambda/(2 MAX_baseline) :
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg] [ u_max = %.6f , delta_u = %.8f , n_pixels = %d ] (using OLD)\n",pixscale_zenith_deg,pixscale_deg_version2,u_max,delta_u,n_pixels);
// 2024-03-16 - after pixscale param added     m_PixscaleAtZenith = pixscale_deg_version2; // pixscale_zenith_deg;
     if( m_bCompareToMiriad && false ){ // 2023-12-17 temporary disabled as WCS is not correctly saved then !!! see 20231215_repeat_processing_of_small_part_of_20230709.odt
        PRINTF_WARNING("WARNING : MIRIAD-like option -> pixscale set to %.8f\n",pixscale_zenith_deg);
        m_PixscaleAtZenith = pixscale_zenith_deg;
     }

     
     /*double pixscale_check_rad = (1/delta_u)/n_pixels;
     printf("CHECK : pixscale_check_deg = %.6f [deg] vs. %.6f [deg]\n",pixscale_check_rad*(180.00/M_PI),(1/(2.00*u_max))*(180.00/M_PI));
     m_PixscaleAtZenith = ((180.00)/n_pixels);*/
     
     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  // is it ok to chose the UV plane center based on this:  
//  double u_center = (u_min + u_max)/2.00;  
//  double v_center = (v_min + v_max)/2.00;  

  // Limits of UVW :
  // double GetStat( double& mean, double& rms, double& minval, double& maxval, 

  int center_x = int( n_pixels / 2 );
  int center_y = int( n_pixels / 2 );
  int is_odd_x = 0 , is_odd_y = 0;
  if( (n_pixels % 2) == 1 ){
     is_odd_x = 1;
  }
  if( (n_pixels % 2) == 1 ){
     is_odd_y = 1;
  }


  // simple gridding :
  uv_grid_real.SetValue( 0.00 );
  uv_grid_imag.SetValue( 0.00 );
  uv_grid_counter.SetValue( 0.00 );
      
  int added=0, high_value=0;
  for( int ant1 = 0 ; ant1 < fits_vis_real.GetXSize(); ant1++ ){
     for( int ant2 = 0 ; ant2 < fits_vis_real.GetXSize(); ant2++ ){
        if( ant1 > ant2 || (m_bIncludeAutos && ant1==ant2) ){ // was ant1 > ant2 
           if( m_MetaData.m_AntennaPositions.size() > 0 ){
              if( m_MetaData.m_AntennaPositions[ant1].flag > 0 || m_MetaData.m_AntennaPositions[ant2].flag > 0 ){
                 // skip flagged antennas 
                 // printf("Fast flagging used\n");
                 printf("DEBUG : baseline %d-%d flagged\n",ant1,ant2);
                 continue;
              }
           }else{
              if( m_FlaggedAntennas.size() > 0 ){
                 // WARNING : this is much less efficient so better to have antenna positions and check flags there
                 if( find_value( m_FlaggedAntennas, ant1 ) >=0 || find_value( m_FlaggedAntennas, ant2 ) >=0 ){
                    printf("DEBUG : baseline %d-%d flagged [2]\n",ant1,ant2);
                    continue;
                 }
              }
           }                     
        
           double re = fits_vis_real.getXY(ant1,ant2);
           double im = fits_vis_imag.getXY(ant1,ant2);
           
           if( ant1==ant2 ){
              PRINTF_DEBUG("Auto-correlation debug values %.4f / %.4f\n",re,im);
           }
           
           if( !isnan(re) && !isnan(im) ){
              if ( fabs(re) < MAX_VIS && fabs(im) < MAX_VIS ){
                 // TODO convert [m] -> wavelength 
                 double u = fits_vis_u.getXY(ant1,ant2) / wavelength_m;
                 
                 // 2022-09-24 : - removed for a test on MWA data
                 // 2022-09-09 - for now sticking to - sign here to have back comatible test EDA2 data
                 //  but looks like this (-) should not be here at least does not match UV coverage from WSCEAN (it is flipped then see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                 double v = fits_vis_v.getXY(ant1,ant2) / wavelength_m; // the - sign here fixes the Y flip, but I am not sure why needed ??? check RTS : imagefromuv.c , LM_CopyFromFFT
                                                                        // where some interesting re-shuffling happens too !
                                                                        // - was for EDA2 data, but + is ok for MWA data - it may be consistent once I start using TMS equation 4.1 consistently for 
                                                                        // both EDA2 and MWA  
                 double w = fits_vis_w.getXY(ant1,ant2) / wavelength_m;
                 double uv_distance = sqrt(u*u + v*v);
                 
                 // check W values against limits :
                 if( fabs(w) < m_ImagerParameters.m_MinW || fabs(w) > m_ImagerParameters.m_MaxW ){
                    // if w is outside limits -> skip visibility point :
                    continue;
                 }
                 
                 
/* this is in WSCLEAN, but here seems to have no effect ...
                 if (w < 0.0 ) { // && !_isComplex
                    u = -u;
                    v = -v;
                    w = -w;
                    im = -im;
                 }
*/ 
                 
                 if( ant1==ant2 ){
                    PRINTF_DEBUG("Auto-correlation debug2 values %.4f / %.4f , uv_distance = %.8f vs. min_uv = %.8f (u = %.8f , v = %.8f , wavelength_m = %.8f [m])\n",re,im,uv_distance,min_uv,u,v,wavelength_m);
                 }
              
                 if( uv_distance > min_uv ){ // check if larger than minimum UV distance 
//                 int u_index = round( (u - u_min)/delta_u );
//                 int v_index = round( (v - v_min)/delta_v );
                    double u_pix = round( u/delta_u );
                    double v_pix = round( v/delta_v );
                    int u_index = u_pix + n_pixels/2; // was u - u_center
                    int v_index = v_pix + n_pixels/2; // was v - v_center                                       
                    
                    // calculate position (x,y) on UV grid as expected by fftw (no fft_unshift required)
                    int x_grid = 0, y_grid = 0;
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }
                    
                    
                 
//                 printf("DEBUG : u_index %d vs. %d ( u = %.2f , u_min = %.2f , delta_u = %.2f , u_center =%.2f)\n",u,u_index1,u_index,u_min,delta_u,u_center);
              
                 // Using CELL averaging method or setXY ?
                    uv_grid_real.addXY( x_grid, y_grid, re );
                    uv_grid_imag.addXY( x_grid, y_grid, im );
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );
                    
                    if( m_ImagerDebugLevel > 101 ){
                         if( x_grid <= 5 && y_grid <= 5 ){
                            printf("DEBUG ADDING : PACER_gridding [1] at (x_grid,y_grid) = (%d,%d) (u_pix,v_pix) = ( %.1f %.1f ) += %.8f + j %.8f = %.8f + j ??? , (u,v,w) = ( %.4f , %.4f , %.4f ), ants %d-%d\n",x_grid,y_grid,u_pix,v_pix,re,im,uv_grid_real.getXY(x_grid,y_grid),u,v,w,ant1,ant2);
                         }
//                       printf("DEBUG gridding : (u,v) = (%.8f,%.8f) -> index = (%d,%d) or (%.2f,%.2f), value = %.8f + j %.8f\n",u,v,u_index,v_index,u_pix,v_pix,re,im);
                    }
              
                 // add conjugates :
//                 u_index = round( (-u - u_min)/delta_u );
//                 v_index = round( (-v - v_min)/delta_v );
                    u_index = -u_pix + n_pixels/2; // was round( (-u - u_center)/delta_u ) + ...
                    v_index = -v_pix + n_pixels/2; // was round( (-v - v_center)/delta_v ) + ...
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }                                       
                    
                    uv_grid_real.addXY( x_grid, y_grid, re );
                    uv_grid_imag.addXY( x_grid, y_grid, -im );
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );

                    if( m_ImagerDebugLevel > 101 ){
                         if( x_grid <= 5 && y_grid <= 5 ){
                            printf("DEBUG ADDING : PACER_gridding [2] at (x_grid,y_grid) = (%d,%d) (u_pix,v_pix) = ( %.1f %.1f ) += %.8f + j %.8f = %.8f + j ??? , (u,v,w) = ( %.4f , %.4f , %.4f ), ants %d-%d\n",x_grid,y_grid,u_pix,v_pix,re,-im,uv_grid_real.getXY(x_grid,y_grid),u,v,w,ant1,ant2);
                         }
//                       printf("DEBUG gridding : (u,v) = (%.8f,%.8f) -> index = (%d,%d) or (%.2f,%.2f), value = %.8f + j %.8f\n",u,v,u_index,v_index,u_pix,v_pix,re,im);
                    }
                    
                    if( ant1==ant2 ){
                       PRINTF_DEBUG("Auto-correlation added values %.4f / %.4f\n",re,im);
                    }
                    
                    added++;
                 }
              }else{
                 PRINTF_DEBUG("DEBUG : visibility value %e +j%e higher than limit %e -> skipped\n",re,im,MAX_VIS);
                 high_value++;
              }
           }else{
              PRINTF_DEBUG("DEBUG : nan_detected for visibility from antenna pair (%d,%d)\n",ant1,ant2);
           }
        }
     }
  }  
  
  // This division is in fact UNIFORM weighting !!!! Not CELL-avareging 
  // normalisation to make it indeed CELL-averaging :
  if( strcmp(weighting, "U" ) == 0 ){
     uv_grid_real.Divide( uv_grid_counter );
     uv_grid_imag.Divide( uv_grid_counter );
  }  
  PACER_PROFILER_END("gridding (no I/O) took")
  PRINTF_DEBUG("DEBUG : added %d UV points to the grid, %d too high values skipped\n",added,high_value);
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
     char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
     sprintf(uv_grid_re_name,"%s/uv_grid_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_im_name,"%s/uv_grid_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_counter_name,"%s/uv_grid_counter_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     
     printf("DEBUG gridding UV_grid_real(4,0) = %.20f\n",uv_grid_real.getXY(4,0));
    
     if( uv_grid_real.WriteFits( uv_grid_re_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_re_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_re_name);
     }

     if( uv_grid_imag.WriteFits( uv_grid_im_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_im_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_im_name);
     }
  
     if( uv_grid_counter.WriteFits( uv_grid_counter_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_counter_name);
     }
  }

}

void CPacerImager::gridding_fast( Visibilities& xcorr, 
               int time_step, 
               int fine_channel,
               CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
               CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
               double frequency_mhz, 
               int n_pixels,
               double min_uv /*=-1000*/,
               const char* weighting /*="" weighting : U for uniform (others not implemented) */
            )
{
  // debug :
  // fits_vis_real.WriteFits("fits/test_re.fits");
  // fits_vis_imag.WriteFits("fits/test_im.fits");
  // fits_vis_u.WriteFits("fits/test_u.fits");
  PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f\n",min_uv);

  PACER_PROFILER_START
  bool bStatisticsCalculated = false;
  if( m_bPrintImageStatistics ){ // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
     fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );
  
     // V : 
     fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );
  
     // W : 
     fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );

     // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
     // was -35 / +35 
     u_min = -u_max;
     //  u_max = +35;  
     v_min = -v_max;
     //  v_max = +35;
     
     bStatisticsCalculated = true;
  }


// calculate using CASA formula from image_tile_auto.py :
// synthesized_beam=(lambda_m/max_baseline)*(180.00/math.pi)*60.00 # in arcmin
// lower=synthesized_beam/5.00
// higher=synthesized_beam/3.00
// cellside_float=(lower+higher)*0.5
// NEW :
//   double alpha = 4.00/15.00; // from (1/5+1/3)*0.5 = 4/15
//   double wrong_factor = 1.00; // 4.00 factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_u = wrong_factor*alpha*(u_max-u_min)/(n_pixels); // factor for due to U range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 )
//   double delta_v = wrong_factor*alpha*(v_max-v_min)/(n_pixels);
//   int freq_channel = 204;
//   double frequency_mhz = freq_channel*(400.00/512.00);
//   if( gFrequencyMHz > 0 ){
//      frequency_mhz = gFrequencyMHz;
//   }
   double frequency_hz = frequency_mhz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;
   // UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c
   // delta_u = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   // delta_v = (VEL_LIGHT/frequency_Hz)/(gFOV_degrees*M_PI/180.);
   PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);
// OLD :
//  double delta_u = (u_max - u_min) / n_pixels;  
//  double delta_v = (v_max - v_min) / n_pixels;

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){  
     // see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220909_image_sizes_verification.odt
     // TODO : try this in the future or so
     // it should rathr be : 1.00/(2 U_max) = Lambda/(2 MAX_baseline) :
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg] , u_max = %.6f (using OLD)\n",pixscale_zenith_deg,pixscale_deg_version2,u_max);
 // 2024-03-16 - after pixscale param added   m_PixscaleAtZenith = pixscale_deg_version2; // pixscale_zenith_deg;
     if( m_bCompareToMiriad && false ){ // 2023-12-17 temporary disabled as WCS is not correctly saved then !!! see 20231215_repeat_processing_of_small_part_of_20230709.odt
        PRINTF_WARNING("WARNING : MIRIAD-like option -> pixscale set to %.8f\n",pixscale_zenith_deg);
        m_PixscaleAtZenith = pixscale_zenith_deg;
     }

     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  // is it ok to chose the UV plane center based on this:  
//  double u_center = (u_min + u_max)/2.00;  
//  double v_center = (v_min + v_max)/2.00;  

  // Limits of UVW :
  // double GetStat( double& mean, double& rms, double& minval, double& maxval, 

  int center_x = int( n_pixels / 2 );
  int center_y = int( n_pixels / 2 );
  int is_odd_x = 0 , is_odd_y = 0;
  if( (n_pixels % 2) == 1 ){
     is_odd_x = 1;
  }
  if( (n_pixels % 2) == 1 ){
     is_odd_y = 1;
  }


  // simple gridding :
  uv_grid_real.SetValue( 0.00 );
  uv_grid_imag.SetValue( 0.00 );
  uv_grid_counter.SetValue( 0.00 );

  int n_ant = xcorr.obsInfo.nAntennas;      
  int added=0, high_value=0;
  for( int ant1 = 0 ; ant1 < n_ant; ant1++ ){
     // for( int ant2 = 0 ; ant2 < n_ant; ant2++ ){
     for( int ant2 = 0 ; ant2 <= ant1; ant2++ ){
        if( ant1 > ant2 || (m_bIncludeAutos && ant1==ant2) ){ // was ant1 > ant2 
           if( m_MetaData.m_AntennaPositions.size() > 0 ){
              if( m_MetaData.m_AntennaPositions[ant1].flag > 0 || m_MetaData.m_AntennaPositions[ant2].flag > 0 ){
                 // skip flagged antennas 
                 // printf("Fast flagging used\n");
                 continue;
              }
           }else{
              if( m_FlaggedAntennas.size() > 0 ){
                 // WARNING : this is much less efficient so better to have antenna positions and check flags there
                 if( find_value( m_FlaggedAntennas, ant1 ) >=0 || find_value( m_FlaggedAntennas, ant2 ) >=0 ){
                    continue;
                 }
              }
           }                     
        
           std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, ant1, ant2 ); // was ant1, ant2 , but ant2,ant1 does not fix the orientation of the final image either ...
           double re = vis->real(); // fits_vis_real.getXY(ant1,ant2);
           double im = vis->imag(); // fits_vis_imag.getXY(ant1,ant2);
           
           if( ant1==ant2 ){
              PRINTF_DEBUG("Auto-correlation debug values %.4f / %.4f\n",re,im);
           }
           
           if( !isnan(re) && !isnan(im) ){
              if ( fabs(re) < MAX_VIS && fabs(im) < MAX_VIS ){
                 // TODO convert [m] -> wavelength 
                 double u = -fits_vis_u.getXY(ant1,ant2) / wavelength_m;
                 
                 // 2022-09-24 : - removed for a test on MWA data
                 // 2022-09-09 - for now sticking to - sign here to have back comatible test EDA2 data
                 //  but looks like this (-) should not be here at least does not match UV coverage from WSCEAN (it is flipped then see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                 double v = -fits_vis_v.getXY(ant1,ant2) / wavelength_m; // the - sign here fixes the Y flip, but I am not sure why needed ??? check RTS : imagefromuv.c , LM_CopyFromFFT
                                                                        // where some interesting re-shuffling happens too !
                                                                        // - was for EDA2 data, but + is ok for MWA data - it may be consistent once I start using TMS equation 4.1 consistently for 
                                                                        // both EDA2 and MWA  
                 double w = -fits_vis_w.getXY(ant1,ant2) / wavelength_m;
                 double uv_distance = sqrt(u*u + v*v);
                 
/* this is in WSCLEAN, but here seems to have no effect ...
                 if (w < 0.0 ) { // && !_isComplex
                    u = -u;
                    v = -v;
                    w = -w;
                    im = -im;
                 }
*/ 
                 
                 if( ant1==ant2 ){
                    PRINTF_DEBUG("Auto-correlation debug2 values %.4f / %.4f , uv_distance = %.8f vs. min_uv = %.8f (u = %.8f , v = %.8f , wavelength_m = %.8f [m])\n",re,im,uv_distance,min_uv,u,v,wavelength_m);
                 }
              
                 if( uv_distance > min_uv ){ // check if larger than minimum UV distance 
//                 int u_index = round( (u - u_min)/delta_u );
//                 int v_index = round( (v - v_min)/delta_v );
                    double u_pix = round( u/delta_u );
                    double v_pix = round( v/delta_v );
                    int u_index = u_pix + n_pixels/2; // was u - u_center
                    int v_index = v_pix + n_pixels/2; // was v - v_center                                       
                    
                    // calculate position (x,y) on UV grid as expected by fftw (no fft_unshift required)
                    int x_grid = 0, y_grid = 0;
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }
                    
                    
                 
//                 printf("DEBUG : u_index %d vs. %d ( u = %.2f , u_min = %.2f , delta_u = %.2f , u_center =%.2f)\n",u,u_index1,u_index,u_min,delta_u,u_center);
              
                 // Using CELL averaging method or setXY ?
                    uv_grid_real.addXY( x_grid, y_grid, re );
                    uv_grid_imag.addXY( x_grid, y_grid, im );
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );
                    
                    if( m_ImagerDebugLevel > 101 ){
                         printf("DEBUG ADDING : PACER_gridding at (x,y) =  %.1f %.1f += %.8f + j %.8f (u,v,w) = ( %.4f , %.4f , %.4f )\n",u_pix,v_pix,re,im,u,v,w);
//                       printf("DEBUG gridding : (u,v) = (%.8f,%.8f) -> index = (%d,%d) or (%.2f,%.2f), value = %.8f + j %.8f\n",u,v,u_index,v_index,u_pix,v_pix,re,im);
                    }
              
                 // add conjugates :
//                 u_index = round( (-u - u_min)/delta_u );
//                 v_index = round( (-v - v_min)/delta_v );
                    u_index = -u_pix + n_pixels/2; // was round( (-u - u_center)/delta_u ) + ...
                    v_index = -v_pix + n_pixels/2; // was round( (-v - v_center)/delta_v ) + ...
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }                                       
                    
                    uv_grid_real.addXY( x_grid, y_grid, re );
                    uv_grid_imag.addXY( x_grid, y_grid, -im );
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );
                    
                    if( ant1==ant2 ){
                       PRINTF_DEBUG("Auto-correlation added values %.4f / %.4f\n",re,im);
                    }
                                                                           
                    added++;
                 }
              }else{
                 PRINTF_DEBUG("DEBUG : visibility value %e +j%e higher than limit %e -> skipped\n",re,im,MAX_VIS);
                 high_value++;
              }
           }
        }
     }
  }  
  
  // This division is in fact UNIFORM weighting !!!! Not CELL-avareging 
  // normalisation to make it indeed CELL-averaging :
  if( strcmp(weighting, "U" ) == 0 ){
     uv_grid_real.Divide( uv_grid_counter );
     uv_grid_imag.Divide( uv_grid_counter );
  }  
  PACER_PROFILER_END("gridding (no I/O) took")
  PRINTF_DEBUG("DEBUG : added %d UV points to the grid, %d too high values skipped\n",added,high_value);
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
     char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
     sprintf(uv_grid_re_name,"%s/uv_grid_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_im_name,"%s/uv_grid_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_counter_name,"%s/uv_grid_counter_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
    
     if( uv_grid_real.WriteFits( uv_grid_re_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_re_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_re_name);
     }

     if( uv_grid_imag.WriteFits( uv_grid_im_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_im_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_im_name);
     }
  
     if( uv_grid_counter.WriteFits( uv_grid_counter_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_counter_name);
     }
  }

}

void CPacerImager::gridding_fast( void* visibilities_param, 
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  void* gridded_visibilities_param, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv /*=-1000*/,    // minimum UV 
                  const char* weighting /*=""*/ // weighting : U for uniform (others not implemented)
                )
{
  int n_ant = fits_vis_u.GetXSize();
  PRINTF_DEBUG("DEBUG : gridding_fast(fftw_complex) : min_uv = %.4f, n_ant = %d\n",min_uv,n_ant);
  fftw_complex* visibilities = (fftw_complex*)visibilities_param;
  fftw_complex* gridded_visibilities = (fftw_complex*)gridded_visibilities_param;

  PACER_PROFILER_START
  bool bStatisticsCalculated = false;
  if( m_bPrintImageStatistics ){ // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
     fits_vis_u.GetStat( u_mean, u_rms, u_min, u_max );
  
     // V : 
     fits_vis_v.GetStat( v_mean, v_rms, v_min, v_max );
  
     // W : 
     fits_vis_w.GetStat( w_mean, w_rms, w_min, w_max );

     // Bacause we are also including conjugates at (-u,-v) UV point in gridding u_min = -u_max and v_min = -v_max :
     // was -35 / +35 
     u_min = -u_max;
     //  u_max = +35;  
     v_min = -v_max;
     //  v_max = +35;
     
     bStatisticsCalculated = true;
  }


   double frequency_hz = frequency_mhz*1e6;
   double wavelength_m = VEL_LIGHT / frequency_hz;
   PRINTF_DEBUG("DEBUG : wavelength = %.4f [m] , frequency = %.4f [MHz]\n",wavelength_m,frequency_mhz);

  if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){  
     // see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220909_image_sizes_verification.odt
     // TODO : try this in the future or so
     // it should rathr be : 1.00/(2 U_max) = Lambda/(2 MAX_baseline) :
     double pixscale_zenith_deg = (1.00/(n_pixels*delta_u))*(180.00/M_PI); // in degrees 
     double pixscale_radians = 1.00/(2.00*u_max);
     double pixscale_deg_version2 = pixscale_radians*(180.00/M_PI);
     printf("DEBUG : pixscale old = %.8f [deg] vs. NEW = %.8f [deg] , u_max = %.6f (using OLD)\n",pixscale_zenith_deg,pixscale_deg_version2,u_max);
// 2024-03-16 - after pixscale param added     m_PixscaleAtZenith = pixscale_deg_version2; // pixscale_zenith_deg;
     if( m_bCompareToMiriad && false ){ // 2023-12-17 temporary disabled as WCS is not correctly saved then !!! see 20231215_repeat_processing_of_small_part_of_20230709.odt
        PRINTF_WARNING("WARNING : MIRIAD-like option -> pixscale set to %.8f\n",pixscale_zenith_deg);
        m_PixscaleAtZenith = pixscale_zenith_deg;
     }

     if( bStatisticsCalculated ){
        printf("DEBUG : U limits %.8f - %.8f , delta_u = %.8f -> pixscale at zenith = %.8f [deg]\n",u_min, u_max , delta_u , m_PixscaleAtZenith );
        printf("DEBUG : V limits %.8f - %.8f , delta_v = %.8f\n", v_min, v_max , delta_v );
        printf("DEBUG : W limits %.8f - %.8f\n", w_min, w_max );
     }
  }

  int center_x = int( n_pixels / 2 );
  int center_y = int( n_pixels / 2 );
  int is_odd_x = 0 , is_odd_y = 0;
  if( (n_pixels % 2) == 1 ){
     is_odd_x = 1;
  }
  if( (n_pixels % 2) == 1 ){
     is_odd_y = 1;
  }


  // simple gridding :
  uv_grid_counter.SetValue( 0.00 );

  int added=0, high_value=0;
  for( int ant1 = 0 ; ant1 < n_ant; ant1++ ){
     // for( int ant2 = 0 ; ant2 < n_ant; ant2++ ){
     for( int ant2 = 0 ; ant2 <= ant1; ant2++ ){
        if( ant1 > ant2 || (m_bIncludeAutos && ant1==ant2) ){ // was ant1 > ant2 
           if( m_MetaData.m_AntennaPositions.size() > 0 ){
              if( m_MetaData.m_AntennaPositions[ant1].flag > 0 || m_MetaData.m_AntennaPositions[ant2].flag > 0 ){
                 // skip flagged antennas 
                 // printf("Fast flagging used\n");
                 continue;
              }
           }else{
              if( m_FlaggedAntennas.size() > 0 ){
                 // WARNING : this is much less efficient so better to have antenna positions and check flags there
                 if( find_value( m_FlaggedAntennas, ant1 ) >=0 || find_value( m_FlaggedAntennas, ant2 ) >=0 ){
                    continue;
                 }
              }
           }                     

           double re = visibilities[ant2*n_ant+ant1][0];
           double im = visibilities[ant2*n_ant+ant1][1];
           
           if( ant1==ant2 ){
              PRINTF_DEBUG("Auto-correlation debug values %.4f / %.4f\n",re,im);
           }
           
           if( !isnan(re) && !isnan(im) ){
              if ( fabs(re) < MAX_VIS && fabs(im) < MAX_VIS ){
                 // TODO convert [m] -> wavelength 
                 double u = -fits_vis_u.getXY(ant1,ant2) / wavelength_m;
                 
                 // 2022-09-24 : - removed for a test on MWA data
                 // 2022-09-09 - for now sticking to - sign here to have back comatible test EDA2 data
                 //  but looks like this (-) should not be here at least does not match UV coverage from WSCEAN (it is flipped then see : /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                 double v = -fits_vis_v.getXY(ant1,ant2) / wavelength_m; // the - sign here fixes the Y flip, but I am not sure why needed ??? check RTS : imagefromuv.c , LM_CopyFromFFT
                                                                        // where some interesting re-shuffling happens too !
                                                                        // - was for EDA2 data, but + is ok for MWA data - it may be consistent once I start using TMS equation 4.1 consistently for 
                                                                        // both EDA2 and MWA  
                 double w = -fits_vis_w.getXY(ant1,ant2) / wavelength_m;
                 double uv_distance = sqrt(u*u + v*v);
                 
/* this is in WSCLEAN, but here seems to have no effect ...
                 if (w < 0.0 ) { // && !_isComplex
                    u = -u;
                    v = -v;
                    w = -w;
                    im = -im;
                 }
*/ 
                 
                 if( ant1==ant2 ){
                    PRINTF_DEBUG("Auto-correlation debug2 values %.4f / %.4f , uv_distance = %.8f vs. min_uv = %.8f (u = %.8f , v = %.8f , wavelength_m = %.8f [m])\n",re,im,uv_distance,min_uv,u,v,wavelength_m);
                 }
              
                 if( uv_distance > min_uv ){ // check if larger than minimum UV distance 
//                 int u_index = round( (u - u_min)/delta_u );
//                 int v_index = round( (v - v_min)/delta_v );
                    double u_pix = round( u/delta_u );
                    double v_pix = round( v/delta_v );
                    int u_index = u_pix + n_pixels/2; // was u - u_center
                    int v_index = v_pix + n_pixels/2; // was v - v_center                                       
                    
                    // calculate position (x,y) on UV grid as expected by fftw (no fft_unshift required)
                    int x_grid = 0, y_grid = 0;
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }
                    
                    
                 
//                 printf("DEBUG : u_index %d vs. %d ( u = %.2f , u_min = %.2f , delta_u = %.2f , u_center =%.2f)\n",u,u_index1,u_index,u_min,delta_u,u_center);
              
                 // Using CELL averaging method or setXY ?
                    int pos = y_grid*n_pixels + x_grid;
                    gridded_visibilities[pos][0] += re;
                    gridded_visibilities[pos][1] += im;
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );
                    
                    if( m_ImagerDebugLevel > 101 ){
                         printf("DEBUG ADDING : PACER_gridding at (x,y) =  %.1f %.1f += %.8f + j %.8f (u,v,w) = ( %.4f , %.4f , %.4f )\n",u_pix,v_pix,re,im,u,v,w);
//                       printf("DEBUG gridding : (u,v) = (%.8f,%.8f) -> index = (%d,%d) or (%.2f,%.2f), value = %.8f + j %.8f\n",u,v,u_index,v_index,u_pix,v_pix,re,im);
                    }
              
                 // add conjugates :
//                 u_index = round( (-u - u_min)/delta_u );
//                 v_index = round( (-v - v_min)/delta_v );
                    u_index = -u_pix + n_pixels/2; // was round( (-u - u_center)/delta_u ) + ...
                    v_index = -v_pix + n_pixels/2; // was round( (-v - v_center)/delta_v ) + ...
                    if( u_index < center_x ){
                       x_grid = center_x + u_index + is_odd_x;
                    }else{
                       x_grid = u_index - center_x;
                    }                    
                    if( v_index < center_y ){
                       y_grid = center_y + v_index + is_odd_y;
                    }else{
                       y_grid = v_index - center_y;
                    }                                       
                    
                    pos = y_grid*n_pixels + x_grid;
                    gridded_visibilities[pos][0] += re;
                    gridded_visibilities[pos][1] += -im;
                    uv_grid_counter.addXY( x_grid, y_grid , 1.00 );
                    
                    if( ant1==ant2 ){
                       PRINTF_DEBUG("Auto-correlation added values %.4f / %.4f\n",re,im);
                    }                    
                    
                    added++;
                 }
              }else{
                 PRINTF_DEBUG("DEBUG : visibility value %e +j%e higher than limit %e -> skipped\n",re,im,MAX_VIS);
                 high_value++;
              }
           }
        }
     }
  }  
  
  // This division is in fact UNIFORM weighting !!!! Not CELL-avareging 
  // normalisation to make it indeed CELL-averaging :
  // if( strcmp(weighting, "U" ) == 0 ){
  //   uv_grid_real.Divide( uv_grid_counter );
  //   uv_grid_imag.Divide( uv_grid_counter );
  // }  
  PACER_PROFILER_END("gridding (no I/O) took")
  PRINTF_DEBUG("DEBUG : added %d UV points to the grid, %d too high values skipped\n",added,high_value);
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
     char uv_grid_re_name[1024],uv_grid_im_name[1024],uv_grid_counter_name[1024];
     sprintf(uv_grid_re_name,"%s/uv_grid_real_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_im_name,"%s/uv_grid_imag_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
     sprintf(uv_grid_counter_name,"%s/uv_grid_counter_%dx%d.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),n_pixels,n_pixels);
    
/*     if( uv_grid_real.WriteFits( uv_grid_re_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_re_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_re_name);
     }

     if( uv_grid_imag.WriteFits( uv_grid_im_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_im_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_im_name);
     }*/
  
     if( uv_grid_counter.WriteFits( uv_grid_counter_name ) ){
        printf("ERROR : could not write output file %s\n",uv_grid_counter_name);
     }else{
        PRINTF_INFO("INFO : saved file %s\n",uv_grid_counter_name);
     }
  }


}                



void CPacerImager::gridding_imaging( Visibilities& xcorr, 
                                     int time_step, 
                                     int fine_channel,
                                     CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                                     double delta_u, double delta_v,
                                     double frequency_mhz,
                                     int    n_pixels,
                                     double min_uv /*=-1000*/,    // minimum UV 
                                     const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
                                     const char* szBaseOutFitsName /*=NULL*/, 
                                     bool do_gridding,
                                     bool do_dirty_image,
                                     const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
                                     const char* in_fits_file_uv_im, /*=""*/ // gridded visibilities can be provided externally
                                     bool bSaveIntermediate /*=false*/, bool bSaveImaginary /*=true*/
                )
{
  printf("DEBUG : gridding_imaging( Visibilities& xcorr ) in pacer_imager.cpp\n");
  // allocates data structures for gridded visibilities:
  AllocGriddedVis( n_pixels, n_pixels );
  
  if( do_gridding ){
     gridding_fast( xcorr, time_step, fine_channel, fits_vis_u, fits_vis_v, fits_vis_w, *m_uv_grid_real, *m_uv_grid_imag, *m_uv_grid_counter, delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting );
/*  }else{
     if( strlen(in_fits_file_uv_re) && strlen(in_fits_file_uv_im) ){
        uv_grid_counter.SetValue(1.00);
        
        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_re);
        if( uv_grid_real.ReadFits( in_fits_file_uv_re, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_re);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_re);
        }

        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_im);
        if( uv_grid_real.ReadFits( in_fits_file_uv_im, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_im);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_im);
        }

     }else{
        printf("ERROR : when gridding is disabled (-g 0) options -r and -i with REAL and IMAG FITS file names must be specified -> cannot continue !\n");
        exit(-1);
     }*/
  }

  if( do_dirty_image ){
     // dirty image :  
     PRINTF_INFO("PROGRESS : executing dirty image\n");
     dirty_image( *m_uv_grid_real, *m_uv_grid_imag, *m_uv_grid_counter, 
                  true, // do not save intermediate FITS files
                  szBaseOutFitsName, // output filename template
                  true,   // save imaginary image
                  false    // do FFTunshift true when gridding() , false for gridding_fast()
                );             
  }


}

void CPacerImager::gridding_imaging( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                                     double delta_u, double delta_v,
                                     double frequency_mhz,
                                     int    n_pixels,
                                     double min_uv /*=-1000*/,    // minimum UV 
                                     const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
                                     const char* szBaseOutFitsName /*=NULL*/, 
                                     bool do_gridding,
                                     bool do_dirty_image,
                                     const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
                                     const char* in_fits_file_uv_im, /*=""*/ // gridded visibilities can be provided externally
                                     bool bSaveIntermediate /*=false*/, bool bSaveImaginary /*=true*/
                )
{
  // TODO: make it member variable to optimise and do not allocate every time !!!
  CBgFits uv_grid_counter( n_pixels, n_pixels ),uv_grid_real( n_pixels, n_pixels ) , uv_grid_imag( n_pixels, n_pixels );  
  
  if( do_gridding ){
     gridding_fast( fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, uv_grid_real, uv_grid_imag, uv_grid_counter, delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting );
  }else{
     if( strlen(in_fits_file_uv_re) && strlen(in_fits_file_uv_im) ){
        uv_grid_counter.SetValue(1.00);
        
        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_re);
        if( uv_grid_real.ReadFits( in_fits_file_uv_re, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_re);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_re);
        }

        PRINTF_INFO("Reading fits file %s ...\n",in_fits_file_uv_im);
        if( uv_grid_imag.ReadFits( in_fits_file_uv_im, 0, 1, 1 ) ){
           printf("ERROR : could not read visibility FITS file %s\n",in_fits_file_uv_im);
           exit(-1); 
        }else{
           PRINTF_INFO("OK : fits file %s read ok\n",in_fits_file_uv_im);
        }

     }else{
        printf("ERROR : when gridding is disabled (-g 0) options -r and -i with REAL and IMAG FITS file names must be specified -> cannot continue !\n");
        exit(-1);
     }
  }

  if( do_dirty_image ){
     // dirty image :  
     PRINTF_INFO("PROGRESS : executing dirty image\n");
     dirty_image( uv_grid_real, uv_grid_imag, uv_grid_counter, 
                  true, // do not save intermediate FITS files
                  szBaseOutFitsName, // output filename template
                  true,   // save imaginary image
                  false    // do FFTunshift true when gridding() , false for gridding_fast()
                );             
  }


}                

bool CPacerImager::ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz, int time_step, int fine_channel )
{
   int n_ant = xcorr.obsInfo.nAntennas;
   double frequency_hz = frequency_mhz*1e6;
   printf("DEBUG (ApplyGeometricCorrections(xcorr)) : freq = %.8f [Hz]\n",frequency_hz);
   
   CBgFits fits_vis_real(n_ant,n_ant);
   CBgFits fits_vis_imag(n_ant,n_ant);   
   
   for(int ant1=0;ant1<n_ant;ant1++){
      for(int ant2=0;ant2<=ant1;ant2++){ // loop over ant2
      // for(int ant2=0;ant2<n_ant;ant2++){
         double w = fits_vis_w.getXY(ant1,ant2); // or ant2,ant1 , was (+1)
         double angle = +2.0*M_PI*w*frequency_hz / SPEED_OF_LIGHT; // TODO : + or - here ??? In brute force branch was -
         
         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, ant1, ant2 );
         
         double re = vis[0].real();
         double im = vis[0].imag();
         // double re = fits_vis_real.getXY(ant1,ant2);
         // double im = fits_vis_imag.getXY(ant1,ant2);
         
         if( ant1 == 1 && ant2 == 0 ){
            printf("CPU : (%d,%d) , angle = %.8f , w = %.8f\n",ant1,ant2,angle,w);
         }

         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;
         
         // change in xcorr :
         std::complex<double> vis_new ( re_prim, im_prim );
         (*vis) = vis_new;
         
         fits_vis_real.setXY(ant1,ant2,re_prim);
         fits_vis_imag.setXY(ant1,ant2,-im_prim);
         fits_vis_real.setXY(ant2,ant1,re_prim);
         fits_vis_imag.setXY(ant2,ant1,im_prim);
         
         printf("\t%d-%d : w = %.4f [m] -> angle = %.8f (sin = %.8f , cos = %.8f) * re=%.8f im=%.8f -> %.4f + j%.4f, freq = %.8f Hz\n",ant1,ant2,w,angle,sin_angle,cos_angle,re,im,re_prim,im_prim,frequency_hz);
      }
   }
   
   char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
   sprintf(szOutPutFitsRE,"%s/vis_re_geom_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   sprintf(szOutPutFitsIM,"%s/vis_im_geom_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());   
   fits_vis_real.WriteFits( szOutPutFitsRE );
   fits_vis_imag.WriteFits( szOutPutFitsIM );

   return true;
}

bool CPacerImager::ApplyGeometricCorrections( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz )
{
   int n_ant = fits_vis_real.GetXSize();
   double frequency_hz = frequency_mhz*1e6;
   printf("DEBUG (ApplyGeometricCorrections) : freq = %.8f [Hz]\n",frequency_hz);
   
   for(int ant1=0;ant1<n_ant;ant1++){
      for(int ant2=0;ant2<n_ant;ant2++){
         double w = fits_vis_w.getXY(ant1,ant2); // or ant2,ant1 , was (+1)
         double angle = -2.0*M_PI*w*frequency_hz / SPEED_OF_LIGHT;
         
         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         double re = fits_vis_real.getXY(ant1,ant2);
         double im = fits_vis_imag.getXY(ant1,ant2);
         
         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;
         fits_vis_real.setXY(ant1,ant2,re_prim);
         fits_vis_imag.setXY(ant1,ant2,im_prim);
         
         if( ant1 == 1 && ant2 == 0 ){
            printf("CPU : (%d,%d) , angle = %.8f , w = %.8f\n",ant1,ant2,angle,w);
         }
         
         printf("\t%d-%d : w = %.4f [m] -> angle = %.8f (sin = %.8f , cos = %.8f) * re=%.8f im=%.8f -> %.4f + j%.4f, freq = %.8f Hz\n",ant1,ant2,w,angle,sin_angle,cos_angle,re,im,re_prim,im_prim,frequency_hz);
      }
   }

   char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
   sprintf(szOutPutFitsRE,"%s/vis_re_geom_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   sprintf(szOutPutFitsIM,"%s/vis_im_geom_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());   
   fits_vis_real.WriteFits( szOutPutFitsRE );
   fits_vis_imag.WriteFits( szOutPutFitsIM );
   
   return true;
}

bool CPacerImager::ApplyCableCorrections( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz )
{
   int n_ant = fits_vis_real.GetXSize();
   double frequency_hz = frequency_mhz*1e6;
   printf("DEBUG (ApplyCableCorrections) : freq = %.8f [Hz]\n",frequency_hz);
   
   for(int ant1=0;ant1<n_ant;ant1++){
      InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant1];
            
      for(int ant2=0;ant2<n_ant;ant2++){
         InputMapping& ant2_info = m_MetaData.m_AntennaPositions[ant2];
         
         double cableDeltaLen = ant1_info.cableLenDelta - ant2_info.cableLenDelta; // was ant2 - ant1 
         double angle = -2.0 * M_PI * cableDeltaLen * frequency_hz / SPEED_OF_LIGHT;
                  
         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         double re = fits_vis_real.getXY(ant1,ant2);
         double im = fits_vis_imag.getXY(ant1,ant2);
         
         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;
         
         if( ant2 == 0 ){
            printf("\t%d-%d : cable = %.4f [m] -> angle = %.8f (sin = %.8f , cos = %.8f) * re=%.8f im=%.8f -> corrected re_prim/im_prim = %.8f/%.8f\n",ant1,ant2,cableDeltaLen,angle,sin_angle,cos_angle,re,im,re_prim,im_prim);
         }

         
         fits_vis_real.setXY(ant1,ant2,re_prim);
         fits_vis_imag.setXY(ant1,ant2,im_prim);         
      }
   }
   
   char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
   sprintf(szOutPutFitsRE,"%s/vis_re_cable_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   sprintf(szOutPutFitsIM,"%s/vis_im_cable_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   fits_vis_real.WriteFits( szOutPutFitsRE );
   fits_vis_imag.WriteFits( szOutPutFitsIM );

   return true;
}

bool CPacerImager::ApplyCableCorrections( Visibilities& xcorr, double frequency_mhz, int time_step, int fine_channel )
{
   int n_ant = xcorr.obsInfo.nAntennas;

   double frequency_hz = frequency_mhz*1e6;
   printf("DEBUG (ApplyCableCorrections) : freq = %.8f [Hz]\n",frequency_hz);
   
   CBgFits fits_vis_real(n_ant,n_ant);
   CBgFits fits_vis_imag(n_ant,n_ant);   

   for(int ant1=0;ant1<n_ant;ant1++){
      InputMapping& ant1_info = m_MetaData.m_AntennaPositions[ant1];

      for(int ant2=0;ant2<=ant1;ant2++){ // loop over ant2 
      // for(int ant2=0;ant2<n_ant;ant2++){
         InputMapping& ant2_info = m_MetaData.m_AntennaPositions[ant2];

         double cableDeltaLen = ant2_info.cableLenDelta - ant1_info.cableLenDelta; // was ant2 - ant1 
         double angle = -2.0 * M_PI * cableDeltaLen * frequency_hz / SPEED_OF_LIGHT;

         double sin_angle,cos_angle;
         sincos(angle, &sin_angle, &cos_angle);
         
         std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, ant1, ant2 );
  
         double re = vis[0].real();
         double im = vis[0].imag(); // TODO : why do I need 1 here ???

//         double re = fits_vis_real.getXY(ant1,ant2);
//         double im = fits_vis_imag.getXY(ant1,ant2);
         
         double re_prim = re*cos_angle - im*sin_angle;
         double im_prim = im*cos_angle + re*sin_angle;
         
         if( ant2 == 0 ){
            printf("\t%d-%d : cable = %.4f [m] -> angle = %.8f (sin = %.8f , cos = %.8f) * re=%.8f im=%.8f -> corrected re/im = %.8f / %.8f for %.8f [Hz]\n",ant1,ant2,cableDeltaLen,angle,sin_angle,cos_angle,re,im,re_prim,im_prim,frequency_hz);
         }

         
         // change in xcorr :
         std::complex<double> vis_new ( re_prim, im_prim );
         (*vis) = vis_new;
         
//         double re_prim = re;
//         double im_prim = im;
         fits_vis_real.setXY(ant1,ant2,re_prim);
         fits_vis_imag.setXY(ant1,ant2,-im_prim);         
         fits_vis_real.setXY(ant2,ant1,re_prim);
         fits_vis_imag.setXY(ant2,ant1,im_prim);         
      }
   }
   
   char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
   sprintf(szOutPutFitsRE,"%s/vis_re_cable_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   sprintf(szOutPutFitsIM,"%s/vis_im_cable_corr.fits",m_ImagerParameters.m_szOutputDirectory.c_str());   
   fits_vis_real.WriteFits( szOutPutFitsRE );
   fits_vis_imag.WriteFits( szOutPutFitsIM );
   
   return true;
}

bool CPacerImager::ApplySolutions(  Visibilities& xcorr, double frequency_mhz, CCalSols& calsol, int time_step, int fine_channel, const char* szPol )
{
   printf("DEBUG : CPacerImager::ApplySolutions( xcorr - 20230928 ) frequency = %.4f [MHz]\n",frequency_mhz);
   double freq_diff = fabs( frequency_mhz - calsol.m_frequency_mhz );
   
   if( freq_diff > 5.00 ){
      printf("WARNING : frequency of calibration solutions is different by more than 1 MHz from the actual data -> no calibration applied (data frequency %.1f vs. calsol frequency %.1f)\n",frequency_mhz,calsol.m_frequency_mhz);
      
      return false;
   }
   
   if( xcorr.obsInfo.nAntennas != calsol.size() ){
      printf("WARNING : wrong number of calibration solutions (%d vs. required %d)\n",int(calsol.size()),xcorr.obsInfo.nAntennas);
      
      return false;
   }
   
   int pol_idx = 0;
   if ( strcmp(szPol,"Y") == 0 ){
      pol_idx = 3;
   }
   
   int n_ant = xcorr.obsInfo.nAntennas;
   CBgFits fits_vis_real(n_ant,n_ant);
   CBgFits fits_vis_imag(n_ant,n_ant);  
   
   for(int i=0;i<n_ant;i++){ // loop over ant1          
      for(int j=0;j<=i;j++){ // loop over ant2 
         std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, i, j );
         std::complex<double> z( vis->real(), vis->imag() );

//         double re = fits_vis_real.getXY(x,y);
//         double im = fits_vis_imag.getXY(x,y);         
//         std::complex<double> z(re,im);
         
         std::complex<double> g1 = calsol[j].m_cal[pol_idx]; // or y ? TODO , TBD 
         std::complex<double> g2_c = std::conj( calsol[i].m_cal[pol_idx] ); // or x ? TODO , TBD
         
// 1/g :         
//         std::complex<double> z_cal = (1.00/g1)*z*(1.00/g2_c); // or g1 , g2_c without 1.00/ ??? TBD vs. --invert amplitude in read cal. sol. from MCCS DB 
// g :
         std::complex<double> z_cal = (g1)*z*(g2_c);
         if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
            if( i<=5 && j<=5 ){
               printf("DEBUG : z = %.4f + i%4f -> %.4f + i%4f\n",z.real(),z.imag(),z_cal.real(),z_cal.imag());
            }
         }
         
         (*vis) = z_cal;
         fits_vis_real.setXY(i,j,z_cal.real());
         fits_vis_imag.setXY(i,j,-z_cal.imag());
         fits_vis_real.setXY(j,i,z_cal.real());
         fits_vis_imag.setXY(j,i,z_cal.imag());
         
//         fits_vis_real.setXY(x,y,z_cal.real());
//         fits_vis_imag.setXY(x,y,z_cal.imag());
      }
   }

   // may be temporary :
   char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
   sprintf(szOutPutFitsRE,"%s/xcorr_re.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   fits_vis_real.WriteFits( szOutPutFitsRE );
   sprintf(szOutPutFitsIM,"%s/xcorr_im.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
   fits_vis_imag.WriteFits( szOutPutFitsIM );
   printf("DEBUG : saved %s and %s\n",szOutPutFitsRE,szOutPutFitsIM);
   
   return true;
   
}

bool CPacerImager::ApplySolutions( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz, const char* szPol /*="X"*/ )
{
   return ApplySolutions( fits_vis_real, fits_vis_imag, frequency_mhz, m_CalibrationSolutions, szPol );
}

bool CPacerImager::ApplySolutions( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz, CCalSols& calsol, const char* szPol )
{
   printf("DEBUG : CPacerImager::ApplySolutions frequency = %.4f [MHz]\n",frequency_mhz);
   double freq_diff = fabs( frequency_mhz - calsol.m_frequency_mhz );
   
   if( freq_diff > 5.00 ){
      printf("WARNING : frequency of calibration solutions is different by more than 1 MHz from the actual data -> no calibration applied (data frequency %.1f vs. calsol frequency %.1f)\n",frequency_mhz,calsol.m_frequency_mhz);
      
      return false;
   }
   
   if( fits_vis_real.GetXSize() != calsol.size() ){
      printf("WARNING : wrong number of calibration solutions (%d vs. required %ld)\n",int(calsol.size()),fits_vis_real.GetXSize());
      
      return false;
   }
   
   int pol_idx = 0;
   if ( strcmp(szPol,"Y") == 0 ){
      pol_idx = 3;
   }
   
   int n_ant = fits_vis_real.GetXSize();
   for(int y=0;y<n_ant;y++){
      for(int x=0;x<n_ant;x++){
         double re = fits_vis_real.getXY(x,y);
         double im = fits_vis_imag.getXY(x,y);
         
         std::complex<double> z(re,im);
         
         std::complex<double> g1 = calsol[x].m_cal[pol_idx]; // or y ? TODO , TBD 
         std::complex<double> g2_c = std::conj( calsol[y].m_cal[pol_idx] ); // or x ? TODO , TBD
         
// 1/g :         
//         std::complex<double> z_cal = (1.00/g1)*z*(1.00/g2_c); // or g1 , g2_c without 1.00/ ??? TBD vs. --invert amplitude in read cal. sol. from MCCS DB 
// g :
         std::complex<double> z_cal = (g1)*z*(g2_c);
         if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){
            if( x<=5 && y<=5 ){
               printf("DEBUG : z = %.4f + i%4f -> %.4f + i%4f\n",z.real(),z.imag(),z_cal.real(),z_cal.imag());
            }
         }
         
         fits_vis_real.setXY(x,y,z_cal.real());
         fits_vis_imag.setXY(x,y,z_cal.imag());
      }
   }
   
   // may be temporary :
   if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG ){
      char szOutPutFitsRE[1024],szOutPutFitsIM[1024];
      sprintf(szOutPutFitsRE,"%s/cal_re.fits",m_ImagerParameters.m_szOutputDirectory.c_str());
      sprintf(szOutPutFitsIM,"%s/cal_im.fits",m_ImagerParameters.m_szOutputDirectory.c_str());   
      fits_vis_real.WriteFits( szOutPutFitsRE );
      fits_vis_imag.WriteFits( szOutPutFitsIM );
   }

   
   return true;
}


bool CPacerImager::run_imager( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                               double frequency_mhz, 
                               int n_pixels,
                               double FOV_degrees,
                               double min_uv,                  /*=-1000,*/
                               bool do_gridding,               /*=true*/
                               bool do_dirty_image,            /*=true*/
                               const char* weighting,          /* ="" */   // weighting : U for uniform (others not implemented)
                               const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
                               const char* in_fits_file_uv_im, /*=""*/ // gridded visibilities can be provided externally
                               const char* szBaseOutFitsName   /*=NULL*/
)
{
  // ensures initalisation of object structures 
  Initialise();
  PRINTF_DEBUG("Number of flagged antennas = %d (CPacerImager::run_imager)\n",int(m_FlaggedAntennas.size()));
  
  if( m_ImagerParameters.m_fUnixTime <= 0.0001 ){
     m_ImagerParameters.m_fUnixTime = get_dttm_decimal();
     PRINTF_WARNING("Time of the data not specified -> setting current time %.6f\n",m_ImagerParameters.m_fUnixTime);
  }

  PACER_PROFILER_START
  
  // based on RTS : UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c  
  double frequency_hz = frequency_mhz*1e6;
  double wavelength_m = VEL_LIGHT / frequency_hz;
  double FoV_radians = FOV_degrees*M_PI/180.;
  // WARNING: it actually cancels out to end up being 1/PI :
  // TODO simplity + 1/(2FoV) !!! should be 
//  double delta_u = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in meters (NOT WAVELENGHTS)
//  double delta_v = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in meters (NOT WAVELENGHTS)
// TODO : 
  double delta_u = 1.00/(FoV_radians); // should be 2*FoV_radians - see TMS etc 
  double delta_v = 1.00/(FoV_radians); // Rick Perley page 16 : /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
  
  if( true ){
     // see : Rick Perley page 16 : /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
     //       and /home/msok/Desktop/PAWSEY/PaCER/logbook/20240315_Jishnu_test.odt
     
     // WARNING : why this is not u_max/wavelength_m - when I use this the image because Frequency dependent. But it does not make sense to have 
     //           delta_u in meters and then UV grid in wavelengths ...
     // MWA TEST:
     delta_u = 2.00*(u_max)/n_pixels; // Looks like this is what it should be NOT u_max/wavelength_m . So delta_u must be in meters here !!! It may all depend on calculation if u_index 
                                      // see discussion in /home/msok/Desktop/PAWSEY/PaCER/logbook/20240320_gridding_delta_u_meters_vs_wavelengths.odt     
     delta_v = 2.00*(v_max)/n_pixels; // and it's not ok because then delta_u is different for both of them, which causes exp/shrink with freq     
     
     // automatic calculation of pixel size in radians 1/(2u_max) - see Rick Perley or just Lambda/B_max divide by 2 for oversampling.
     m_ImagerParameters.m_PixsizeInRadians = 1.00/(2.00*u_max); // does this one need to be /wavelength or not ???
     
     // NEW : based on desired image resolution 
     // double delta_theta = (wavelength_m/35.0)/2.00; // 2 times oversampled
     // double delta_theta = ((230.0/300.00)/(2.00*35.00)); // use maximum resolution oversampled by a factor of 2., at 230 MHz -> Lambda ~1.3043m
     double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
// MWA TEST 
//     delta_u = 1.00/(n_pixels*delta_theta);
//     delta_v = 1.00/(n_pixels*delta_theta);
     PRINTF_DEBUG("delta_u = %.8f (u_max = %.8f), delta_v = %.8f (v_max = %.8f), calculated as 1/FoV = 1/(%d pixels * %.5f rad), delta_theta = %.5f [deg]\n",delta_u,u_max,delta_v,v_max,n_pixels,delta_theta,delta_theta*(180.00/M_PI));
     
     // PRINTF_DEBUG("delta_u = %.8f , delta_v = %.8f , calculated as 2.00*u_max/n_pixels, u_max = %.8f, n_pixels = %d\n",delta_u,delta_v,u_max,n_pixels);
  }
    
  // PRINTF_DEBUG("delta_u = %.8f , delta_v = %.8f , calculated as 1/FoV_radians = 1/%.8f (NOT THIS WAY)\n",delta_u,delta_v,FoV_radians);
  
  // test forced to this value (based on pixscale in MIRIAD):
  if( m_bCompareToMiriad ){
     // Brute force comparison to MIRIAD assuming pixscale from the final image FITS FILE = 0.771290200761 degree 
     if( fabs(frequency_mhz-159.375) <= 0.5 ){
        delta_u = 0.412697967; // at ch=204
        delta_v = 0.412697967; // at ch=204
     }
     if( fabs(frequency_mhz-246.09375) <= 0.5 ){
        delta_u = .63624180895350226815; // at ch=315
        delta_v = .63624180895350226815; // at ch=315
     }
  }
  
  // Apply Calibration if provided :
  if( m_CalibrationSolutions.size() > 0 )
  {
     if( m_CalibrationSolutions.size() == fits_vis_real.GetXSize() )
     {
        printf("INFO : applying calibration solutions (for %d antennas)\n",int(m_CalibrationSolutions.size()));
        ApplySolutions( fits_vis_real, fits_vis_imag, frequency_mhz, m_CalibrationSolutions );
     }else
     {
        printf("WARNING : wrong number of calibration solutions (%d vs. required %ld)\n",int(m_CalibrationSolutions.size()),fits_vis_real.GetXSize());
     }
  }else{
     printf("WARNING : no calibration solutions (size=%d) -> no applied!\n",int(m_CalibrationSolutions.size()));
  }

  // 2024-04-16 - hopefully not applied double time !
  if( m_ImagerParameters.m_bApplyGeomCorr ){
     printf("CPacerImager::run_imager(FITS) - ApplyGeometricCorrections call\n");
     ApplyGeometricCorrections( fits_vis_real, fits_vis_imag, m_U, m_V, m_W, frequency_mhz );
  }

  if( m_ImagerParameters.m_bApplyCableCorr ){
     printf("CPacerImager::run_imager(FITS) - ApplyCableCorrections call\n");
     ApplyCableCorrections( fits_vis_real, fits_vis_imag, frequency_mhz );
  }

  if( do_gridding || do_dirty_image ){
     // virtual function calls gridding and imaging in GPU/HIP version it is overwritten and 
     // both gridding and imaging are performed on GPU memory :
     gridding_imaging( fits_vis_real, fits_vis_imag, fits_vis_u, fits_vis_v, fits_vis_w, 
                       delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting,
                       szBaseOutFitsName, do_gridding, do_dirty_image, in_fits_file_uv_re, in_fits_file_uv_im 
                     );
  }

  // increse image counter:
  m_SkyImageCounter++;
  if( m_SkyImageCounter >= INT_MAX )
  {
    PRINTF_WARNING("WARNING : image counter reached maximum value for int = %d -> reset to zero\n",INT_MAX);
    m_SkyImageCounter = 0;
  }
  
  PACER_PROFILER_END("full imaging (gridding + dirty image) took")

  return  true;   
}


//-----------------------------------------------------------------------------------------------------------------------------
// Wrapper to run_imager as before : run_imager( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w ...)
// it uses AstroIO class Visibility 
// Reads FITS files and executes overloaded function run_imager ( as above )
//-----------------------------------------------------------------------------------------------------------------------------
bool CPacerImager::run_imager( Visibilities& xcorr, 
                 int time_step, 
                 int fine_channel,
                 CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,  // UVW
                 double frequency_mhz,
                 int    n_pixels,
                 double FOV_degrees,
                 double min_uv,                      // =-1000, minimum UV 
                 bool   do_gridding,                 // =true, excute gridding  (?)
                 bool   do_dirty_image,              // =true, form dirty image (?)
                 const char* weighting,              // ="",  weighting : U for uniform (others not implemented)
                 const char* in_fits_file_uv_re,     // ="", gridded visibilities can be provided externally
                 const char* in_fits_file_uv_im,     // ="",  gridded visibilities can be provided externally
                 const char* szBaseOutFitsName       // =NULL
                  )
{
  // ensures initalisation of object structures 
  Initialise();
  
  if( m_ImagerParameters.m_fUnixTime <= 0.0001 ){
     m_ImagerParameters.m_fUnixTime = get_dttm_decimal();
     PRINTF_WARNING("Time of the data not specified -> setting current time %.6f\n",m_ImagerParameters.m_fUnixTime);
  }

  PACER_PROFILER_START
  
  // based on RTS : UV pixel size as function FOVtoGridsize in  /home/msok/mwa_software/RTS_128t/src/gridder.c  
  double frequency_hz = frequency_mhz*1e6;
  double wavelength_m = VEL_LIGHT / frequency_hz;
  double FoV_radians = FOV_degrees*M_PI/180.;
  // WARNING: it actually cancels out to end up being 1/PI :
  // TODO simplity + 1/(2FoV) !!! should be 
//  double delta_u = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in meters (NOT WAVELENGHTS)
//  double delta_v = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in meters (NOT WAVELENGHTS)
// TODO : 
  double delta_u = 1.00/(FoV_radians); // should be 2*FoV_radians - see TMS etc 
  double delta_v = 1.00/(FoV_radians); // Rick Perley page 16 : /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
  
  printf("DEBUG : frequency = %.6f [MHz] -> %.2f [Hz]\n",frequency_mhz,frequency_hz);

  if( true ){
     // see : Rick Perley page 16 : /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
     //       and /home/msok/Desktop/PAWSEY/PaCER/logbook/20240315_Jishnu_test.odt
     
     // WARNING : why this is not u_max/wavelength_m - when I use this the image because Frequency dependent. But it does not make sense to have 
     //           delta_u in meters and then UV grid in wavelengths ...
     // MWA TEST:
     delta_u = 2.00*(u_max)/n_pixels; // delta_u = 2.00*(u_max)/n_pixels; // Looks like this is what it should be NOT u_max/wavelength_m . So delta_u must be in meters here !!! It may all depend on calculation if u_index 
                                      // see discussion in /home/msok/Desktop/PAWSEY/PaCER/logbook/20240320_gridding_delta_u_meters_vs_wavelengths.odt     
     delta_v = 2.00*(v_max)/n_pixels; // and it's not ok because then delta_u is different for both of them, which causes exp/shrink with freq     
     
     // automatic calculation of pixel size in radians 1/(2u_max) - see Rick Perley or just Lambda/B_max divide by 2 for oversampling.
     m_ImagerParameters.m_PixsizeInRadians = 1.00/(2.00*u_max); // does this one need to be /wavelength or not ???
     
     // NEW : based on desired image resolution 
     // double delta_theta = (wavelength_m/35.0)/2.00; // 2 times oversampled
     // double delta_theta = ((230.0/300.00)/(2.00*35.00)); // use maximum resolution oversampled by a factor of 2., at 230 MHz -> Lambda ~1.3043m
     // double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
     double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
// MWA TEST 
//     delta_u = 1.00/(n_pixels*delta_theta);
//     delta_v = 1.00/(n_pixels*delta_theta);
     PRINTF_DEBUG("delta_u = %.8f (u_max = %.8f), delta_v = %.8f (v_max = %.8f), calculated as 1/FoV = 1/(%d pixels * %.5f rad), delta_theta = %.5f [deg]\n",delta_u,u_max,delta_v,v_max,n_pixels,delta_theta,delta_theta*(180.00/M_PI));
     
     // PRINTF_DEBUG("delta_u = %.8f , delta_v = %.8f , calculated as 2.00*u_max/n_pixels, u_max = %.8f, n_pixels = %d\n",delta_u,delta_v,u_max,n_pixels);
  }
  
  // test forced to this value (based on pixscale in MIRIAD):
  if( m_bCompareToMiriad ){
     // Brute force comparison to MIRIAD assuming pixscale from the final image FITS FILE = 0.771290200761 degree 
     if( fabs(frequency_mhz-159.375) <= 0.5 ){
        delta_u = 0.412697967; // at ch=204
        delta_v = 0.412697967; // at ch=204
     }
     if( fabs(frequency_mhz-246.09375) <= 0.5 ){
        delta_u = .63624180895350226815; // at ch=315
        delta_v = .63624180895350226815; // at ch=315
     }
  }
  
  // Apply Calibration if provided :
  // TODO: Calibration needs to be implemented for AstroIO Visibility library, currently it only takes CBgFits class as an input:
  if( m_CalibrationSolutions.size() > 0 )
  {
     if( m_CalibrationSolutions.size() == xcorr.obsInfo.nAntennas )
     {
        printf("INFO : applying calibration solutions (for %d antennas)\n",int(m_CalibrationSolutions.size()));
        // ApplySolutions( fits_vis_real, fits_vis_imag, frequency_mhz, m_CalibrationSolutions );
        ApplySolutions( xcorr, frequency_mhz, m_CalibrationSolutions, time_step, fine_channel );
        
     }else
     {
        printf("WARNING : wrong number of calibration solutions (%d vs. required %d)\n",int(m_CalibrationSolutions.size()),xcorr.obsInfo.nAntennas);
     }
//     printf("ERROR : application of calibration is not implemented in the current version\n");
//     exit(-1);
  }

  if( do_gridding || do_dirty_image ){
     // virtual function calls gridding and imaging in GPU/HIP version it is overwritten and 
     // both gridding and imaging are performed on GPU memory :
// TODO : create gridding_imaging( Visibilities& xcorr, ... )
     gridding_imaging( xcorr, time_step, fine_channel,
                       fits_vis_u, fits_vis_v, fits_vis_w, 
                       delta_u, delta_v, frequency_mhz, n_pixels, min_uv, weighting,
                       szBaseOutFitsName, do_gridding, do_dirty_image, in_fits_file_uv_re, in_fits_file_uv_im 
                     );
  }


  // increse image counter:
  m_SkyImageCounter++;
  if( m_SkyImageCounter >= INT_MAX )
  {
    PRINTF_WARNING("WARNING : image counter reached maximum value for int = %d -> reset to zero\n",INT_MAX);
    m_SkyImageCounter = 0;
  }
  
  PACER_PROFILER_END("full imaging (gridding + dirty image) took")

  return  true;   
   
}                  


//-----------------------------------------------------------------------------------------------------------------------------
// Wrapper to run_imager :
// Reads FITS files and executes overloaded function run_imager ( as above )
//-----------------------------------------------------------------------------------------------------------------------------
bool CPacerImager::run_imager( const char* basename, const char* szPostfix,
                               bool   do_gridding             /*=true*/, // excute gridding  (?)
                               bool   do_dirty_image          /*=true*/, // form dirty image (?)
                               const char* in_fits_file_uv_re /*=""*/,   // gridded visibilities can be provided externally
                               const char* in_fits_file_uv_im /*=""*/,   // gridded visibilities can be provided externally                    
                               const char* szBaseOutFitsName  /*=NULL*/
                             )
{
    return run_imager( basename, szPostfix, 
                       m_ImagerParameters.m_fCenterFrequencyMHz,
                       m_ImagerParameters.m_ImageSize,
                       m_ImagerParameters.m_ImageFOV_degrees,
                       m_ImagerParameters.m_fMinUV,
                       do_gridding,
                       do_dirty_image,
                       m_ImagerParameters.m_szWeighting.c_str(),
                       in_fits_file_uv_re, in_fits_file_uv_im,
                       szBaseOutFitsName
                     );
}                             
                             
                             
bool CPacerImager::run_imager( const char* basename, const char* szPostfix,
                               double frequency_mhz, 
                               int n_pixels,
                               double FOV_degrees,
                               double min_uv,                  /*=-1000,*/
                               bool do_gridding,               /*=true*/
                               bool do_dirty_image,            /*=true*/
                               const char* weighting,          /* ="" */   // weighting : U for uniform (others not implemented)
                               const char* in_fits_file_uv_re, /*=""*/ // gridded visibilities can be provided externally
                               const char* in_fits_file_uv_im, /*=""*/ // gridded visibilities can be provided externally                               
                               const char* szBaseOutFitsName   /*=NULL*/
                  )
{
   // ensures initalisation of object structures 
   Initialise();  

   // read input data (correlation matrix and UVW) :
   CBgFits fits_vis_real, fits_vis_imag;
   if( read_corr_matrix( basename, fits_vis_real, fits_vis_imag, szPostfix ) )
   { // also included reading or calculation of UVW 
      PRINTF_INFO("OK : input files read ok\n");
   }else
   {
      printf("ERROR : could not read one of the input files\n");
      return false;
   }
   
   // read calibration solutions (if specified) :
   if( m_CalibrationSolutions.read_calsolutions() > 0 )
   {
      m_CalibrationSolutions.show();
   }

   // TODO : once application of calibration is implemented in the xcorr-path of the code -> comment/remove the below call and uncomment the code later:
   bool ret = run_imager( fits_vis_real, fits_vis_imag, m_U, m_V, m_W,
                         frequency_mhz, 
                         n_pixels, 
                         FOV_degrees, 
                         min_uv,
                         do_gridding, 
                         do_dirty_image, 
                         weighting, 
                         in_fits_file_uv_re, in_fits_file_uv_im, 
                         szBaseOutFitsName
                        );

   // READY TO GO CODE FOR THE FUTURE - please do not remove - see comment above
   // WARNING : this version already works ok, but does not pass test 1 because application of calibration is not implemented yet !
   // TODO : implemented application of calibration for xcorr-version !
   // TODO : New version, which requires allocation of xcorr.data (some new function in astroio Visibilities):
   // convert from CBgFits re/im to xcorr to use the same function:

/*   printf("DEBUG : this is experimental version of code CBgFits -> xcorr -> run_imager( xcorr )\n");fflush(stdout);
   Visibilities* xcorr = ConvertFits2XCorr( fits_vis_real, fits_vis_imag );
   printf("DEBUG : this is experimental version after ConvertFits2XCorr\n");fflush(stdout);
   // 0,0 -> because there is only single timestamp and frequency for CBgFits input :
   bool ret = run_imager(  *xcorr, 0, 0, frequency_mhz, n_pixels, FOV_degrees, min_uv, do_gridding, do_dirty_image, weighting, szBaseOutFitsName );

   // remove xcorr which was allocated (new) inside ConvertFits2XCorr   
   delete xcorr;
*/

   return ret;               
}

//-----------------------------------------------------------------------------------------------------------------------------
// Wrapper to run_imager - executes overloaded function run_imager ( as above ) :
//  INPUT : pointer to data
//-----------------------------------------------------------------------------------------------------------------------------
bool CPacerImager::run_imager( float* data_real, 
                               float* data_imag,
                               int n_ant, 
                               int n_pol,
                               double frequency_mhz, 
                               int n_pixels,
                               double FOV_degrees,
                               double min_uv,                   /*=-1000,*/
                               bool do_gridding,                /*=true  */
                               bool do_dirty_image,             /*=true  */
                               const char* weighting,           /* =""   */   // weighting : U for uniform (others not implement
                               const char* szBaseOutFitsName,   /* =NULL */
                               bool bCornerTurnRequired         /* =true , TODO : change default to false and perform corner-turn in eda2 imager code using correlation matrix from xGPU correlator */
                             )
{
  // ensures initalisation of object structures 
  Initialise();

  CBgFits fits_vis_real( n_ant, n_ant ), fits_vis_imag( n_ant, n_ant );
  
  if( bCornerTurnRequired )
  {
     fits_vis_real.SetValue(0.00);
     fits_vis_imag.SetValue(0.00);

     // ~corner-turn operation which can be quickly done on GPU: 
     int counter=0;
     for(int i=0;i<n_ant;i++)
     {
        for(int j=0;j<(i + 1);j++)
        {
           fits_vis_real.setXY( i, j , data_real[counter] );
           fits_vis_real.setXY( j, i , data_real[counter] );
           counter++;
        }
     }

     counter=0;
     for(int i=0;i<n_ant;i++)
     {
        for(int j=0;j<(i + 1);j++)
        {
           fits_vis_imag.setXY( i, j , data_imag[counter] );
           fits_vis_imag.setXY( j, i , -(data_imag[counter]) );
           counter++;
        }
     }
  }else
  {
     fits_vis_real.SetData( data_real );
     fits_vis_imag.SetData( data_imag );
  }
  
  // const char* get_filename(  time_t ut_time , char* out_buffer, int usec=0, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix=""
  if( m_ImagerParameters.m_fUnixTime <= 0.0001 ){
     m_ImagerParameters.m_fUnixTime = get_dttm_decimal();
     PRINTF_WARNING("Time of the data not specified -> setting current time %.6f\n",m_ImagerParameters.m_fUnixTime);
  }
  
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_INFO ){
     char outDirtyImageReal[1024];
     // sprintf(outDirtyImageReal,"%s/%s_vis_real.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
     get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageReal, m_ImagerParameters.m_szOutputDirectory.c_str(), "visibility_", "_real" );
   
     fits_vis_real.WriteFits( outDirtyImageReal );
     PRINTF_DEBUG("Saved real file to %s\n",outDirtyImageReal);
  }
     
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_INFO ){
     char outDirtyImageImag[1024];
     // sprintf(outDirtyImageImag,"%s/%s_vis_imag.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
     get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageImag, m_ImagerParameters.m_szOutputDirectory.c_str(), "visibility_", "_imag" );
        
     fits_vis_imag.WriteFits( outDirtyImageImag );
     PRINTF_DEBUG("Saved imag file to %s\n",outDirtyImageImag);
  }

  // calculate UVW (if required)
  CalculateUVW();
  
  if( m_ImagerParameters.m_bApplyGeomCorr ){
     ApplyGeometricCorrections( fits_vis_real, fits_vis_imag, m_U, m_V, m_W, frequency_mhz );
  }
  
  if( m_ImagerParameters.m_bApplyCableCorr ){
     ApplyCableCorrections( fits_vis_real, fits_vis_imag, frequency_mhz );
  }
  
  bool ret = run_imager( fits_vis_real, fits_vis_imag, m_U, m_V, m_W,
                         frequency_mhz, 
                         n_pixels, 
                         FOV_degrees, 
                         min_uv,
                         do_gridding, 
                         do_dirty_image, 
                         weighting, "","", szBaseOutFitsName );

   return ret;    
}

// TODO : required to convert from FITS to xcorr Visibility structure and avoid duplication of code 
/*Visibilities* CPacerImager::ConvertFits2XCorr( CBgFits& vis_re, CBgFits& vis_im )
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
   
   ObservationInfo obsInfo;
   obsInfo.nAntennas = n_ant;
   obsInfo.nFrequencies = 1;
   obsInfo.nPolarizations = 2;
   obsInfo.nTimesteps = 1;
   const size_t matrixSize = n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations;
   const size_t nIntervals  = (obsInfo.nTimesteps); // TODO  + voltages.nIntegrationSteps - 1) / voltages.nIntegrationSteps;
   unsigned int nChannelsToAvg = 1; // TODO : verify
   const size_t nOutFrequencies = obsInfo.nFrequencies / nChannelsToAvg;
   const size_t nValuesInTimeInterval = matrixSize * nOutFrequencies;
   const size_t outSize = nValuesInTimeInterval * nIntervals;
   unsigned int avgCh = fine_channel / nChannelsToAvg;
   
   std::complex<VISIBILITY_TYPE>* data = new std::complex<VISIBILITY_TYPE>[matrixSize]; // was n_ant*n_ant*obsInfo.nPolarizations
   Visibilities* xcorr = new Visibilities( data, obsInfo, 1, 1 );
   xcorr->nIntegrationSteps = 1;
   xcorr->nAveragedChannels = 1;
   xcorr->nFrequencies = 1;
   printf("DEBUG : CPacerImager::ConvertFits2XCorr - allocations ok, visibility array size = %d\n",int(matrixSize));fflush(stdout);


   for(int i=0;i<n_ant;i++){ // loop over ant1          
     for(int j=0;j<=i;j++){ // loop over ant2 
        std::complex<VISIBILITY_TYPE>* vis = xcorr->at( time_step, fine_channel, i, j );
        
        double re = vis_re.getXY( j , i );
        double im = vis_im.getXY( j , i );

        std::complex<VISIBILITY_TYPE> vis_value( re, im );
        (*vis) = vis_value;
  
// not required as only one side of the matrix is stored ?       
//        std::complex<VISIBILITY_TYPE>* vis_conj = xcorr->at( time_step, fine_channel, j , i );
//        (*vis_conj) = conj( vis_value );
     }
   }
   printf("DEBUG : CPacerImager::ConvertFits2XCorr end\n");

   return xcorr;
}*/


void CPacerImager::ConvertXCorr2Fits( Visibilities& xcorr, CBgFits& vis_re, CBgFits& vis_im, int time_step, int fine_channel, const char* szBaseFitsName )
{
   int n_ant = xcorr.obsInfo.nAntennas;
   int n_corrs = 4; // 4 correlation products : XX XY YX YY 
   int n_baselines = n_ant*(n_ant+1)/2;
   ObservationInfo& obsInfo = xcorr.obsInfo;
   const size_t matrixSize = n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations;
   const size_t nIntervals  = (obsInfo.nTimesteps); // TODO  + voltages.nIntegrationSteps - 1) / voltages.nIntegrationSteps;
   unsigned int nChannelsToAvg = 1; // TODO : verify
   const size_t nOutFrequencies = obsInfo.nFrequencies / nChannelsToAvg;
   const size_t nValuesInTimeInterval = matrixSize * nOutFrequencies;
   const size_t outSize = nValuesInTimeInterval * nIntervals;
   unsigned int avgCh = fine_channel / nChannelsToAvg;

   // size_t outIndex {interval * nValuesInTimeInterval + avgCh * matrixSize + baseline * obsInfo.nPolarizations * obsInfo.nPolarizations
   //                             + p1*obsInfo.nPolarizations + p2};

  // TODO !!!

  // assuming single timestep for now :  
  // assuming ordering of data as Cristian told me during the meeting :
  // Ant11        | Ant 12       | Ant 13       | ...
  // XX XY YX YY  | XX XY YX YY  | XX XY YX YY  | ...
  int index = 0;

  vis_re.SetNaN();
  vis_im.SetNaN();

  printf("DEBUG : CPacerImager::ConvertXCorr2Fits n_ant = %d\n",n_ant);
  // using at :
  // Complex<float> *at_float(Visibilities& vis, unsigned int interval, unsigned int frequency, unsigned int a1, unsigned a2)
  for(int i=0;i<n_ant;i++){ // loop over ant1          
     for(int j=0;j<=i;j++){ // loop over ant2 
        // auto& vis = xcorr.data[idx];
//        std::complex<double>* vis = at( xcorr, time_step, fine_channel, i, j );
        std::complex<VISIBILITY_TYPE>* vis = xcorr.at( time_step, fine_channel, i, j );

        vis_re.setXY(j,i,float(vis[0].real()));
        vis_im.setXY(j,i,float(vis[0].imag()));

        vis_re.setXY(i,j,float(vis[0].real()));
        vis_im.setXY(i,j,-float(vis[0].imag()));
     }
     index += (n_ant-i);
  }


  char szReFits[1024],szImFits[1024];
  sprintf(szReFits,"%s/%s_re.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseFitsName);
  sprintf(szImFits,"%s/%s_im.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseFitsName);
  vis_re.WriteFits( szReFits );
  vis_im.WriteFits( szImFits );
  printf("DEBUG : saved %s and %s\n",szReFits,szImFits);
}


bool CPacerImager::run_imager( Visibilities& xcorr, 
                               int time_step, 
                               int fine_channel,
                               double frequency_mhz, 
                               int n_pixels,
                               double FOV_degrees,
                               double min_uv /*=-1000*/,      // minimum UV
                               bool do_gridding /*=true*/,    // excute gridding  (?)
                               bool do_dirty_image /*=true*/, // form dirty image (?)
                               const char* weighting /*=""*/, // weighting : U for uniform (others not implemented)
                               const char* szBaseOutFitsName /*=NULL*/,
                               bool bCornerTurnRequired /*=true*/ // changes indexing of data "corner-turn" from xGPU structure to continues (FITS image-like)
                             )
{
  // ensures initalisation of object structures 
  Initialise();
  
  // m_ImagerParameters.m_ImageFOV_degree
  printf("DEBUG : frequency = %.6f [MHz] vs. m_ImagerParameters.m_fCenterFrequencyMHz = %.6f [MHz] , FoV = %.4f [deg] (CPacerImager::run_imager)\n",frequency_mhz,m_ImagerParameters.m_fCenterFrequencyMHz,FOV_degrees);
  
  int n_ant = xcorr.obsInfo.nAntennas;
  int n_pol = xcorr.obsInfo.nPolarizations;
  m_ImagerParameters.m_fCenterFrequencyMHz = frequency_mhz; // 2024-03-24 - otherwise everything was wrong - freq = 0 !!!
  m_ImagerParameters.m_ImageFOV_degrees = FOV_degrees;      // 2024-03-24 - reasons as above
  

/*  CBgFits fits_vis_real( n_ant, n_ant ), fits_vis_imag( n_ant, n_ant ), data_real_tmp( n_ant, n_ant ), data_imag_tmp( n_ant, n_ant );
  
  // temporary code for testing, later uncomment the line above:
  // step 1/ use line above to convert xcorr to FITS 
  // step 2/ remove ConvertXCorr2Fits altogether and propagate Visibility structure further down the stream 
  ConvertXCorr2Fits( xcorr, data_real_tmp, data_imag_tmp, time_step, fine_channel );
  float* data_real = data_real_tmp.get_data();
  float* data_imag = data_imag_tmp.get_data();
  
  
  if( bCornerTurnRequired )
  {
     fits_vis_real.SetValue(0.00);
     fits_vis_imag.SetValue(0.00);

     // ~corner-turn operation which can be quickly done on GPU: 
     int counter=0;
     for(int i=0;i<n_ant;i++)
     {
        for(int j=0;j<(i + 1);j++)
        {
           fits_vis_real.setXY( i, j , data_real[counter] );
           fits_vis_real.setXY( j, i , data_real[counter] );
           counter++;
        }
     }

     counter=0;
     for(int i=0;i<n_ant;i++)
     {
        for(int j=0;j<(i + 1);j++)
        {
           fits_vis_imag.setXY( i, j , data_imag[counter] );
           fits_vis_imag.setXY( j, i , -(data_imag[counter]) );
           counter++;
        }
     }
  }
  else
  {
     fits_vis_real.SetData( data_real );
     fits_vis_imag.SetData( data_imag );
  }*/
  
  // const char* get_filename(  time_t ut_time , char* out_buffer, int usec=0, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix=""
  if( m_ImagerParameters.m_fUnixTime <= 0.0001 ){
     m_ImagerParameters.m_fUnixTime = get_dttm_decimal();
     PRINTF_WARNING("Time of the data not specified -> setting current time %.6f\n",m_ImagerParameters.m_fUnixTime);
  }
  
/*  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_INFO ){
     char outDirtyImageReal[1024];
     // sprintf(outDirtyImageReal,"%s/%s_vis_real.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
     get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageReal, m_ImagerParameters.m_szOutputDirectory.c_str(), "visibility_", "_real" );
   
     fits_vis_real.WriteFits( outDirtyImageReal );
     PRINTF_DEBUG("Saved real file to %s\n",outDirtyImageReal);
  }
     
  if( CPacerImager::m_SaveFilesLevel >= SAVE_FILES_INFO ){
     char outDirtyImageImag[1024];
     // sprintf(outDirtyImageImag,"%s/%s_vis_imag.fits",m_ImagerParameters.m_szOutputDirectory.c_str(),szBaseOutFitsName);
     get_filename( m_ImagerParameters.m_fUnixTime, outDirtyImageImag, m_ImagerParameters.m_szOutputDirectory.c_str(), "visibility_", "_imag" );
        
     fits_vis_imag.WriteFits( outDirtyImageImag );
     PRINTF_DEBUG("Saved imag file to %s\n",outDirtyImageImag);
  }*/

  // calculate UVW (if required)
  CalculateUVW();
  
  if( m_ImagerParameters.m_bApplyGeomCorr ){
     if( IsGPU() ){
        printf("DEBUG : this is object CPacerImagerHIP -> geometric correction will be applied by the GPU Kernel (skipped here to avoid double-correction)\n");
     }else{
        printf("DEBUG : xcorr path test of ApplyGeometricCorrections\n");
        ApplyGeometricCorrections( xcorr, m_U, m_V, m_W, frequency_mhz, time_step, fine_channel  );
     }
  }
  
  if( m_ImagerParameters.m_bApplyCableCorr ){
     if( IsGPU() ){
        printf("DEBUG : this is object CPacerImagerHIP -> cable correction will be applied by the GPU Kernel (skipped here to avoid double-correction)\n");
     }else{
        printf("DEBUG : xcorr path test of ApplyCableCorrections\n");
        ApplyCableCorrections( xcorr, frequency_mhz, time_step, fine_channel  );
     }
  }
  
  if( IsGPU() ){
     printf("DEBUG : this is object CPacerImagerHIP -> cable/geo. corrections not executed here -> not saving control corr. matrix imager_test_vis_re.fits / imager_test_vis_im.fits\n");
     CBgFits re( xcorr.obsInfo.nAntennas, xcorr.obsInfo.nAntennas) ,im( xcorr.obsInfo.nAntennas, xcorr.obsInfo.nAntennas );
     ConvertXCorr2Fits( xcorr, re, im, time_step, fine_channel, "corrmatrix_before_run_imager_nocorrections_gpu" );
  }else{
     printf("DEBUG : saving control corr-matrix just before imaging and after applying (or not) corrections and cal (if required)\n");
     CBgFits re( xcorr.obsInfo.nAntennas, xcorr.obsInfo.nAntennas) ,im( xcorr.obsInfo.nAntennas, xcorr.obsInfo.nAntennas );
     ConvertXCorr2Fits( xcorr, re, im, time_step, fine_channel, "corrmatrix_after_cable_and_geo_corr_cpu" );
  }

  printf("DEBUG : just before run_imager(time_step=%d, fine_channel=%d )\n",time_step, fine_channel);fflush(stdout);  
  bool ret = // run_imager( fits_vis_real, fits_vis_imag, m_U, m_V, m_W,
             run_imager( xcorr, time_step, fine_channel, 
                         m_U, m_V, m_W,
                         frequency_mhz, 
                         n_pixels, 
                         FOV_degrees, 
                         min_uv,
                         do_gridding, 
                         do_dirty_image, 
                         weighting, "","", szBaseOutFitsName );

   return ret;    
}
                             
void CPacerImager::test()                             
{
   AllocGriddedVis(180,180);
}
                             