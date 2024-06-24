#include <stdlib.h>
#include <myfile.h>
#include <mystring.h>
#include <bg_fits.h>
#include "image_analyser.h"

#include "pacer_imager.h"
#include "pacer_imager_defs.h"
#include "pacer_common.h"

CComplexImage::CComplexImage( int n_pixels )
{
   m_pImageReal = new CBgFits( n_pixels, n_pixels );
   m_pImageImag = new CBgFits( n_pixels, n_pixels );
}

CComplexImage::~CComplexImage()
{
   if( m_pImageReal ){
      delete m_pImageReal;
   }
   if( m_pImageImag ){
      delete m_pImageImag;
   }
}

void CComplexImage::SetValue( float value )
{
   m_pImageReal->SetValue( value );
   m_pImageImag->SetValue( value );
}

CComplexImage& CComplexImage::operator=( CComplexImage& right )
{
   if( m_pImageReal->GetXSize() != right.m_pImageReal->GetXSize() || m_pImageReal->GetYSize() != right.m_pImageReal->GetYSize() ){
       printf("ERROR in code CComplexImage::operator= Real image sizes are not the same (%ld,%ld) != (%ld,%ld)\n",m_pImageReal->GetXSize(),m_pImageReal->GetYSize(),right.m_pImageReal->GetXSize(),right.m_pImageReal->GetYSize());
       throw string("ERROR in code CComplexImage::operator= Real image sizes are not the same");
   }

   if( m_pImageImag->GetXSize() != right.m_pImageImag->GetXSize() || m_pImageImag->GetYSize() != right.m_pImageImag->GetYSize() ){
       printf("ERROR in code CComplexImage::operator= Imag image sizes are not the same (%ld,%ld) != (%ld,%ld)\n",m_pImageImag->GetXSize(),m_pImageImag->GetYSize(),right.m_pImageImag->GetXSize(),right.m_pImageImag->GetYSize());
       throw string("ERROR in code CComplexImage::operator= Imag image sizes are not the same");
   }

   (*m_pImageReal) = (*right.m_pImageReal);
/*   for(int y=0;y<m_pImageReal->GetYSize();y++){
      for(int x=0;x<m_pImageReal->GetXSize();x++){
         m_pImageReal->setXY( x, y, m_pImageReal->getXY(x,y) );
      }
   }*/

   (*m_pImageImag) = (*right.m_pImageImag);
/*   for(int y=0;y<m_pImageImag->GetYSize();y++){
      for(int x=0;x<m_pImageImag->GetXSize();x++){
         m_pImageImag->setXY( x, y, m_pImageImag->getXY(x,y) );
      }
   }*/
   
   return (*this);
}

CImageAnalyser::CImageAnalyser( int n_pixels , double freq_mhz, int channel, bool bIncludeAutos /*=false*/  )
: m_pPacerImager(NULL), m_pImageNew(NULL), m_pImageOld(NULL), m_FreqMHz(freq_mhz), m_Channel(channel), m_nPixels(n_pixels), m_ThresholdSNR(5.00)
{
   m_pImage1 = new CComplexImage( n_pixels );
   m_pImage2 = new CComplexImage( n_pixels );
   m_pDifferenceImage = new CBgFits( n_pixels , n_pixels );

   m_pPacerImager = new CPacerImager();
   if( bIncludeAutos ){
      m_pPacerImager->m_bIncludeAutos = bIncludeAutos;
   }
}


CImageAnalyser::~CImageAnalyser()
{
   if( m_pPacerImager ){
      delete m_pPacerImager;
   }

   if( m_pImage1 ){
      delete m_pImage1;
   }
   if( m_pImage2 ){
      delete m_pImage2;
   }
   
   if( m_pDifferenceImage ){
      delete m_pDifferenceImage;
   }
}

void CImageAnalyser::SetOutputDir( const char* szOutDir /*="./"*/ )
{
   (m_pPacerImager->m_ImagerParameters).m_szOutputDirectory = szOutDir;
}

void CImageAnalyser::SetParameters( const char* szAntennaPositionsFile, bool bConstantUVW /*=false*/ , const char* szOutDir /*="./"*/, const char* szCalSolFile /*=NULL*/ )
{
   (m_pPacerImager->m_ImagerParameters).SetGlobalParameters( szAntennaPositionsFile, bConstantUVW );

   // set default values specific to this program :
   CPacerImager::m_SaveFilesLevel = SAVE_FILES_ALL; // SAVE_FILES_NONE
   // but save every 10th files (1 per second) :
   // TODO : in the future add saving mean 1-second file too
   CPacerImager::m_SaveControlImageEveryNth = 10;

   // set output directory :   
   // char szCH_Dir[256];
   // sprintf(szCH_Dir,"%s/ch05%d/",szOutDir,m_Channel);
   // SetOutputDir( szCH_Dir );
   (m_pPacerImager->m_ImagerParameters).m_szOutputDirectory = szOutDir;   
   
   if( szCalSolFile && strlen(szCalSolFile) ){
      (m_pPacerImager->m_CalibrationSolutions).m_filename = szCalSolFile;
      
      if( !(m_pPacerImager->m_CalibrationSolutions).read_calsolutions() ){
         printf("ERROR : could not read calibration solutions from files %s !!!\n",szCalSolFile);
         exit(-1);
      }
   }   
}

bool CImageAnalyser::Image( CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, double unix_time /*=0.00*/ )
{
   bool bFirstImage = false;
   
   if( !m_pImageNew ){
      m_pImageNew = m_pImage1;      
      m_pImageOld = m_pImage2;
      m_pImage1->SetValue(0.00);
      m_pImage2->SetValue(0.00);
      
      bFirstImage = true;
   }else{
      // Swap NEW AND OLD pointers to avoid copying data:
      CComplexImage* pTmpPtr = m_pImageOld;
      m_pImageOld = m_pImageNew;
      m_pImageNew = pTmpPtr;
   }
      
   if( !m_pPacerImager ){
      printf("ERROR in code (CImageAnalyser::Image) : pointer m_pPacerImager not set -> cannot continue\n");
      return false;
   }
   
   if( !m_pImageNew ){
      printf("ERROR in code (CImageAnalyser::Image) : pointer m_pImageNew not set -> cannot continue\n");
      return false;
   }
   
   if( !m_pImageNew->m_pImageReal ){
      printf("ERROR in code (CImageAnalyser::Image) : pointer m_pImageNew->m_pImageReal not set -> cannot continue\n");
      return false;      
   }
   
   if( !m_pImageNew->m_pImageImag ){
      printf("ERROR in code (CImageAnalyser::Image) : pointer m_pImageNew->m_pImageImag not set -> cannot continue\n");
      return false;
   }
      
   // set buffer for output image :
   m_pPacerImager->SetOutputImagesExternal( m_pImageNew->m_pImageReal, m_pImageNew->m_pImageImag );
      
   int n_ant = corr_matrix_real.GetXSize(); // size of correlation matrix == number of antennas 
   int n_pixels = m_pImageNew->m_pImageReal->GetXSize();  // number of pixels as the array size for image
   double FOV_deg = 180.00;
   
   // enable test MIRIAD parameters :
   CPacerImager::m_bCompareToMiriad = true;
   // printf("WARNING : changing CPacerImager::m_bCompareToMiriad = true;\n");
   
   // set unixtime parameter:
   (m_pPacerImager->m_ImagerParameters).SetUnixTime( unix_time );
   
   // apply calibration if provided :
//   if( strlen( (m_pPacerImager->m_CalibrationSolutions).m_filename.c_str() ) ){
//      PRINTF_DEBUG("DEBUG : applying calibration as in file %s\n",(m_pPacerImager->m_CalibrationSolutions).m_filename.c_str());
//      m_pPacerImager->ApplySolutions( corr_matrix_real, corr_matrix_imag, m_FreqMHz, "X" );
//   }

/*   char szBaseOutFitsName[256],szCH[16],szCH_Dir[256];
   sprintf(szCH,"_ch%05d",m_Channel);
   sprintf(szCH_Dir,"%s/ch05%d/",(m_pPacerImager->m_ImagerParameters).m_szOutputDirectory.c_str(),m_Channel);
   get_filename_base( unix_time, szBaseOutFitsName, szCH_Dir, "dirty_image_", szCH );*/
   m_pPacerImager->run_imager( corr_matrix_real.get_data(), corr_matrix_imag.get_data(), n_ant, 2, m_FreqMHz, n_pixels, FOV_deg, -1000, true, true, "N", NULL, false );  // Natural or Uniform ?
   
   
   if( bFirstImage ){
      // set m_pImageOld to the same value, there is no operator= in CBgFits yet, so use += operator (both initialised to ZERO so will be the same);
      (*m_pImage2) = (*m_pImage1);
     
      printf("INFO : first image, difference image analysis skipped\n");
   }else{   
      // analyse images (calculate difference image etc)
      AnalyseDifferenceImage();
   }

   return true;
}

void CImageAnalyser::SubtractImages()
{
   for(int y=0;y<m_nPixels;y++){
      for(int x=0;x<m_nPixels;x++){
         double diff = m_pImageNew->m_pImageReal->getXY(x,y) - m_pImageOld->m_pImageReal->getXY(x,y);
      
         m_pDifferenceImage->setXY( x, y, diff );   
      }
   }
}

void CImageAnalyser::AnalyseDifferenceImage()
{
   SubtractImages();

   // double GetStat( double& mean, double& rms, double& minval, double& maxval, int x_start=0, int y_start=0, int x_end=-1, int y_end=-1 );
   double mean,rms,minval,maxval;
   int x_start = m_nPixels/2 - 10;
   int x_end = m_nPixels/2 + 10;
   int y_start = x_start;
   int y_end = x_end;
   
   m_pDifferenceImage->GetStat( mean, rms, minval, maxval, x_start, y_start, x_end, y_end );
   printf("STAT in window (%d,%d)-(%d,%d) : channel %.4f MHz image, has mean = %.6f [Jy], rms = %.6f [Jy], minval = %.6f [Jy], maxval = %.6f [Jy]\n",x_start,y_start,x_end,y_end,m_FreqMHz,mean,rms,minval,maxval);

   char szLogFilename[1024];
   int freq_hz = int(m_FreqMHz*1000000.00);
   sprintf(szLogFilename,"%s/candidates_%06dHz.txt",(m_pPacerImager->m_ImagerParameters).m_szOutputDirectory.c_str(),freq_hz);
   
   bool bNewFile=true;
   if( MyFile::DoesFileExist(szLogFilename) ){
      bNewFile = false;
   }
   
   FILE* out_f = fopen(szLogFilename,"a+");
   if( bNewFile ){
      fprintf(out_f,"# FILENAME X Y SNR F[Jy] Thresh[Jy]\n");
      fprintf(out_f,"# Frequency = %.3f MHz\n",m_FreqMHz);
   }
   
   int n_transient_candidates = 0;
   double threshold = mean + m_ThresholdSNR*rms;
   for(int y=0;y<m_nPixels;y++){
      for(int x=0;x<m_nPixels;x++){
         double val = m_pDifferenceImage->getXY(x,y);
         
         if( val > threshold ){
            double snr = (val-mean)/rms;
            
            printf("CAND at %.3f MHz %s : %d %d SNR = %.3f , flux = %.6f [Jy] > %.6f\n",m_FreqMHz,m_pImageNew->m_pImageReal->GetFileName(),x,y,snr,val,threshold);
            fprintf(out_f,"%s %d %d %.3f %.6f %.6f\n",m_pImageNew->m_pImageReal->GetFileName(),x,y,snr,val,threshold);
            
            n_transient_candidates++;
         }
      }
   }
   printf("STAT : number of candidates in different image exceeding threshold %.3f is %d ( %s )\n",threshold,n_transient_candidates,m_pImageNew->m_pImageReal->GetFileName());
   
   if( n_transient_candidates ){
      char szDiffOutFile[256];
      // mystring getfname( const mystring& name, mystring& szExt)
      mystring szFullFileName = m_pImageNew->m_pImageReal->GetFileName();
      mystring szFileNameBase,szExt;
      szFileNameBase = getfname( szFullFileName, szExt );
      
      sprintf(szDiffOutFile,"%s/diff/%s.fits", (m_pPacerImager->m_ImagerParameters).m_szOutputDirectory.c_str() , szFileNameBase.c_str() );
      if( m_pDifferenceImage->WriteFits( szDiffOutFile ) ){
         printf("WARNING : could not write output different image %s\n",szDiffOutFile);
      }else{
         printf("INFO : difference image saved to file %s ( %d transient candidates )\n",szDiffOutFile,n_transient_candidates);
      }
   }
   
   fclose(out_f);
}

