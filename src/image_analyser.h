#ifndef _IMAGE_ANALYSER_H__
#define _IMAGE_ANALYSER_H__

class CPacerImager;

class CBgFits;

class CComplexImage
{
public :
   CComplexImage( int n_pixels );
   ~CComplexImage();
   
   void SetValue( float value );
   CComplexImage& operator=( CComplexImage& right );

   CBgFits* m_pImageReal;
   CBgFits* m_pImageImag;   
};


class CImageAnalyser
{
public :
   // frequency of the channel in MHz 
   double m_FreqMHz;
   int    m_Channel;
   int    m_nPixels; // number of pixels in the image
   double m_ThresholdSNR;

   CImageAnalyser( int n_pixels , double freq_mhz, int channel, bool bIncludeAutos=false );
   ~CImageAnalyser();
   
   CComplexImage* m_pImage1; // image buffer 1
   CComplexImage* m_pImage2; // image buffer 2
   
   // pointers to New and previous image:
   // previous time stamp images - for not just one last image is kept, maybe later it will be :
   //   a/ many or average
   //   b/ homeopatic average image 
   //   c/ just average image of N last images
   CComplexImage* m_pImageNew;
   CComplexImage* m_pImageOld;
   CBgFits* m_pDifferenceImage;

   // Imager to actually do images :
   CPacerImager* m_pPacerImager;
   
   //-------------------------------------------------------------------------------------------------------
   // Section with analysis functions :
   bool Image( CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, double unix_time=0.00 );
 
   // subtract m_pImageNew - m_pImageOld
   void SubtractImages();
   
   void AnalyseDifferenceImage();
   
   void SetParameters( const char* szAntennaPositionsFile, bool bConstantUVW=false, 
                       const char* szOutDir="./",
                       const char* szCalSolFile=NULL                       
                     );
                     
   void SetOutputDir( const char* szOutDir="./" );
};

#endif
