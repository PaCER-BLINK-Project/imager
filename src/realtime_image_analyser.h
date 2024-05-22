#ifndef _REALTIME_IMAGE_ANALYSER_H__
#define _REALTIME_IMAGE_ANALYSER_H__

#include <vector>
using namespace std;

// homeopatic formula :
// n images
// alfa
// I_i - i-th image
// R_i - i-th reference image
// R_n = I_0(1-alfa)^n + Sum_i=0_to_i=(n-1) (1-alfa)^i I_(n-1)
// (1-alfa)^n <= Epsilon -> alfa >= 1 - Epsilon^(1/n)
//

class CImageAnalyser;
class CBgFits;

class CRealTimeImageAnalyser 
{
public :
   // current images per channel (not historical):
   // m_ChannelImagers.size() == number of channels (no need to keep another int m_nChannels)
   vector<CImageAnalyser*> m_ChannelImagers;
   vector<int> m_SelectedChannels;
         
   CRealTimeImageAnalyser(int n_channels,int n_pixels, double center_freq_mhz, double bw_mhz, bool bIncludeAutos=false);
   ~CRealTimeImageAnalyser();
   
   // set parameters :
   void SetParameters( const char* szAntennaPositionsFile, bool bConstantUVW=false, const char* szOutDir="./", const char* szCalSolFile=NULL );
   void SetOutputDir( const char* szOutDir="./" );
      
   // imager :
   bool Image( int channel, CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, double unix_time=0.00 );
   
   // return imager object for a specific frequency channel : 
   CImageAnalyser* GetChannelImager( int channel );

   // set flagged antennas :
   void SetFlaggedAntennas( vector<int>& flagged_antennas );
};

#endif
