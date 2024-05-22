#include <stdlib.h>
#include <stdio.h>
#include <myutil.h>
#include "realtime_image_analyser.h"
#include "image_analyser.h"
#include "pacer_imager.h"
#include "pacer_imager_defs.h"

CRealTimeImageAnalyser::CRealTimeImageAnalyser(int n_channels, int n_pixels, double center_freq_mhz, double bw_mhz, bool bIncludeAutos /*=false*/ )
{
   double channel_bw_mhz = bw_mhz/n_channels;
   double start_freq_mhz = center_freq_mhz - bw_mhz/2.00;
   double first_channel_center_freq_mhz = start_freq_mhz + channel_bw_mhz/2.00;

   for(int i=0;i<n_channels;i++){
      double freq_mhz = first_channel_center_freq_mhz + i*channel_bw_mhz;
   
      CImageAnalyser* pNewAnalyser = new CImageAnalyser( n_pixels, freq_mhz, i, bIncludeAutos  );      
      
      // update parameters :
      if( pNewAnalyser->m_pPacerImager ){
         pNewAnalyser->m_pPacerImager->m_ImagerParameters.m_fCenterFrequencyMHz = center_freq_mhz; 
      }

      m_ChannelImagers.push_back( pNewAnalyser );
      printf("DEBUG : initialised imager for channel = %d\n",i);
   }
   
   // TODO : just temporary for testing : only image channel 0 
//   m_SelectedChannels.push_back(0);
}

CRealTimeImageAnalyser::~CRealTimeImageAnalyser()
{
   for(int i=0;i<m_ChannelImagers.size();i++){
      CImageAnalyser* pImageAnalyser = m_ChannelImagers[i];
      
      delete pImageAnalyser;
   }
}

CImageAnalyser* CRealTimeImageAnalyser::GetChannelImager( int channel )
{
   if( channel<0 || channel >= m_ChannelImagers.size() ){
      return NULL;
   }

   return m_ChannelImagers[channel];
}

void CRealTimeImageAnalyser::SetOutputDir( const char* szOutDir )
{
   for(int i=0;i<m_ChannelImagers.size();i++){
      CImageAnalyser* pImageAnalyser = m_ChannelImagers[i];
      char szFullOutDir[1024];
//      sprintf(szFullOutDir,"%s/%05d/",szOutDir,i);
      sprintf(szFullOutDir,"%s/",szOutDir);
   
      pImageAnalyser->SetOutputDir( szFullOutDir );
   }
}

void CRealTimeImageAnalyser::SetParameters( const char* szAntennaPositionsFile, bool bConstantUVW /*=false*/ , const char* szOutDir /*="./"*/, const char* szCalSolFile /*=NULL*/ )
{
   for(int i=0;i<m_ChannelImagers.size();i++){
      CImageAnalyser* pImageAnalyser = m_ChannelImagers[i];
      char szFullOutDir[1024];
      sprintf(szFullOutDir,"%s/%05d/",szOutDir,i);
//      sprintf(szFullOutDir,"%s/",szOutDir);
   
      pImageAnalyser->SetParameters( szAntennaPositionsFile, bConstantUVW, szFullOutDir, szCalSolFile ); // 2023-04-25 : this is now set in Image function
   }
}

bool CRealTimeImageAnalyser::Image( int channel, CBgFits& corr_matrix_real, CBgFits& corr_matrix_imag, double unix_time /*=0.00*/ )
{
   if( channel >= 0 && channel < m_ChannelImagers.size() ){
      bool bProcessChannel=true;
      if( m_SelectedChannels.size() > 0 ){
         if( find_value( m_SelectedChannels, channel ) < 0 ){
            bProcessChannel=false;
         }
      }
   
   
      if( bProcessChannel ){      
         clock_t start = clock();
         (m_ChannelImagers[channel])->Image( corr_matrix_real, corr_matrix_imag, unix_time );
         printf("INFO : imaging channel %d is required and took %ld [usec]\n",channel,(clock()-start));
      }else{
         printf("WARNING : channel %d is not required to be imaged\n",channel);
      }
   }else{
      printf("ERROR in code : cannot image channel %d , CRealTimeImageAnalyser initialised for %d channels\n",channel,int(m_ChannelImagers.size()));
   }
   
   return false;
}


// set flagged antennas :
void CRealTimeImageAnalyser::SetFlaggedAntennas( vector<int>& flagged_antennas)
{
   for(int i=0;i<m_ChannelImagers.size();i++){
      CImageAnalyser* pImageAnalyser = m_ChannelImagers[i];
      
      if( pImageAnalyser && pImageAnalyser->m_pPacerImager ){
         pImageAnalyser->m_pPacerImager->SetFlaggedAntennas( flagged_antennas );
      }else{
         PRINTF_ERROR("Error in code CRealTimeImageAnalyser::SetFlaggedAntennas one of the pointers is not set !\n");         
      }
   }
}
