#ifndef _PACER_IMAGER_PARAMETERS_H__
#define _PACER_IMAGER_PARAMETERS_H__

// #include "antenna_positions.h"
// #include "apply_calibration.h"
// #include <bg_fits.h>

#include <string>

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
   
   // UV range :
   double m_fMinUV;
   
   // output image size :
   int m_ImageSize;
   
   // Image Field-of-View :
   double m_ImageFOV_degrees;
   double m_PixsizeInRadians; // pixel size in radians 
   
   std::string m_szWeighting; // U for uniform, N for natural etc 
   
   // W-range 
   // bool m_bW_range_specified;
   double m_MinW; // -INF
   double m_MaxW; // +INF
   

   CImagerParameters();
   CImagerParameters( const CImagerParameters& right );
   static void SetGlobalParameters( const char* szAntennaPositionsFile, bool bConstantUVW=false );
   CImagerParameters& operator=( const CImagerParameters& right );
   
   inline void SetConstantUVW(){ m_bConstantUVW = true; }
   inline void SetUnixTime( double fUnixTime ){ m_fUnixTime = fUnixTime; printf("DEBUG : unixtime set to %.6f\n",m_fUnixTime); }
   
//   void SetAntPositionsFile( const char* szAntennaPositionsFile );
   
};

#endif 
