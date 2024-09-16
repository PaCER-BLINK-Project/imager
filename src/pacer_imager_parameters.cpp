#include "pacer_imager_parameters.h"
#include <string.h>

#include <limits>

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

