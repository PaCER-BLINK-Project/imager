#ifndef _OBSERVATION_METADATA_H__
#define _OBSERVATION_METADATA_H__

#include <vector>
#include <cstring>
#include <string>
#include <map>
#include <fitsio.h>

#include "antenna_positions.h"
#include "pacer_imager_defs.h"

// struct fitsfile;

// based on mwaconfig.h(cpp) and metafitsfile.h(cpp) in cotter 

#define DEFUALT_STRING_VALUE "Unknown"

/*
#define MWA_LATTITUDE -26.703319        // Array latitude. degrees North
#define MWA_LONGITUDE 116.67081         // Array longitude. degrees East
#define MWA_HEIGHT 377               // Array altitude. meters above sea level
*/


class CObsMetadata 
{
public :
   CObsMetadata( const char* filename=DEFUALT_STRING_VALUE );
   
   std::string m_filename; // Meta data file (txt,metafits etc).
   
   // 
   fitsfile *_fptr;
   bool m_bHasMetaFits;
   
   // Antenna positions :   
   CAntennaPositions m_AntennaPositions;

   size_t nInputs;                 // the number of inputs to the correlator
   size_t nScans;                  // the number of time samples
   size_t nChannels;               // the number of spectral channels
   enum CorrType { None, CrossCorrelation, AutoCorrelation, BothCorrelations } correlationType;
   double integrationTime;         // per time sample, in seconds
   double centralFrequencyMHz;     // observing central frequency and 
   double bandwidthMHz;            // bandwidth (MHz)
   double haHrs, raHrs, decDegs;          // ra,dec of phase centre.
   double haHrsStart;              // the HA of the phase center at the start of the integration
//   double refEl, refAz;            // NOT USED / TO BE REMOVED : the el/az of the normal to the plane of the array (radian)
   double AzimDeg;                 // Pointing Azimuth - used
   double ElevDeg;                 // Pointing Elevatoion - used
   int    year, month, day;        // date/time in UTC.
   int    refHour, refMinute;
   double refSecond;
   //bool   invertFrequency;         // flag to indicate that freq decreases with increasing channel number.
   bool   conjugate;               // conjugate the final vis to correct for any sign errors
   bool   geomCorrection;          // apply geometric phase correction
   std::string fieldName;
   std::string polProducts;

   // time :
   double startUnixTime;
   double dateFirstScanMJD;        // The central time of the first time step. Not in header file, derived from it.
   
   // array location:
   double geo_long;
   double geo_lat;
   double geo_height;
   
   // part from MWAHeaderExt
   int gpsTime;
   std::string observerName, projectName, gridName, mode, mwaPyVersion, mwaPyDate, metaDataVersion;
   std::string filename; // name of file in the METADATA
   int delays[16], subbandGains[24], subbandNumbers[24];
   bool hasCalibrator, hasGlobalSubbandGains;
   int centreSBNumber;
   //double fiberFactor;
   double tilePointingRARad, tilePointingDecRad;
   double dateRequestedMJD;
   
   // Telescope:
   std::string szTelescopeName;
   eTelescopName_Type eTelescopeName;
   
   bool HasMetaFits(){
      return m_bHasMetaFits;
   }

   double GetDateLastScanMJD() const
   {
      return dateFirstScanMJD + (integrationTime/86400.0)*nScans;
   }
   
   // TODO :
//   double GetStartDateMJD() const;
   double GetDateFirstScanFromFields() const;
   
   //------------------------------------------------------------------------------------------------------------------------------
   // Read metadata 
   // reading meta data from text file : 
   bool ReadMetaDataTxt( const char* filename );

   // if obsid > 0 is passed RA,DEC and other metadata are updated automatically    
   bool ReadMetaFitsFile( const char* filename);

   // if obsid > 0 is passed RA,DEC and other metadata are updated automatically    
   bool ReadMetaData( const char* filename);
   
   // reading of antenna positions :
   int ReadAntPositions();

   // parsing Metafits file information :
   bool parseKeyword( const std::string& keyName, const std::string& keyValue );
   bool parseFitsString(const char* valueStr, std::string& outString );
   std::string stripBand(const std::string& input);
   bool parseFitsDate(const char* valueStr, int& year, int& month, int& day, int& hour, int& min, double& sec);
   double parseFitsDateToMJD(const char* valueStr);
   bool parseIntArray(const char* valueStr, int* values, size_t count);
   bool parseBool(const char* valueStr);

   // error checking :
   bool checkStatus(int status, const char* szErrorMsg );
   
   // updating pointing information and time based on current time and (Az,Elev) pointing :
   // based on python scripts in SMART imaging pipeline : see ~/github/smart/scripts/fix_metafits_time_radec_all.py
   bool fix_metafits( double obsid, double inttime_sec=1.00 ); // OBSID is in fact GPS time
   

   // auxiliary functions :
   static double ux2gps( double uxtime );   
   static double gps2ux( double gpstime );

   // void Validate(bool lockPointing);
   private:
      CObsMetadata(const CObsMetadata &) { }
      void operator=(const CObsMetadata &) { }
      
      static bool isDigit(char c) { return c >= '0' && c <= '9'; }
};


#endif
