#ifndef _PACER_APPLY_CALIBRATION_H__
#define _PACER_APPLY_CALIBRATION_H__

// #include <math.h>

#include <vector>
#include <complex>
#include <string>

class CCalSol
{
public :
   CCalSol( int ant, std::complex<double> _xx, std::complex<double> _xy, std::complex<double> _yx, std::complex<double> _yy );
 
   int m_ant;
   std::complex<double> m_cal[4]; // 0 - xx , 1 - xy , 2 - yx , 3 - yy
};

// calibration solutions for all antennas :
class CCalSols : public std::vector<CCalSol> 
{
public :
   CCalSols();
   CCalSols(const CCalSols& right );
   
   CCalSols& operator=( const CCalSols& right );

   // reads calibration solutions for a single frequency channel from a text file with a structure :
   // # frequency channel 204 = 159.3750 MHz
   // # AMP_X PHASE_X AMP_Y PHASE_Y AMP_XY PHASE_XY AMP_YX PHASE_YX
   // see ~/aavs-calibration/station/calibration.py 
   int read_calsolutions_text( const char* filename , bool bForce=false, double sign_phase = -1 , bool bInverseAmp=false );   

   // read using filename from member variable m_filename    
   int read_calsolutions_text( double sign_phase = -1, bool bInverseAmp=false, bool bForce=false );
   
   // read using filename from member variable m_filename and chosing a different reading function according to file extension :
   int read_calsolutions( double sign_phase = -1, bool bInverseAmp=false, bool bForce=false );

   //---------------------------------------------------------------------------------------------------------
   // show calibration solutions 
   //---------------------------------------------------------------------------------------------------------
   void show();

   std::string m_filename;   
   int m_freq_channel;
   double m_frequency_mhz;
};

class CApplyCalibration
{
public :
   CApplyCalibration(){}
   ~CApplyCalibration(){}
};

#endif 
