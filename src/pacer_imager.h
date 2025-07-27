#ifndef _PACER_IMAGER_H__
#define _PACER_IMAGER_H__

#define COTTER_COMPATIBLE true

#include "pacer_imager_parameters.h"
#include "observation_metadata.h"
#include "pacer_imager_defs.h"
#include "images.hpp"
#include <bg_fits.h>

#include <string>
using namespace std;

// FFTW, math etc :
// #include <fftw3.h>

// AstroIO for Visibilities class :
#include <astroio.hpp>
#include <memory_buffer.hpp>
// #define VISIBILITY_TYPE double // visibility type currently double but may be float 
// forward declaration of class  Visibilities, the header file is included in pacer_imager.cpp
// class  Visibilities;


// FFT shift to convert from DC in bin 0 of the FFTW output array to DC in the center bin :
void fft_shift(std::complex<float>* image, size_t image_x_side, size_t image_y_side);




class CPacerImager {

protected:
    MemoryBuffer<double> frequencies;
    MemoryBuffer<float> grids_counters;
    MemoryBuffer<std::complex<float>> grids;

public :
   // TODO: decide if this should be static or member variables
   // debug level for the whole library / program :
   static int m_ImagerDebugLevel; // see pacer_imager_defs.h for defines IMAGER_DEBUG_LEVEL etc
   
   // level of saving intermediate and test files , see pacer_imager_defs.h for defines SAVE_FILES_NONE
   static int m_SaveFilesLevel;
   
   // can also save control files every N-th file 
   static int m_SaveControlImageEveryNth;
   
   // show statistics of :
   //   - final imags
   //   - uv gridded visibilities 
   // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
   static bool m_bPrintImageStatistics;
   
   // TESTING and VALIDATION :
   static bool m_bCompareToMiriad; // should always be FALSE and set to TRUE only to test against MIRIAD 
   
   // include auto-correlations in the imaging :
   bool m_bIncludeAutos;

   // parameters :
   // WARNING : I am assuming they are the same for all objects of this class.
   //           we shall see and potentially remove "static"
   CImagerParameters m_ImagerParameters;

   // Antenna positions :   
//   CAntennaPositions m_AntennaPositions;
   
   // meta data :
   CObsMetadata m_MetaData;
   
   // Flagged antennas , if list m_AntennaPositions is filled it will also be updated (field flag)
   vector<int> m_FlaggedAntennas;
   
   // flag if initialisation already performed
   bool m_bInitialised;
   
   
   // UVW for SKA-Low station zenith phase-centered all-sky imaging :
   CBgFits m_U;
   CBgFits m_V;
   CBgFits m_W;
   double  u_mean, u_rms, u_min, u_max;
   double  v_mean, v_rms, v_min, v_max;
   double  w_mean, w_rms, w_min, w_max;

   int m_Baselines; // number of calculated baselines, also indicator if m_U, m_V and m_W have been initialised 
   int m_nAntennas;
   
   // MemoryBuffer<double> frequencies;
   
   // values calculated for the current image :
   double m_PixscaleAtZenith;
   

   bool     m_bLocalAllocation;


   CPacerImager(double unixTime, const std::string metadata_file);
   ~CPacerImager();
   void Initialise(double frequency_hz);
   
   void update_metadata(const std::string& metadata_file); // implement initialisation of object here, read antenna positions, calculate UVW if constant etc 
 
   // Set / Get functions :
   //-----------------------------------------------------------------------------------------------------------------------------
   // verbosity level 
   //-----------------------------------------------------------------------------------------------------------------------------
   static void SetDebugLevel( int debug_level );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // File save level
   //-----------------------------------------------------------------------------------------------------------------------------
   static void SetFileLevel( int filesave_level );
   
   
   void SetFlaggedAntennas( vector<int>& flagged_antennas);
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // read antenna positions and do whatever else is necessary (update flags etc)
   //-----------------------------------------------------------------------------------------------------------------------------
   int ReadAntennaPositions( bool bConvertToXYZ );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // update antenna flags in m_MetaData.m_AntennaPositons object based on m_FlaggedAntennas list 
   //-----------------------------------------------------------------------------------------------------------------------------
   int UpdateFlags();

   //-----------------------------------------------------------------------------------------------------------------------------
   // calculates UVW and also checks if it is required at all 
   // U,V,W are calculated in meters (not wavelengths) -> have to be divided later by wavelength
   //-----------------------------------------------------------------------------------------------------------------------------
   bool CalculateUVW(double frequency_hz, bool bForce=false, bool bInitialise=true);
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Re-calculate parameters based on new UVW :
   //-----------------------------------------------------------------------------------------------------------------------------
   bool UpdateParameters(double frequency_hz);
   

   //-----------------------------------------------------------------------------------------------------------------------------
   // 1st version producing a dirty image (tested on both MWA and SKA-Low).
   // TODO : Test cases can be found in PaCER documentation 
   //-----------------------------------------------------------------------------------------------------------------------------
   void dirty_image(MemoryBuffer<std::complex<float>>& grids, MemoryBuffer<float>& grids_counters,
     int grid_side, int n_integration_intervals, int n_frequencies, MemoryBuffer<std::complex<float>>& images_buffer);

   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 
   // GOAL   : fast version of gridding code which does not require fft_unshift
   //
   // INPUT  : 
   //          fits_vis_real, fits_vis_imag : visibilities (REAL and IMAG 2D arrays as FITS class) 
   //          fits_vis_u, fits_vis_v, fits_vis_w : UVW (real values baselines in units of wavelength - see TMS)
   //          delta_u, delta_v : size of the UV cell 
   //          frequency_mhz : frequency in MHz
   //
   // OUTPUT : 
   //          - uv_grid_real, uv_grid_imag : visibilities on UV grid (real and imag arrays)
   //          - uv_grid_counter : visibility counter and 
   //-----------------------------------------------------------------------------------------------------------------------------


   // same as above 
   // TODO : only use the new one gridding_fast( Visibilities& xcorr, ... ) and use ConvertFits2XCorr in gridding_fast( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, ... ) to
   //        call this one :
   void gridding_fast( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  MemoryBuffer<std::complex<float>>& grids, MemoryBuffer<float>& grids_counters, double delta_u, double delta_v,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );

   
   // same as above but using AstroIO Visibility as input
   // TODO : at some point this one should stay and the other one should be a wrapper using ConvertFits2XCorr and calling this one:
   virtual Images gridding_imaging( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,                   
                  double delta_u, double delta_v,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL
                );



   Images run_imager( Visibilities& xcorr, 
                    int time_step, 
                    int fine_channel,
                    int n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,      // minimum UV
                    bool do_gridding=true,    // excute gridding  (?)
                    bool do_dirty_image=true, // form dirty image (?)
                    const char* weighting="", // weighting : U for uniform (others not implemented)
                    const char* szBaseOutFitsName=NULL,
                    bool bCornerTurnRequired=true // changes indexing of data "corner-turn" from xGPU structure to continues (FITS image-like)
                  );


   //-----------------------------------------------------------------------------------------------------------------------------
   // Apply phase corrections : geometrical correction, cable correction etc   
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual void ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_w, MemoryBuffer<double>& frequencies);
   
   virtual void ApplyCableCorrections( Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies);

   //-----------------------------------------------------------------------------------------------------------------------------
   // Function saving output FITS files: 
   //-----------------------------------------------------------------------------------------------------------------------------
   bool SaveSkyImage( const char* outFitsName, CBgFits* pFits, double unixtime=0.00 );

   double get_frequency_hz(const Visibilities& xcorr, int fine_channel, bool cotter_compatible);           
};

#endif 
