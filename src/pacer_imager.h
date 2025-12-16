#ifndef _PACER_IMAGER_H__
#define _PACER_IMAGER_H__

#define COTTER_COMPATIBLE true

#include "observation_metadata.h"
#include "pacer_imager_defs.h"
#include <images.hpp>

#include <string>
using namespace std;

// FFTW, math etc :
// #include <fftw3.h>

// AstroIO for Visibilities class :
#include <astroio.hpp>
#include <memory_buffer.hpp>

#include "gridding.hpp"
// #define VISIBILITY_TYPE double // visibility type currently double but may be float 
// forward declaration of class  Visibilities, the header file is included in pacer_imager.cpp
// class  Visibilities;


// FFT shift to convert from DC in bin 0 of the FFTW output array to DC in the center bin :
void fft_shift(std::complex<float>* image, size_t image_x_side, size_t image_y_side);


class CPacerImager {

protected:
   MemoryBuffer<double> cable_lengths;
   MemoryBuffer<double> frequencies;
   MemoryBuffer<float> grids_counters;
   MemoryBuffer<std::complex<float>> grids;
   size_t n_gridded_channels {0};
   size_t n_gridded_intervals {0};
   MemoryBuffer<float> u_cpu;
   MemoryBuffer<float> v_cpu;
   MemoryBuffer<float> w_cpu;

public :
   // TODO: decide if this should be static or member variables
   // debug level for the whole library / program :
   static int m_ImagerDebugLevel; // see pacer_imager_defs.h for defines IMAGER_DEBUG_LEVEL etc
   
   // level of saving intermediate and test files , see pacer_imager_defs.h for defines SAVE_FILES_NONE
   static int m_SaveFilesLevel;
   
   // show statistics of :
   //   - final imags
   //   - uv gridded visibilities 
   // TODO : ? may require a separate flag in the future, for now just using a single Statistics switch ON/OFF flag
  static bool m_bPrintImageStatistics;
   
   // include auto-correlations in the imaging :
   bool m_bIncludeAutos {false};

  bool autofix_metadata {true};  // automatically recalculate RA,DEC,TIME using just standard METAFITS file (no need to create special METAFITS using fix_metafits_time_radec_all.py )
  bool constant_uvw; // default false and only true to zenith phase centered images (all-sky from EDA2)
  bool apply_geom_correction {true};
  bool apply_cable_correction {true};
  bool average_images {false};
  float oversampling_factor {2.0f};

   // meta data :
   CObsMetadata m_MetaData;
   
   // Flagged antennas , if list m_AntennaPositions is filled it will also be updated (field flag)
   vector<int> m_FlaggedAntennas;
   vector<int> m_FlaggedBaselines;
   
   // UVW for SKA-Low station zenith phase-centered all-sky imaging :
   double u_min {0}, u_max {0};
   double v_min {0}, v_max {0};
   double w_min {0}, w_max {0};
   double delta_u {0}, delta_v {0};

   int m_Baselines;
   int n_pixels;
   double min_uv {-1000};
   const char *weighting = "";
   // MemoryBuffer<double> frequencies;
   
   // values calculated for the current image :
   double m_PixscaleAtZenith {0};
   double pixsize_in_radians {0};
   
   Polarization pol_to_image {Polarization::XX};
   // TEMPORARY SOLUTION for imaging single frequency:
   // single frequency testing so that it is possible to pass frequency externally and do not calculate from channel (one channel per freq)
   // TODO : later use array of frequencies in MHz instead of channels to avoid calculation of frequency inside imager (get list as parameter)
   double m_fFrequencyMHz {-1};
   

   CPacerImager(const std::string metadata_file, int n_pixels, const std::vector<int>& flagged_antennas, bool average_images = false,
      Polarization pol_to_image = Polarization::XX, float oversampling_factor = 2.0f, double min_uv=-1000, const char* weighting="");
   
   // Set / Get functions :
   //-----------------------------------------------------------------------------------------------------------------------------
   // verbosity level 
   //-----------------------------------------------------------------------------------------------------------------------------
   static void SetDebugLevel( int debug_level );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // File save level
   //-----------------------------------------------------------------------------------------------------------------------------
   static void SetFileLevel( int filesave_level );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // update antenna flags in m_MetaData.m_AntennaPositons object based on m_FlaggedAntennas list 
   //-----------------------------------------------------------------------------------------------------------------------------
   int UpdateFlags();

   //-----------------------------------------------------------------------------------------------------------------------------
   // calculates UVW and also checks if it is required at all 
   // U,V,W are calculated in meters (not wavelengths) -> have to be divided later by wavelength
   //-----------------------------------------------------------------------------------------------------------------------------
   bool CalculateUVW();

   virtual void gridding(Visibilities& xcorr);

/** 
    @brief run gridding + fft to produce an image.
    
    @param xcorr: Visibilities to be imaged.
*/
   Images run(Visibilities& xcorr);
   void grid(Visibilities& xcorr);
   virtual Images image(ObservationInfo& obsInfo);


   //-----------------------------------------------------------------------------------------------------------------------------
   // Apply phase corrections : geometrical correction, cable correction etc   
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual void ApplyGeometricCorrections( Visibilities& xcorr, MemoryBuffer<float>& w, MemoryBuffer<double>& frequencies);
   
   virtual void ApplyCableCorrections( Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies);

   double get_frequency_hz(const Visibilities& xcorr, int fine_channel, bool cotter_compatible);           


   //-----------------------------------------------------------------------------------------------------------------------------
   // START of functions and members for convolution kernel - taken from the msok_devel branch which still uses CBgFits class etc   
   //-----------------------------------------------------------------------------------------------------------------------------
   // GOAL : convolve with a specified 
   //-----------------------------------------------------------------------------------------------------------------------------
   void convolve_with_kernel( std::complex<float>* current_grid, int x_grid, int y_grid, double re, double im, double u, double v, double delta_u, double delta_v, int n_pixels, double im_sign, int oversampling=10 );
//   void convolve_with_kernel( CBgFits& uv_grid_real, CBgFits& uv_grid_imag, int x_grid, int y_grid, double re, double im, double u, double v, double delta_u, double delta_v, int n_pixels, double im_sign, int oversampling=10 );
   
   //-------------------------------------------------------------------------------------------------------------
   // Anti-aliasing and kernel convolution member functions. Currently taken from WSCLEAN code - to be re-implemented 
   //-------------------------------------------------------------------------------------------------------------
   static void calc_xy_grid2( double u_pix, double v_pix, double delta_u, double delta_v, int n_pixels, int& x_grid, int& y_grid, int im_sign, int is_odd_x, int is_odd_y );
   
   void makeKernels( int _kernelSize, int _overSamplingFactor=1023 );
   
   static void makeKaiserBesselKernel( std::vector<double> &kernel, double alpha, size_t overSamplingFactor, bool withSinc);
   
   //-------------------------------------------------------------------------------------------------------------
   // Anti-aliasing and kernel convolution member variables. Currently taken from WSCLEAN code - to be re-implemented 
   //-------------------------------------------------------------------------------------------------------------
   std::vector<double> _1dKernel;
   std::vector<std::vector<double>> _griddingKernels;
   GridMode _gridMode {GridMode::KaiserBesselKernel};
   int m_nConvolvingKernelSize {0}; // size of the gridding kernel (used to be a parameter in m_ImagerParameters but there is no such object anymore
                                       // I have no clue where parameter are now 
   // END OF functions and members for convolution kernel
};

#endif 
