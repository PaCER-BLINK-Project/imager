#ifndef _PACER_IMAGER_H__
#define _PACER_IMAGER_H__

#include "apply_calibration.h"
#include "pacer_imager_parameters.h"
#include "observation_metadata.h"
#include "pacer_imager_defs.h"
#include <bg_fits.h>

#include <string>
using namespace std;

// FFTW, math etc :
// #include <fftw3.h>

// AstroIO for Visibilities class :
#include <astroio.hpp>
// #define VISIBILITY_TYPE double // visibility type currently double but may be float 
// forward declaration of class  Visibilities, the header file is included in pacer_imager.cpp
// class  Visibilities;

class CPacerImager
{
// simple functions required for some internal calculations : date/time/filename etc :
protected : 
   // memory buffers for internal processing 
   // FFTW / cuFFT memory :
//   fftw_complex* m_in_buffer;
   void*         m_in_buffer; // using void* and casting to (fftw_complex*) later in .cpp file to avoid including fftw3.h which conflicts with cufft.h (cufftComplex !!!)
   int           m_in_size;
//   fftw_complex* m_out_buffer;
   void*         m_out_buffer; // using void* and casting to (fftw_complex*) later in .cpp file to avoid including fftw3.h which conflicts with cufft.h (cufftComplex !!!)
   int           m_out_size;
 

protected :
   // FFT shift to convert from DC in bin 0 of the FFTW output array to DC in the center bin :
   void fft_shift( CBgFits& dirty_image, CBgFits& out_image );

   // FFT unshift converts from DC in the center bin to DC in bin 0 (as expected input to complex FFTW)
   void fft_unshift( CBgFits& dirty_image, CBgFits& out_image );

   // check if image size is as required and re-alloc if not:
   // Return : 
   //    false - if image size was ok and nothing was required
   //    true  - if image size was changed
   bool CheckSize( CBgFits& image, int sizeX, int sizeY );
   
   // allocate memory if not yet done (NULL pointers) or check image size and re-size if needed :
   bool AllocOutPutImages( int sizeX, int sizeY );   
   bool AllocGriddedVis( int sizeX, int sizeY );
   
   // clean local variables allocated locally 
   void CleanLocalAllocations();
   
   // TEMPORARY FOR TESTING - TO BE REMOVED LATER:
   void ConvertXCorr2Fits( Visibilities& xcorr, CBgFits& vis_re, CBgFits& vis_im, int time_step, int fine_channel, const char* szBaseFitsName );
   
   // TODO : required to convert from FITS to xcorr Visibility structure and avoid duplication of code 
   // WARNING : currently only converts single timestep and frequency channel :
   // It also allocated memory inside, which needs to be released later !
   // Visibilities* ConvertFits2XCorr( CBgFits& vis_re, CBgFits& vis_im );

   // Internal functions which are not exposed to the public usage :
   // run_imager with all the parameters, but the one exposed to the public has much fewer parameters and the rest are taken from the m_ImagerParameters object.
   bool run_imager( const char* basename, const char* szPostfix,
                    double frequency_mhz,
                    int    n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,        // minimum UV 
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* weighting="",   // weighting : U for uniform (others not implemented)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="", // gridded visibilities can be provided externally                    
                    const char* szBaseOutFitsName=NULL
                  );


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
   
   // calibration solutions for a single frequency channel (assuming one object of CPacerImager is per frequency channel)
   CCalSols m_CalibrationSolutions;
   
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
   
   // gridded visibilities :
   CBgFits* m_uv_grid_counter;
   CBgFits* m_uv_grid_real;
   CBgFits* m_uv_grid_imag;
   
   // values calculated for the current image :
   double m_PixscaleAtZenith;
   
   //
   // Resulting sky images :
   // WARNING : be careful with using this object as global !
   CBgFits* m_pSkyImageReal;
   CBgFits* m_pSkyImageImag;
   CBgFits* m_pSkyImageRealTmp;
   CBgFits* m_pSkyImageImagTmp;
   bool     m_bLocalAllocation;

   //-------------------------------------------------------------------------------------------------------------
   // image counter : counts how many sky images have been already created
   //-------------------------------------------------------------------------------------------------------------
   int      m_SkyImageCounter;
   

   //-------------------------------------------------------------------------------------------------------------
   // Anti-aliasing and kernel convolution member variables. Currently taken from WSCLEAN code - to be re-implemented 
   //-------------------------------------------------------------------------------------------------------------
   static std::vector<double> _1dKernel;
   static std::vector<std::vector<double>> _griddingKernels;
   static GridMode _gridMode;

   


   CPacerImager();
   ~CPacerImager();
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Initialisation of data structures in the object, which are required for imaging 
   //-----------------------------------------------------------------------------------------------------------------------------   
   void Initialise(); // implement initialisation of object here, read antenna positions, calculate UVW if constant etc 
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // update calibration solutions file when processing multiple channels (in BLINK_PIPELINE) and different files are required for different channels
   //-----------------------------------------------------------------------------------------------------------------------------
   void UpdateCalSolFile( const char* calsol_filename );
 
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
   bool CalculateUVW( bool bForce=false, bool bInitialise=true );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Re-calculate parameters based on new UVW :
   //-----------------------------------------------------------------------------------------------------------------------------
   bool UpdateParameters();
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Reads UVW from FITS files with files names basename + szPostfix + U/V/W.fits
   // or calculates using antenna positions
   //-----------------------------------------------------------------------------------------------------------------------------
   bool ReadOrCalcUVW( const char* basename, const char* szPostfix );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // Set external buffers for output images (Real/Imag) 
   //-----------------------------------------------------------------------------------------------------------------------------
   void SetOutputImagesExternal( CBgFits* pSkyImageRealExt, CBgFits* pSkyImageImagExt );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 1st version producing a dirty image (tested on both MWA and SKA-Low).
   // TODO : Test cases can be found in PaCER documentation 
   //-----------------------------------------------------------------------------------------------------------------------------
   void dirty_image( CBgFits& uv_grid_real_param, CBgFits& uv_grid_imag_param, CBgFits& uv_grid_counter, 
                     double delta_u, double delta_v,
                     bool bSaveIntermediate=false, const char* szBaseOutFitsName=NULL, bool bSaveImaginary=true, bool bFFTUnShift=true );

   //-----------------------------------------------------------------------------------------------------------------------------
   // read input data correlation matrix and UVW from FITS files (this is the first type of input )
   // TODO : add other input formats (e.g. binary data)
   //-----------------------------------------------------------------------------------------------------------------------------
   bool read_corr_matrix( const char* basename, CBgFits& fits_vis_real, CBgFits& fits_vis_imag, const char* szPostfix );

   //-----------------------------------------------------------------------------------------------------------------------------
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
   void gridding( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // GOAL : gridding correction - divide by Fourier Transform of the convolution kernel. See s98book.pdf (Taylor, Carilli et al., 1998)
   //-----------------------------------------------------------------------------------------------------------------------------
   void gridding_correction( CBgFits& image, double delta_u, double delta_v );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // GOAL : convolve with a specified 
   //-----------------------------------------------------------------------------------------------------------------------------
   void convolve_with_kernel( CBgFits& uv_grid_real, CBgFits& uv_grid_imag, int x_grid, int y_grid, double re, double im, double u, double v, double delta_u, double delta_v, int n_pixels, double im_sign, int oversampling=10 );

   //-----------------------------------------------------------------------------------------------------------------------------
   // Calculate coordintes withing in the UV grid array 
   //-----------------------------------------------------------------------------------------------------------------------------   
   void calc_xy_grid( double u, double v, double delta_u, double delta_v, int n_pixels, int& x_grid, int& y_grid, int im_sign, int is_odd_x, int is_odd_y );
   
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
   void gridding_fast( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );

   // same as above 
   // TODO : only use the new one gridding_fast( Visibilities& xcorr, ... ) and use ConvertFits2XCorr in gridding_fast( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, ... ) to
   //        call this one :
   void gridding_fast( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );

   void gridding_fast( void* visibilities_param,  // cast to fftw_complex
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
                  void* gridded_visibilities_param, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="" // weighting : U for uniform (others not implemented)
                );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // 
   // GOAL   : function calling both gridding_fast and imaging. This function is virtual and will be overwritten in HIP version to do GPU specific code
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
   virtual void gridding_imaging( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,
//                  CBgFits& uv_grid_real, CBgFits& uv_grid_imag, CBgFits& uv_grid_counter, double delta_u, double delta_v, 
                  double delta_u, double delta_v,
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL,
                  bool do_gridding=true,                  
                  bool do_dirty_image=true,
                  const char* in_fits_file_uv_re = "", // gridded visibilities can be provided externally
                  const char* in_fits_file_uv_im = "", // gridded visibilities can be provided externally
                  bool bSaveIntermediate=false, bool bSaveImaginary=true
                );

   // same as above but using AstroIO Visibility as input
   // TODO : at some point this one should stay and the other one should be a wrapper using ConvertFits2XCorr and calling this one:
   virtual void gridding_imaging( Visibilities& xcorr, 
                  int time_step, 
                  int fine_channel,
                  CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w,                   
                  double delta_u, double delta_v,
                  double frequency_mhz,
                  int    n_pixels,
                  double min_uv=-1000,    // minimum UV 
                  const char* weighting="", // weighting : U for uniform (others not implemented)
                  const char* szBaseOutFitsName=NULL,
                  bool do_gridding=true,                  
                  bool do_dirty_image=true,
                  const char* in_fits_file_uv_re = "", // gridded visibilities can be provided externally
                  const char* in_fits_file_uv_im = "", // gridded visibilities can be provided externally
                  bool bSaveIntermediate=false, bool bSaveImaginary=true
                );



   //-----------------------------------------------------------------------------------------------------------------------------
   // Executes imager:
   // INPUT  : 
   //          fits_vis_real, fits_vis_imag : visibilities (REAL and IMAG 2D arrays as FITS class) 
   //          fits_vis_u, fits_vis_v, fits_vis_w : UVW (real values baselines in units of wavelength - see TMS)     
   //-----------------------------------------------------------------------------------------------------------------------------
   bool run_imager( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                    double frequency_mhz,
                    int    n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,        // minimum UV 
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* weighting="",   // weighting : U for uniform (others not implemented)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="",  // gridded visibilities can be provided externally
                    const char* szBaseOutFitsName=NULL
                  );

   // same as above, but uses Visibility from the AstroIO library :
   bool run_imager( Visibilities& xcorr, 
                    int time_step, 
                    int fine_channel,
                    CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                    double frequency_mhz,
                    int    n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,        // minimum UV 
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* weighting="",   // weighting : U for uniform (others not implemented)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="",  // gridded visibilities can be provided externally
                    const char* szBaseOutFitsName=NULL
                  );

   //-----------------------------------------------------------------------------------------------------------------------------
   // Wrapper to run_imager :
   // Reads FITS files and executes overloaded function run_imager ( as above )
   //-----------------------------------------------------------------------------------------------------------------------------
   /*bool run_imager( const char* basename, const char* szPostfix,
                    double frequency_mhz,
                    int    n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,        // minimum UV 
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* weighting="",   // weighting : U for uniform (others not implemented)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="", // gridded visibilities can be provided externally                    
                    const char* szBaseOutFitsName=NULL
                  );*/

   bool run_imager( const char* basename, const char* szPostfix,
                    bool   do_gridding=true,    // excute gridding  (?)
                    bool   do_dirty_image=true, // form dirty image (?)
                    const char* in_fits_file_uv_re="", // gridded visibilities can be provided externally
                    const char* in_fits_file_uv_im="", // gridded visibilities can be provided externally                    
                    const char* szBaseOutFitsName=NULL
                  );

                  
   bool run_imager( float* data_real, 
                    float* data_imag,
                    int n_ant, 
                    int n_pol,
                    double frequency_mhz, 
                    int n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,      // minimum UV
                    bool do_gridding=true,    // excute gridding  (?)
                    bool do_dirty_image=true, // form dirty image (?)
                    const char* weighting="", // weighting : U for uniform (others not implemented)
                    const char* szBaseOutFitsName=NULL,
                    bool bCornerTurnRequired=true // changes indexing of data "corner-turn" from xGPU structure to continues (FITS image-like)
                  );

   // overloaded function using astroio class Visibilities as an input:
   bool run_imager( Visibilities& xcorr, 
                    int time_step, 
                    int fine_channel,
                    double frequency_mhz, 
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
   // Apply calibration solutions :
   //-----------------------------------------------------------------------------------------------------------------------------
   bool ApplySolutions( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz, CCalSols& calsol, const char* szPol="X" );
   
   bool ApplySolutions( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz, const char* szPol="X" );
   
   bool ApplySolutions(  Visibilities& xcorr, double frequency_mhz, CCalSols& calsol, int time_step, int fine_channel, const char* szPol="X" );

   //-----------------------------------------------------------------------------------------------------------------------------
   // Apply phase corrections : geometrical correction, cable correction etc   
   //-----------------------------------------------------------------------------------------------------------------------------
   bool ApplyGeometricCorrections( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz );
   virtual bool ApplyGeometricCorrections( Visibilities& xcorr, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, double frequency_mhz, int time_step, int fine_channel );
   
   bool ApplyCableCorrections( CBgFits& fits_vis_real, CBgFits& fits_vis_imag, double frequency_mhz );
   virtual bool ApplyCableCorrections( Visibilities& xcorr, double frequency_mhz, int time_step, int fine_channel );

   //-----------------------------------------------------------------------------------------------------------------------------
   // Function saving output FITS files: 
   //-----------------------------------------------------------------------------------------------------------------------------
   bool SaveSkyImage( const char* outFitsName, CBgFits* pFits, double unixtime=0.00 );
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // test functions - just to test certain other (sometimes non-public functions)
   //-----------------------------------------------------------------------------------------------------------------------------
   void test();
   
   //-----------------------------------------------------------------------------------------------------------------------------
   // IsGPU() - returns true if CPacerImagerHip object and false here:
   //-----------------------------------------------------------------------------------------------------------------------------
   virtual inline bool IsGPU(){ return false; }
   
   //-------------------------------------------------------------------------------------------------------------
   // Anti-aliasing and kernel convolution member functions. Currently taken from WSCLEAN code - to be re-implemented 
   //-------------------------------------------------------------------------------------------------------------
   static void calc_xy_grid2( double u_pix, double v_pix, double delta_u, double delta_v, int n_pixels, int& x_grid, int& y_grid, int im_sign, int is_odd_x, int is_odd_y );
   
   void makeKernels( int _kernelSize, int _overSamplingFactor=1023 );
   
   static void makeKaiserBesselKernel( std::vector<double> &kernel, double alpha, size_t overSamplingFactor, bool withSinc);
   static double bessel0(double x, double precision);


                  
};

#endif 
