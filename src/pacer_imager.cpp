#include "pacer_imager.h"

#include <bg_fits.h>

#include "pacer_common.h"

// FFTW, math etc :
#include <fftw3.h>
#include <math.h>

// local defines :
#include "pacer_imager_defs.h"

// msfitslib library :
#include <libnova_interface.h>
#include <myfile.h>

#ifdef _PACER_PROFILER_ON_
#include <mydate.h>
#endif

// AstroIO for Visibilities class :
#include <astroio.hpp>
#include <memory_buffer.hpp>
#include "files.hpp"

#include <fstream>
#include <omp.h>

void memdump(char *ptr, size_t nbytes, std::string filename){
    std::ofstream outfile;
    outfile.open(filename, std::ofstream::binary);
    outfile.write(ptr, nbytes);
    if(!outfile){
        std::cerr << "Error while dumping data to " << filename << std::endl;
        exit(1);
    }
    outfile.close();
}


namespace {
    
    void save_fits_file(const std::string filename, float* data, long side_x, long side_y){
        FITS fitsImage;
        FITS::HDU hdu;
        hdu.set_image(data,  side_x, side_y);
        // hdu.add_keyword("TIME", static_cast<long>(obsInfo.startTime), "Unix time (seconds)");
        // hdu.add_keyword("MILLITIM", msElapsed, "Milliseconds since TIME");
        // hdu.add_keyword("INTTIME", integrationTime, "Integration time (s)");
        // hdu.add_keyword("COARSE_CHAN", obsInfo.coarseChannel, "Receiver Coarse Channel Number (only used in offline mode)");
        fitsImage.add_HDU(hdu);
        fitsImage.to_file(filename);
    }

    void compare_xcorr_to_fits_file(Visibilities& xcorr, std::string filename){
        auto vis2 = Visibilities::from_fits_file(filename, xcorr.obsInfo);
        size_t fine_channel {0}, int_time {0};
        size_t n_nans {0};
        size_t total {0};
        for(size_t a1 {0}; a1 < xcorr.obsInfo.nAntennas; a1++){
            for(size_t a2 {0}; a2 < a1; a2++){
                std::complex<float> *p1 = xcorr.at(int_time, fine_channel, a1, a2);
                std::complex<float> *p2 = vis2.at(int_time, fine_channel, a1, a2);
                for(size_t p {0}; p < 4; p++){
                    total++;
                    if(isnan(p1->real()) && isnan(p2->real()) && isnan(p2->imag()) && isnan(p1->imag())){
                        n_nans++;
                        continue;
                    }
                    if(*p1 != *p2){
                        std::cerr << "xcorr differs from " << filename << "!!!!" << std::endl;
                        std::cerr << "[a1 = " << a1 << ", a2 = " << a2 << "] p1 = " << *p1 << ", p2 = " << *p2 << std::endl;
                        exit(1);
                    }
                }
            }
        }
        std::cout << "OKK comparison with " << filename << std::endl;
        std::cout << "Percentage NaNs: " << (static_cast<double>(n_nans) / total * 100.0) << std::endl;
    }


    int calc_fft_shift_marcin(int pos, int side){
        int half = side / 2;
        int is_odd = side % 2;             
        if( pos < half ){
            return half + pos + is_odd;
        }else{
            return pos - half;
        }                    
    }

    int calc_fft_shift_cristian(int pos, int side){
        int is_odd = side % 2;
        return (pos + side/2 + is_odd) % (side);
    }

    int calc_fft_shift(int pos, int side){
        return calc_fft_shift_cristian(pos, side);
    }
}

void Images::to_fits_files(const std::string& directory_path, bool save_as_complex, bool save_imaginary) {
    if(on_gpu()) to_cpu();
    MemoryBuffer<float> img_real(this->image_size(), false, false);
    MemoryBuffer<float> img_imag(this->image_size(), false, false);
    for(size_t interval {0}; interval < this->integration_intervals(); interval++){
        for(size_t fine_channel {0}; fine_channel < this->nFrequencies; fine_channel++){
            std::complex<float> *current_data {this->data() + this->image_size() * this->nFrequencies * interval + fine_channel * this->image_size()}; 
            std::stringstream full_directory;
            full_directory << directory_path << "/" << "start_time_" << obsInfo.startTime << \
                "/" << "int_" << interval << "/coarse_" << obsInfo.coarseChannel << "/fine_ch" << fine_channel;
            std::string full_directory_str = full_directory.str();
            blink::imager::create_directory(full_directory_str);
            if(save_as_complex){
                std::string filename {full_directory_str + "/image.fits"};
                std::complex<float>* p_data = this->at(interval, fine_channel);
                ::save_fits_file(filename, reinterpret_cast<float*>(p_data), this->side_size, this->side_size *2);
            }else{
                for(size_t i {0}; i < this->image_size(); i++){
                    img_real[i] = current_data[i].real();
                    img_imag[i] = current_data[i].imag();
                }
                ::save_fits_file(full_directory_str + "/image_real.fits", img_real.data(), this->side_size, this->side_size);
                ::save_fits_file(full_directory_str + "/image_imag.fits", img_imag.data(), this->side_size, this->side_size);
            }
        }
    }
}


// TEST OPTIONS to compare with MIRIAD image
// see memos : PAWSEY/PaCER/logbook/20220305_pacer_imager_validation.odt, MIRIAD
// natural weighting (sup=0) etc: invert vis=chan_204_20211116T203000_yx.uv
// map=chan_204_20211116T203000_iyx.map imsize=180,180
// beam=chan_204_20211116T203000_iyx.beam  sup=0 options=imaginary stokes=yx
// select='uvrange(0.0,100000)'
bool CPacerImager::m_bCompareToMiriad = false;

// debug level : see pacer_imager_defs.h for SAVE_FILES_NONE etc
int CPacerImager::m_ImagerDebugLevel = IMAGER_INFO_LEVEL; // IMAGER_ALL_MSG_LEVEL;

// level of saving intermediate and test files , see pacer_imager_defs.h for
// defines SAVE_FILES_NONE
int CPacerImager::m_SaveFilesLevel = SAVE_FILES_ALL;

// can also save control files every N-th file
int CPacerImager::m_SaveControlImageEveryNth = -1;

// show final image statistics :
bool CPacerImager::m_bPrintImageStatistics = false; // default disabled to make imaging as fast as possible

void CPacerImager::SetDebugLevel(int debug_level)
{
    CPacerImager::m_ImagerDebugLevel = debug_level;
    gBGPrintfLevel = debug_level;
}

void CPacerImager::SetFileLevel(int filesave_level)
{
    CPacerImager::m_SaveFilesLevel = filesave_level;
}

int CPacerImager::UpdateFlags()
{
    int count_flagged = 0;

    if (m_FlaggedAntennas.size() > 0)
    {
        if (m_MetaData.m_AntennaPositions.size() > 0)
        {
            // flagging antennas in the list :
            for (int i = 0; i < m_FlaggedAntennas.size(); i++)
            {
                int ant_index = m_FlaggedAntennas[i];

                if (ant_index >= 0 && ant_index < m_MetaData.m_AntennaPositions.size())
                {
                    m_MetaData.m_AntennaPositions[ant_index].flag = 1;
                    count_flagged++;
                }
            }
            PRINTF_DEBUG("CPacerImager::SetFlaggedAntennas : Flagged %d antennas in the "
                         "imager object\n",
                         count_flagged);
        }
        else
        {
            PRINTF_DEBUG("CPacerImager::SetFlaggedAntennas : No antennas in object "
                         "m_MetaData.m_AntennaPositions\n");
        }
    }

    return count_flagged;
}

void CPacerImager::SetFlaggedAntennas(vector<int> &flagged_antennas)
{
    m_FlaggedAntennas = flagged_antennas;
    UpdateFlags();
}

CPacerImager::CPacerImager()
    : m_bInitialised(false), m_Baselines(0), m_pSkyImageReal(NULL), m_pSkyImageImag(NULL), m_pSkyImageRealTmp(NULL),
      m_pSkyImageImagTmp(NULL), m_bLocalAllocation(false), m_SkyImageCounter(0), m_in_buffer(NULL), m_out_buffer(NULL),
      m_in_size(-1), m_out_size(-1), m_bIncludeAutos(false), m_uv_grid_counter(NULL), m_uv_grid_real(NULL),
      m_uv_grid_imag(NULL), m_nAntennas(0), u_mean(0.00), u_rms(0.00), u_min(0.00), u_max(0.00), v_mean(0.00),
      v_rms(0.00), v_min(0.00), v_max(0.00), w_mean(0.00), w_rms(0.00), w_min(0.00), w_max(0.00)
{
    m_PixscaleAtZenith = 0.70312500; // deg for ch=204 (159.375 MHz) EDA2
}

CPacerImager::~CPacerImager()
{

}

void CPacerImager::CleanLocalAllocations()
{
    if (m_bLocalAllocation)
    { // only remove when locally allocated, as it can
      // also be passed from outside using
      // SetOutputImagesExternal()
        if (m_pSkyImageReal)
        {
            delete m_pSkyImageReal;
            m_pSkyImageReal = NULL;
        }
        if (m_pSkyImageImag)
        {
            delete m_pSkyImageImag;
            m_pSkyImageImag = NULL;
        }

        m_bLocalAllocation = false;
    }

    // gridded visibilites are always local allocations (cannot be passed from
    // external function)
    if (m_uv_grid_counter)
    {
        delete m_uv_grid_counter;
        m_uv_grid_counter = NULL;
    }
    if (m_uv_grid_real)
    {
        delete m_uv_grid_real;
        m_uv_grid_real = NULL;
    }
    if (m_uv_grid_imag)
    {
        delete m_uv_grid_imag;
        m_uv_grid_imag = NULL;
    }

    // these are always only allocated locally :
    if (m_pSkyImageRealTmp)
    {
        delete m_pSkyImageRealTmp;
        m_pSkyImageRealTmp = NULL;
    }
    if (m_pSkyImageImagTmp)
    {
        delete m_pSkyImageImagTmp;
        m_pSkyImageImagTmp = NULL;
    }
}

bool CPacerImager::AllocOutPutImages(int sizeX, int sizeY)
{
    bool bRet = false;
    if (!m_pSkyImageReal)
    {
        m_pSkyImageReal = new CBgFits(sizeX, sizeY);
        m_bLocalAllocation = true;
        bRet = true;
    }
    else
    {
        CheckSize(*m_pSkyImageReal, sizeX, sizeY);
    }

    if (!m_pSkyImageImag)
    {
        m_pSkyImageImag = new CBgFits(sizeX, sizeY);
        m_bLocalAllocation = true;
        bRet = true;
    }
    else
    {
        CheckSize(*m_pSkyImageImag, sizeX, sizeY);
    }

    // temporary buffers, always allocated locally, but flag m_bLocalAllocation is
    // only for output images which can go out of to external calling routines:
    if (!m_pSkyImageRealTmp)
    {
        m_pSkyImageRealTmp = new CBgFits(sizeX, sizeY);
        bRet = true;
    }
    else
    {
        CheckSize(*m_pSkyImageRealTmp, sizeX, sizeY);
    }

    if (!m_pSkyImageImagTmp)
    {
        m_pSkyImageImagTmp = new CBgFits(sizeX, sizeY);
        bRet = true;
    }
    else
    {
        CheckSize(*m_pSkyImageImagTmp, sizeX, sizeY);
    }

    return bRet;
}

bool CPacerImager::AllocGriddedVis(int sizeX, int sizeY)
{
    if (!m_uv_grid_counter)
    {
        m_uv_grid_counter = new CBgFits(sizeX, sizeY);
    }
    else
    {
        CheckSize(*m_uv_grid_counter, sizeX, sizeY);
    }

    if (!m_uv_grid_real)
    {
        m_uv_grid_real = new CBgFits(sizeX, sizeY);
    }
    else
    {
        CheckSize(*m_uv_grid_real, sizeX, sizeY);
    }

    if (!m_uv_grid_imag)
    {
        m_uv_grid_imag = new CBgFits(sizeX, sizeY);
    }
    else
    {
        CheckSize(*m_uv_grid_imag, sizeX, sizeY);
    }

    return true;
}

void CPacerImager::SetOutputImagesExternal(CBgFits *pSkyImageRealExt, CBgFits *pSkyImageImagExt)
{
    CleanLocalAllocations();

    m_pSkyImageReal = pSkyImageRealExt;
    m_pSkyImageImag = pSkyImageImagExt;
    m_bLocalAllocation = false;
}

int CPacerImager::ReadAntennaPositions(bool bConvertToXYZ)
{
    m_nAntennas = m_MetaData.m_AntennaPositions.ReadAntennaPositions(m_ImagerParameters.m_AntennaPositionsFile.c_str(),
                                                                     bConvertToXYZ);
    UpdateFlags(); // if antenna positions are read now for the first time, the
                   // flags need to be updated
    PRINTF_INFO("INFO : read %d antenna positions from file %s\n", m_nAntennas,
                m_ImagerParameters.m_AntennaPositionsFile.c_str());

    return m_nAntennas;
}

void CPacerImager::Initialise(double frequency_hz)
{
    if (!m_bInitialised)
    {
        m_bInitialised = true;

        // test file with antenna positions can be used to overwrite whatever was in
        // .metafits
        if (strlen(m_ImagerParameters.m_AntennaPositionsFile.c_str()) &&
            MyFile::DoesFileExist(m_ImagerParameters.m_AntennaPositionsFile.c_str()))
        {
            bool bConvertToXYZ = false;
            if (!m_ImagerParameters.m_bConstantUVW)
            { // if non-constant UVW -> non zenith phase
              // centered all-sky image
                bConvertToXYZ = true;
            }
            if (m_ImagerParameters.m_bAntennaPositionsXYZ)
            { // text file already has XYZ in WG54
              // system - for example CASA dump
                printf("INFO : antenna positions already in XYZ coordinate system (WG54) "
                       "no need to convert\n");
                bConvertToXYZ = false;
            }

            // read antenna positions and do whatever else is necessary (update flags
            // etc)
            ReadAntennaPositions(bConvertToXYZ);

            if (/*true ||*/ strlen(m_ImagerParameters.m_MetaDataFile.c_str()) == 0)
            { // only calculate UVW here when Metadata is not required
                // initial recalculation of UVW at zenith (no metadata provided ->
                // zenith): WARNING : bool CPacerImager::CalculateUVW() - could be used,
                // but it also calls this function itself which may cause infinite
                // recursise call m_Baselines =
                // m_MetaData.m_AntennaPositions.CalculateUVW( m_U, m_V, m_W,
                // (CPacerImager::m_SaveFilesLevel>=SAVE_FILES_DEBUG),
                // m_ImagerParameters.m_szOutputDirectory.c_str(), m_bIncludeAutos );
                // UpdateParameters();
                CalculateUVW(frequency_hz, true,
                             false); // bForce=true to force call of
                                     // m_MetaData.m_AntennaPositions.CalculateUVW and ,
                                     // bInitialise=false to avoid call to this
                                     // (CPacerImager::Initialise) function in a recursive way !
                PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines (include Autos "
                            "= %d)\n",
                            m_Baselines, m_bIncludeAutos);
            }
            else
            {
                printf("INFO : non-zenith pointing meta data is required to calculate "
                       "UVW\n");
            }
        }
        else
        {
            PRINTF_WARNING("WARNING : antenna position file %s not specified or does not "
                           "exist\n",
                           m_ImagerParameters.m_AntennaPositionsFile.c_str());
        }

        // read all information from metadata
        if (strlen(m_ImagerParameters.m_MetaDataFile.c_str()) &&
            MyFile::DoesFileExist(m_ImagerParameters.m_MetaDataFile.c_str()))
        {
            PRINTF_INFO("INFO : reading meta data from file %s\n", m_ImagerParameters.m_MetaDataFile.c_str());
            if (!m_MetaData.ReadMetaData(m_ImagerParameters.m_MetaDataFile.c_str()))
            {
                PRINTF_ERROR("ERROR : could not read meta data from file %s\n",
                             m_ImagerParameters.m_MetaDataFile.c_str());
            }
        }
    }
}

bool CPacerImager::CheckSize(CBgFits &image, int sizeX, int sizeY)
{
    if (image.GetXSize() != sizeX || image.GetYSize() != sizeY)
    {
        image.Realloc(sizeX, sizeY);
        PRINTF_INFO("DEBUG : change of image size to (%d,%d) was required\n", sizeX, sizeY);

        return true;
    }

    // if image size was ok and nothing was required
    return false;
}

void fft_shift(std::complex<float>* image, size_t image_x_side, size_t image_y_side){
    
    for (size_t y = 0; y < image_y_side; y++){
        for (size_t x = 0; x < image_x_side/2; x++){
            size_t src =  y * image_x_side + x;
            size_t dst_col = ::calc_fft_shift(x, image_x_side);
            size_t dst = y * image_x_side + dst_col;
            std::swap(image[src], image[dst]);
        }
    }
    for (size_t y = 0; y < image_y_side/2; y++){
        for (size_t x = 0; x < image_x_side; x++){
            size_t src =  y * image_x_side + x;
            size_t dst_row = ::calc_fft_shift(y, image_y_side);
            size_t dst = dst_row * image_x_side + x;
            std::swap(image[src], image[dst]);
        }
    }
}



bool CPacerImager::SaveSkyImage(const char *outFitsName, CBgFits *pFits, double unixtime /*=0.00*/)
{
    if (!pFits)
    {
        PRINTF_ERROR("ERROR in code SaveSkyImage, pFits pointer not set\n");
        return false;
    }

    PRINTF_INFO("INFO : saving image %s\n", outFitsName);
    pFits->SetFileName(outFitsName);

    // fill FITS header :
    // TODO :
    pFits->SetKeyword("TELESCOP", "EDA2");

    // scripts/add_fits_header.py
    // TODO - use properly calculated values :
    // double pixscale = m_PixscaleAtZenith/3; // ??? why divided by 3 seems to be
    // best ???
    //   double pixscale =
    //   m_ImagerParameters.m_ImageFOV_degrees/pFits->GetXSize();
    double pixscale = m_ImagerParameters.m_PixsizeInRadians * (180.00 / M_PI);
    if (m_bCompareToMiriad && false)
    { // 2023-12-17 temporary disabled as WCS is not correctly saved
      // then !!! see
      // 20231215_repeat_processing_of_small_part_of_20230709.odt
        PRINTF_WARNING("WARNING (CPacerImager::SaveSkyImage) : MIRIAD-like option -> changing "
                       "pixscale  %.8f -> %.8f\n",
                       pixscale, m_PixscaleAtZenith);
        pixscale = m_PixscaleAtZenith;
    }

    PRINTF_DEBUG("DEBUG : m_PixscaleAtZenith = %.6f [deg] -> pixscale = %.6f [deg] , FoV "
                 "= %.4f [deg], ImageSize = %ld x %ld\n",
                 m_PixscaleAtZenith, pixscale, m_ImagerParameters.m_ImageFOV_degrees, pFits->GetXSize(),
                 pFits->GetYSize());
    // azh2radec 1581220006 mwa 0 90
    // (RA,DEC) = ( 312.07545047 , -26.70331900 )
    // 20.80503003133333333333
    // double ra_deg = 312.07545047; // = 20.80503003133333333333 hours ; //
    // was 2.13673600000E+02; double dec_deg = -2.67033000000E+01;

    // libnova_interface.h
    // void azh2radec( double az, double alt, time_t unix_time, double
    // geo_long_deg, double geo_lat_deg, double& out_ra, double& out_dec ); for
    // all-sky images at zenith :
    double ra_deg = -1000, dec_deg = -1000;
    if (m_ImagerParameters.m_bConstantUVW)
    {
        PRINTF_DEBUG("DEBUG : calculating RA,DEC from (AZ,EL) = (0,90) [deg] for unixtime = "
                     "%d, geo location (%.4f,%.4f)\n",
                     int(unixtime), m_MetaData.geo_long, m_MetaData.geo_lat);
        ::azh2radec(0.00, 90.00, unixtime, m_MetaData.geo_long, m_MetaData.geo_lat, ra_deg, dec_deg);
    }
    else
    {
        // use RA DEC from metafits or something like this
        ra_deg = m_MetaData.raHrs * 15.00;
        dec_deg = m_MetaData.decDegs;
        PRINTF_DEBUG("DEBUG : using RA,DEC = (%.8f,%.8f) [deg] from metafits\n", ra_deg, dec_deg);
    }

    int crpix1 = int(pFits->GetXSize() / 2) + 1;
    pFits->SetKeyword("CTYPE1", "RA---SIN");
    pFits->SetKeyword("CRPIX1", crpix1);
    pFits->SetKeywordFloat("CDELT1", -pixscale); // WARNING / TODO : this may be related to image
                                                 // flip and can differ for different data -> really
                                                 // need to sort this out ASAP !!!
    pFits->SetKeywordFloat("CRVAL1", ra_deg);    // RA of the centre
    pFits->SetKeyword("CUNIT1", "deg    ");

    int crpix2 = int(pFits->GetYSize() / 2) + 1;
    pFits->SetKeyword("CTYPE2", "DEC--SIN");
    pFits->SetKeyword("CRPIX2", crpix2);
    pFits->SetKeywordFloat("CDELT2", pixscale);
    pFits->SetKeywordFloat("CRVAL2", dec_deg); // RA of the centre
    pFits->SetKeyword("CUNIT2", "deg    ");

    // add UTC - also required for scripts pix2sky.py to work ok !
    struct tm gmtime_tm;
    double fraction_sec = 0.00;
    time_t unix_time_time_t = (time_t)unixtime;
    fraction_sec = unixtime - double(unix_time_time_t);
    if (gmtime_r(&unix_time_time_t, &gmtime_tm))
    {
        char tempstring[64];
        // 2023-06-01T10:42:50.1
        sprintf(tempstring, "%04u%02u%02u_%02u%02u%02u.%03d", gmtime_tm.tm_year + 1900, (gmtime_tm.tm_mon + 1),
                gmtime_tm.tm_mday, gmtime_tm.tm_hour, gmtime_tm.tm_min, gmtime_tm.tm_sec, int(fraction_sec * 1000.00));

        pFits->SetKeyword("DATE-OBS", tempstring);
    }
    else
    {
        printf("ERROR : could not convert unixtime = %.4f to UTC using gmtime_r\n", unixtime);
    }

    pFits->WriteFits(outFitsName);

    return true;
}

// Based on example :
// https://github.com/AccelerateHS/accelerate-examples/blob/master/examples/fft/src-fftw/FFTW.c
void CPacerImager::dirty_image(MemoryBuffer<std::complex<float>>& grids_buffer, MemoryBuffer<float>& grids_counters_buffer,
    int grid_side, int n_integration_intervals, int n_frequencies, MemoryBuffer<std::complex<float>>& images_buffer) {

    // TODO CRISTIAN: add check for overflows!
    int width = grid_side;
    int height = grid_side;
    int grid_size = grid_side * grid_side;
    int n_images = n_frequencies * n_integration_intervals;
    const int n[2] {height, width};

    auto tstart = std::chrono::steady_clock::now();
    fftwf_init_threads();
    std::cout << "dirty_image: n threads used = " << omp_get_max_threads() << std::endl;
    fftwf_plan_with_nthreads(omp_get_max_threads());
    fftwf_plan pFwd = fftwf_plan_many_dft(2, n, n_images, reinterpret_cast<fftwf_complex*>(grids_buffer.data()), NULL,
        1, grid_size, reinterpret_cast<fftwf_complex*>(images_buffer.data()), NULL, 1, grid_size, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(pFwd);
    fftwf_destroy_plan(pFwd);
    fftwf_cleanup_threads();
    
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t time_step = 0; time_step < n_integration_intervals; time_step++)
    {
        for (size_t fine_channel = 0; fine_channel < n_frequencies; fine_channel++)
        {
            std::complex<float>* current_grid = grids_buffer.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
            std::complex<float>* current_image = images_buffer.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
            
            float* current_counter = grids_counters_buffer.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
           
            // WARNING : this is image in l,m = cos(alpha), cos(beta) coordinates and
            // still needs to go to SKY COORDINATES !!!

            //   double fnorm = 1.00/sqrt(size); //normalisation see :
            //   /home/msok/Desktop/PAWSEY/PaCER/logbook/20220119_testing_new_versions_dirty_image_polishing.odt
            double counter_sum {0.0};
            for(size_t i {0}; i < grid_size; i++) counter_sum += current_counter[i];
            
            double fnorm = 1.00 / counter_sum; // see RTS :
                                            // /home/msok/mwa_software/RTS_128t/src/newgridder.cu
                                            // SumVisibilityWeights and gridKernel.c:650 also
                                            // read TMS (Thomson, Moran, Swenson) about this
            PRINTF_DEBUG("DEBUG : size = %lu (%d x %d), fnorm = %e (counter sum = %.8f)\n", grid_size, width, height, fnorm,
                        counter_sum);
            for (size_t i = 0; i < grid_size; i++) current_image[i] *= fnorm;

            // TODO: CRISTIAN: is this needed?
            // if(fine_channel == 0)
            // memdump((char *) current_image, grid_size * sizeof(std::complex<double>), "image_before_shift.bin");
            fft_shift(current_image, grid_side, grid_side);
            // if(fine_channel == 0)
            // memdump((char *) current_image, grid_size * sizeof(std::complex<double>), "image_after_shift.bin");
        }
    }

    auto tstop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(tstop - tstart).count();
    std::cout << "Imaging took " << duration << "s" << std::endl;
    // TODO : re-grid to SKY COORDINATES !!!
    // convert cos(alpha) to alpha - see notes !!!
    // how to do it ???

}

bool CPacerImager::CalculateUVW(double frequency_hz, bool bForce /*=false*/, bool bInitialise /*=true*/)
{
    if (bInitialise)
    { // is not required in one case of call from ::Initialise
      // itself -> to avoid recursive call
        Initialise(frequency_hz);
    }

    bool bRecalculationRequired = false;

    if (m_Baselines <= 0 || m_U.GetXSize() <= 0 || m_V.GetXSize() <= 0 || m_W.GetXSize() <= 0 ||
        !m_ImagerParameters.m_bConstantUVW || bForce)
    {
        bRecalculationRequired = true;
    }

    if (bRecalculationRequired)
    {
        PRINTF_DEBUG("DEBUG : recalculation of UVW is required\n");

        m_Baselines = m_MetaData.m_AntennaPositions.CalculateUVW(
            m_U, m_V, m_W, (CPacerImager::m_SaveFilesLevel >= SAVE_FILES_DEBUG),
            m_ImagerParameters.m_szOutputDirectory.c_str(), m_bIncludeAutos);
        PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines\n", m_Baselines);

        UpdateParameters(frequency_hz); // update parameters after change in UVW
    }

    return (m_Baselines > 0);
}

bool CPacerImager::UpdateParameters(double frequency_hz)
{
    PRINTF_DEBUG("DEBUG : updating parameters (pixscale etc.) based on new UVW\n");
    // U :
    m_U.GetStat(u_mean, u_rms, u_min, u_max);

    // V :
    m_V.GetStat(v_mean, v_rms, v_min, v_max);

    // W :
    m_W.GetStat(w_mean, w_rms, w_min, w_max);

    // Bacause we are also including conjugates at (-u,-v) UV point in gridding
    // u_min = -u_max and v_min = -v_max : was -35 / +35
    u_min = -u_max;
    //  u_max = +35;
    v_min = -v_max;
    //  v_max = +35;

    // recalculate PIXSCALE at zenith :
    double wavelength_m = VEL_LIGHT / frequency_hz;
    double u_max_lambda = u_max / wavelength_m;
    double pixscale_radians = 1.00 / (2.00 * u_max_lambda);
    m_PixscaleAtZenith = pixscale_radians * (180.00 / M_PI);
    PRINTF_INFO("INFO : UVW updated UV range (%.6f,%.6f) - (%.6f,%.6f), pixscale at "
                "zenith calculated as %.6f [deg] ( = %.3f [arcmin] ) for frequency_hz = "
                "%.8f Hz\n",
                u_min / wavelength_m, v_min / wavelength_m, u_max / wavelength_m, v_max / wavelength_m,
                m_PixscaleAtZenith, m_PixscaleAtZenith * 60.00, frequency_hz);

    return true;
}

inline int wrap_index(int i, int side){
    if(i >= 0) return i % side;
    else return (side + i);
}

void CPacerImager::gridding_fast(Visibilities &xcorr, int time_step, int fine_channel, CBgFits &fits_vis_u,
                                 CBgFits &fits_vis_v, CBgFits &fits_vis_w,
                                 MemoryBuffer<std::complex<float>> &grids_buffer,
                                 MemoryBuffer<float> &grids_counters_buffer, double delta_u, double delta_v, int n_pixels,
                                 double min_uv /*=-1000*/,
                                 const char *weighting /*="" weighting : U for uniform (others not implemented) */
)
{
    PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f\n", min_uv);

    bool bStatisticsCalculated = false;
    if (m_bPrintImageStatistics)
    { // TODO : ? may require a separate flag in the
      // future, for now just using a single
      // Statistics switch ON/OFF flag
        fits_vis_u.GetStat(u_mean, u_rms, u_min, u_max);
        fits_vis_v.GetStat(v_mean, v_rms, v_min, v_max);
        fits_vis_w.GetStat(w_mean, w_rms, w_min, w_max);

        // Bacause we are also including conjugates at (-u,-v) UV point in gridding
        // u_min = -u_max and v_min = -v_max : was -35 / +35
        u_min = -u_max;
        //  u_max = +35;
        v_min = -v_max;
        //  v_max = +35;

        bStatisticsCalculated = true;
    }
    memset(grids_buffer.data(), 0, grids_buffer.size() * sizeof(std::complex<float>));
    memset(grids_counters_buffer.data(), 0, grids_counters_buffer.size() * sizeof(float));
    int grid_size = n_pixels * n_pixels;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++)
    {
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++)
        {
            std::complex<float>* current_grid = grids_buffer.data() + time_step * xcorr.nFrequencies * grid_size + fine_channel * grid_size;
            float* current_counter = grids_counters_buffer.data() + time_step * xcorr.nFrequencies * grid_size + fine_channel * grid_size;
            // calculate using CASA formula from image_tile_auto.py :
            // synthesized_beam=(lambda_m/max_baseline)*(180.00/math.pi)*60.00 # in
            // arcmin lower=synthesized_beam/5.00 higher=synthesized_beam/3.00
            // cellside_float=(lower+higher)*0.5
            // NEW :
            //   double alpha = 4.00/15.00; // from (1/5+1/3)*0.5 = 4/15
            //   double wrong_factor = 1.00; // 4.00 factor for due to U range -35 to
            //   +35 instead of ~-17.5/2 to +17.5/2 (factor 4 ) double delta_u =
            //   wrong_factor*alpha*(u_max-u_min)/(n_pixels); // factor for due to U
            //   range -35 to +35 instead of ~-17.5/2 to +17.5/2 (factor 4 ) double
            //   delta_v = wrong_factor*alpha*(v_max-v_min)/(n_pixels); int
            //   freq_channel = 204; double frequency_mhz =
            //   freq_channel*(400.00/512.00); if( gFrequencyMHz > 0 ){
            //      frequency_mhz = gFrequencyMHz;
            //   }
            double frequency_hz = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
            double wavelength_m = VEL_LIGHT / frequency_hz;

            // is it ok to chose the UV plane center based on this:
            //  double u_center = (u_min + u_max)/2.00;
            //  double v_center = (v_min + v_max)/2.00;

            int n_ant = xcorr.obsInfo.nAntennas;
            int added = 0, high_value = 0;
            for (int ant1 = 0; ant1 < n_ant; ant1++)
            {
                for (int ant2 = 0; ant2 <= ant1; ant2++)
                {
                    if (ant1 > ant2 || (m_bIncludeAutos && ant1 == ant2))
                    { // was ant1 > ant2
                        if (m_MetaData.m_AntennaPositions.size() > 0)
                        {
                            if (m_MetaData.m_AntennaPositions[ant1].flag > 0 ||
                                m_MetaData.m_AntennaPositions[ant2].flag > 0)
                            {
                                // skip flagged antennas
                                // printf("Fast flagging used\n");
                                continue;
                            }
                        }
                        else
                        {
                            if (m_FlaggedAntennas.size() > 0)
                            {
                                // WARNING : this is much less efficient so better to have
                                // antenna positions and check flags there
                                if (find_value(m_FlaggedAntennas, ant1) >= 0 ||
                                    find_value(m_FlaggedAntennas, ant2) >= 0)
                                {
                                    continue;
                                }
                            }
                        }

                        std::complex<float> *vis = xcorr.at(time_step, fine_channel, ant1,
                                                            ant2); // was ant1, ant2 , but ant2,ant1 does not fix
                                                                   // the orientation of the final image either ...
                        double re = vis->real();                   // fits_vis_real.getXY(ant1,ant2);
                        double im = vis->imag();                   // fits_vis_imag.getXY(ant1,ant2);

                        if (!isnan(re) && !isnan(im))
                        {
                            if (fabs(re) < MAX_VIS && fabs(im) < MAX_VIS)
                            {
                                // TODO convert [m] -> wavelength
                                double u = fits_vis_u.getXY(ant1, ant2) / wavelength_m;

                                // 2022-09-24 : - removed for a test on MWA data
                                // 2022-09-09 - for now sticking to - sign here to have back
                                // comatible test EDA2 data
                                //  but looks like this (-) should not be here at least does not
                                //  match UV coverage from WSCEAN (it is flipped then see :
                                //  /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                                double v = fits_vis_v.getXY(ant1, ant2) /
                                           wavelength_m; // the - sign here fixes the Y flip, but I am
                                                         // not sure why needed ??? check RTS :
                                                         // imagefromuv.c , LM_CopyFromFFT where some
                                                         // interesting re-shuffling happens too !
                                                         // - was for EDA2 data, but + is ok for MWA
                                                         // data - it may be consistent once I start
                                                         // using TMS equation 4.1 consistently for
                                                         // both EDA2 and MWA
                                double w = fits_vis_w.getXY(ant1, ant2) / wavelength_m;
                                double uv_distance = sqrt(u * u + v * v);

                                /* this is in WSCLEAN, but here seems to have no effect ...
                                                 if (w < 0.0 ) { // && !_isComplex
                                                    u = -u;
                                                    v = -v;
                                                    w = -w;
                                                    im = -im;
                                                 }
                                */

                                if (ant1 == ant2)
                                {
                                    PRINTF_DEBUG("Auto-correlation debug2 values %.4f / %.4f , "
                                                 "uv_distance = %.8f vs. min_uv = %.8f (u = %.8f , v = "
                                                 "%.8f , wavelength_m = %.8f [m])\n",
                                                 re, im, uv_distance, min_uv, u, v, wavelength_m);
                                }

                                if (uv_distance > min_uv)
                                { // check if larger than minimum UV distance
                                    int u_pix = static_cast<int>(round(u / delta_u));
                                    int v_pix = static_cast<int>(round(v / delta_v));


                                    int u_index = wrap_index(u_pix, n_pixels); //+ n_pixels / 2;
                                    int v_index = wrap_index(v_pix, n_pixels); //+ n_pixels / 2;
                                    if(u_index < 0) printf("MODULO ERROR: u_index = %d , u_pix = %d , %d \n", u_index, u_pix, n_pixels);
                                    if(v_index < 0) printf("MODULO ERROR: v_index = %d , v_pix = %d , %d \n", v_index, v_pix, n_pixels);


                                    // now fft shift
                                    //u_index = ::calc_fft_shift(u_index, n_pixels);
                                    //v_index = ::calc_fft_shift(v_index ,n_pixels);

                                    // Using CELL averaging method or setXY ?
                                    current_grid[v_index * n_pixels + u_index].real(current_grid[v_index * n_pixels + u_index].real() + re);
                                    current_grid[v_index * n_pixels + u_index].imag(current_grid[v_index * n_pixels + u_index].imag() + im);
                                    current_counter[v_index * n_pixels + u_index] += 1;

                                    // add conjugates :
                                    u_index = wrap_index(-u_pix, n_pixels);
                                    v_index = wrap_index(-v_pix, n_pixels);
                                    // now fft shift
                                    //u_index = calc_fft_shift(u_index, n_pixels);
                                    //v_index = calc_fft_shift(v_index ,n_pixels);
                                    
                                    current_grid[v_index * n_pixels + u_index].real(current_grid[v_index * n_pixels + u_index].real() + re);
                                    current_grid[v_index * n_pixels + u_index].imag(current_grid[v_index * n_pixels + u_index].imag() - im);
                                    current_counter[v_index * n_pixels + u_index] += 1;
                                }
                            }
                            else
                            {
                                PRINTF_DEBUG("DEBUG : visibility value %e +j%e higher than limit %e -> "
                                             "skipped\n",
                                             re, im, MAX_VIS);
                            }
                        }
                    }
                }
            }
             // This division is in fact UNIFORM weighting !!!! Not CELL-avareging
            // normalisation to make it indeed CELL-averaging :
            if (strcmp(weighting, "U") == 0)
            {
                for (size_t i {0}; i < n_pixels * n_pixels; i++){
                    current_grid[i].real(current_grid[i].real() / current_counter[i]);
                    current_grid[i].imag(current_grid[i].imag() / current_counter[i]);

                }
            }
        }
    }
}

Images CPacerImager::gridding_imaging(Visibilities &xcorr, int time_step, int fine_channel, CBgFits &fits_vis_u,
                                    CBgFits &fits_vis_v, CBgFits &fits_vis_w, double delta_u, double delta_v,
                                    int n_pixels, double min_uv /*=-1000*/, // minimum UV
                                    const char *weighting /*=""*/,          // weighting : U for uniform (others not
                                                                            // implemented)
                                    const char *szBaseOutFitsName /*=NULL*/, bool do_gridding, bool do_dirty_image,
                                    const char *in_fits_file_uv_re,
                                    /*=""*/ // gridded visibilities can be provided externally
                                    const char *in_fits_file_uv_im,
                                    /*=""*/ // gridded visibilities can be provided externally
                                    bool bSaveIntermediate /*=false*/, bool bSaveImaginary /*=true*/
)
{
    printf("DEBUG : gridding_imaging( Visibilities& xcorr ) in pacer_imager.cpp\n");
    // allocates data structures for gridded visibilities:
    if(xcorr.on_gpu()) xcorr.to_cpu();
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/03_after_geo_corrections.fits");
    size_t n_images{xcorr.integration_intervals() * xcorr.nFrequencies};
    size_t buffer_size{n_pixels * n_pixels * n_images};
    MemoryBuffer<float> grids_counters_buffer(buffer_size, false, false);
    MemoryBuffer<std::complex<float>> grids_buffer(buffer_size, false, false);
    // Should be long long int, but keeping float now for compatibility reasons
    if (do_gridding)
    {
        gridding_fast(xcorr, time_step, fine_channel, fits_vis_u, fits_vis_v, fits_vis_w, grids_buffer,
                      grids_counters_buffer, delta_u, delta_v, n_pixels, min_uv, weighting);
    }

    // CBgFits reference_grid, reference_grid_counter;
    // reference_grid.ReadFits("/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data_2/1592584200/133/000/uv_grid_real_8192x8192.fits", 0, 1, 1 );
    // reference_grid_counter.ReadFits("/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data_2/1592584200/133/000/uv_grid_counter_8192x8192.fits", 0, 1, 1 );

    // for(size_t i {0}; i < n_pixels * n_pixels; i++){
    //     if(grids_counters_buffer[i] != reference_grid_counter.getXY(i % n_pixels, i / n_pixels)){
    //         std::cerr << "Error!! Counters are not the same at position " << i << ": " << grids_counters_buffer[i] << " != " << reference_grid_counter.getXY(i % n_pixels, i / n_pixels) << std::endl;
    //         break;
    //     }
    // }
    // for(size_t i {0}; i < n_pixels * n_pixels; i++){
    //     if(std::abs(grids_buffer[i].real() - reference_grid.getXY(i % n_pixels, i / n_pixels)) > 1e-4){
    //         std::cerr << "Error!! Grids are not the same at position " << i << ": " << grids_buffer[i].real() << " != " << reference_grid.getXY(i % n_pixels, i / n_pixels) << std::endl;
    //         exit(1);
    //     }
    // }
    std::cout << "OK!!! GRIDS are fine!" << std::endl;
    // // TODO: Cristian investigate use of single precision fftw
    //     // need this memory allocation just to catch the buffer overflow happening in fft_shift!!
    // // we cannot free this memory as it is corrupted
    // MemoryBuffer<std::complex<double>> images_buffer_double(buffer_size, false, false);
    MemoryBuffer<std::complex<float>> images_buffer_float(buffer_size, false, false);
    if (do_dirty_image)
    {
        // dirty image :
        PRINTF_INFO("PROGRESS : executing dirty image\n");
        dirty_image(grids_buffer, grids_counters_buffer, n_pixels, xcorr.integration_intervals(), xcorr.nFrequencies, images_buffer_float);
    }
    // float *dest {reinterpret_cast<float*>(images_buffer_float.data())};
    // double *src {reinterpret_cast<double*>(images_buffer_double.data())};
    // for(size_t i {0}; i < buffer_size * 2; i++) dest[i] = static_cast<float>(src[i]);
    return {std::move(images_buffer_float), xcorr.obsInfo, xcorr.nIntegrationSteps, xcorr.nAveragedChannels, static_cast<unsigned int>(n_pixels)};
}

bool CPacerImager::ApplyGeometricCorrections(Visibilities &xcorr, CBgFits &fits_vis_w)
{
    if(xcorr.on_gpu()){
        std::cout << "xcorr is on GPU!!! Moving to cpu.." << std::endl;
        xcorr.to_cpu();
    }else{
        std::cout << "xcorr is on CPU.." << std::endl;
    }
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/01_before_geo_corrections.fits");
    int n_ant = xcorr.obsInfo.nAntennas;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++)
    {
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++)
        {
            for (int ant1 = 0; ant1 < n_ant; ant1++)
            {
                for (int ant2 = 0; ant2 <= ant1; ant2++)
                {
                    double w = fits_vis_w.getXY(ant1, ant2); // or ant2,ant1 , was (+1)
                    double angle = - 2.0 * M_PI * w * this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE) /
                                   SPEED_OF_LIGHT; // TODO : + or - here ??? In brute force branch
                                                   // was -

                    double sin_angle, cos_angle;
                    sincos(angle, &sin_angle, &cos_angle);

                    std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2);

                    double re = vis[0].real();
                    double im = vis[0].imag();

                    double re_prim = re * cos_angle - im * sin_angle;
                    double im_prim = im * cos_angle + re * sin_angle;

                    std::complex<double> vis_new(re_prim, im_prim);
                    *vis = vis_new;
                }
            }
        }
    }
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/01_after_geo_corrections.fits");
    return true;
}

bool CPacerImager::ApplyCableCorrections(Visibilities &xcorr)
{
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/02_before_cable_corrections.fits");
    int n_ant = xcorr.obsInfo.nAntennas;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++)
    {
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++)
        {
            for (int ant1 = 0; ant1 < n_ant; ant1++)
            {
                InputMapping &ant1_info = m_MetaData.m_AntennaPositions[ant1];
                for (int ant2 = 0; ant2 <= ant1; ant2++)
                {
                    InputMapping &ant2_info = m_MetaData.m_AntennaPositions[ant2];

                    double cableDeltaLen = ant2_info.cableLenDelta - ant1_info.cableLenDelta; // was ant2 - ant1
                    double angle = -2.0 * M_PI * cableDeltaLen *
                                   this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE) / SPEED_OF_LIGHT;

                    double sin_angle, cos_angle;
                    sincos(angle, &sin_angle, &cos_angle);

                    std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2);

                    double re = vis[0].real();
                    double im = vis[0].imag(); // TODO : why do I need 1 here ???
                    double re_prim = re * cos_angle - im * sin_angle;
                    double im_prim = im * cos_angle + re * sin_angle;

                    std::complex<double> vis_new(re_prim, im_prim);
                    *vis = vis_new;
                }
            }
        }
    }
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/02_after_cable_corrections.fits");
    
    return true;
}

//-----------------------------------------------------------------------------------------------------------------------------
// Wrapper to run_imager as before : run_imager( CBgFits& fits_vis_real,
// CBgFits& fits_vis_imag, CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits&
// fits_vis_w ...) it uses AstroIO class Visibility Reads FITS files and
// executes overloaded function run_imager ( as above )
//-----------------------------------------------------------------------------------------------------------------------------
Images CPacerImager::run_imager(Visibilities &xcorr, int time_step, int fine_channel, CBgFits &fits_vis_u,
                              CBgFits &fits_vis_v, CBgFits &fits_vis_w, // UVW
                              int n_pixels, double FOV_degrees,
                              double min_uv,         // =-1000, minimum UV
                              bool do_gridding,      // =true, excute gridding  (?)
                              bool do_dirty_image,   // =true, form dirty image (?)
                              const char *weighting, // ="",  weighting : U for uniform (others not implemented)
                              const char *in_fits_file_uv_re, // ="", gridded visibilities can be
                                                              // provided externally
                              const char *in_fits_file_uv_im, // ="",  gridded visibilities can be
                                                              // provided externally
                              const char *szBaseOutFitsName   // =NULL
)
{
    //   // based on RTS : UV pixel size as function FOVtoGridsize in
    //   /home/msok/mwa_software/RTS_128t/src/gridder.c double frequency_hz =
    //   frequency_mhz*1e6; double wavelength_m = VEL_LIGHT / frequency_hz;
    double FoV_radians = FOV_degrees * M_PI / 180.;
    // WARNING: it actually cancels out to end up being 1/PI :
    // TODO simplity + 1/(2FoV) !!! should be
    //  double delta_u = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) /
    //  wavelength_m; // in meters (NOT WAVELENGHTS) double delta_v = (
    //  (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in
    //  meters (NOT WAVELENGHTS)
    // TODO :
    double delta_u = 1.00 / (FoV_radians); // should be 2*FoV_radians - see TMS etc
    double delta_v = 1.00 / (FoV_radians); // Rick Perley page 16 :
                                           // /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf

    if (true)
    {
        // see : Rick Perley page 16 :
        // /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
        //       and
        //       /home/msok/Desktop/PAWSEY/PaCER/logbook/20240315_Jishnu_test.odt

        // WARNING : why this is not u_max/wavelength_m - when I use this the image
        // because Frequency dependent. But it does not make sense to have
        //           delta_u in meters and then UV grid in wavelengths ...
        // MWA TEST:
        delta_u =
            2.00 * (u_max) /
            n_pixels; // delta_u = 2.00*(u_max)/n_pixels; // Looks like this is
                      // what it should be NOT u_max/wavelength_m . So delta_u must
                      // be in meters here !!! It may all depend on calculation if
                      // u_index see discussion in
                      // /home/msok/Desktop/PAWSEY/PaCER/logbook/20240320_gridding_delta_u_meters_vs_wavelengths.odt
        delta_v = 2.00 * (v_max) / n_pixels; // and it's not ok because then delta_u is different
                                             // for both of them, which causes exp/shrink with freq

        // automatic calculation of pixel size in radians 1/(2u_max) - see Rick
        // Perley or just Lambda/B_max divide by 2 for oversampling.
        m_ImagerParameters.m_PixsizeInRadians =
            1.00 / (2.00 * u_max); // does this one need to be /wavelength or not ???

        // NEW : based on desired image resolution
        // double delta_theta = (wavelength_m/35.0)/2.00; // 2 times oversampled
        // double delta_theta = ((230.0/300.00)/(2.00*35.00)); // use maximum
        // resolution oversampled by a factor of 2., at 230 MHz -> Lambda ~1.3043m
        // double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
        double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
        // MWA TEST
        //     delta_u = 1.00/(n_pixels*delta_theta);
        //     delta_v = 1.00/(n_pixels*delta_theta);
        PRINTF_DEBUG("delta_u = %.8f (u_max = %.8f), delta_v = %.8f (v_max = %.8f), "
                     "calculated as 1/FoV = 1/(%d pixels * %.5f rad), delta_theta = %.5f "
                     "[deg]\n",
                     delta_u, u_max, delta_v, v_max, n_pixels, delta_theta, delta_theta * (180.00 / M_PI));

        // PRINTF_DEBUG("delta_u = %.8f , delta_v = %.8f , calculated
        // as 2.00*u_max/n_pixels, u_max = %.8f, n_pixels =
        // %d\n",delta_u,delta_v,u_max,n_pixels);
    }


    if (do_gridding || do_dirty_image)
    {
        // virtual function calls gridding and imaging in GPU/HIP version it is
        // overwritten and both gridding and imaging are performed on GPU memory :
        return gridding_imaging(xcorr, time_step, fine_channel, fits_vis_u, fits_vis_v, fits_vis_w, delta_u, delta_v, n_pixels,
                         min_uv, weighting, szBaseOutFitsName, do_gridding, do_dirty_image, in_fits_file_uv_re,
                         in_fits_file_uv_im);
    }
}

void CPacerImager::ConvertXCorr2Fits(Visibilities &xcorr, CBgFits &vis_re, CBgFits &vis_im, int time_step,
                                     int fine_channel, const char *szBaseFitsName)
{
    int n_ant = xcorr.obsInfo.nAntennas;
    int n_corrs = 4; // 4 correlation products : XX XY YX YY
    int n_baselines = n_ant * (n_ant + 1) / 2;
    ObservationInfo &obsInfo = xcorr.obsInfo;
    const size_t matrixSize = n_baselines * obsInfo.nPolarizations * obsInfo.nPolarizations;
    const size_t nIntervals = (obsInfo.nTimesteps); // TODO  + voltages.nIntegrationSteps - 1) /
                                                    // voltages.nIntegrationSteps;
    unsigned int nChannelsToAvg = 1;                // TODO : verify
    const size_t nOutFrequencies = obsInfo.nFrequencies / nChannelsToAvg;
    const size_t nValuesInTimeInterval = matrixSize * nOutFrequencies;
    const size_t outSize = nValuesInTimeInterval * nIntervals;
    unsigned int avgCh = fine_channel / nChannelsToAvg;

    // size_t outIndex {interval * nValuesInTimeInterval + avgCh * matrixSize +
    // baseline * obsInfo.nPolarizations * obsInfo.nPolarizations
    //                             + p1*obsInfo.nPolarizations + p2};

    // TODO !!!

    // assuming single timestep for now :
    // assuming ordering of data as Cristian told me during the meeting :
    // Ant11        | Ant 12       | Ant 13       | ...
    // XX XY YX YY  | XX XY YX YY  | XX XY YX YY  | ...
    int index = 0;

    vis_re.SetNaN();
    vis_im.SetNaN();

    printf("DEBUG : CPacerImager::ConvertXCorr2Fits n_ant = %d\n", n_ant);
    // using at :
    // Complex<float> *at_float(Visibilities& vis, unsigned int interval, unsigned
    // int frequency, unsigned int a1, unsigned a2)
    for (int i = 0; i < n_ant; i++)
    { // loop over ant1
        for (int j = 0; j <= i; j++)
        { // loop over ant2
            // auto& vis = xcorr.data[idx];
            //        std::complex<double>* vis = at( xcorr, time_step, fine_channel,
            //        i, j );
            std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, i, j);

            vis_re.setXY(j, i, float(vis[0].real()));
            vis_im.setXY(j, i, float(vis[0].imag()));

            vis_re.setXY(i, j, float(vis[0].real()));
            vis_im.setXY(i, j, -float(vis[0].imag()));
        }
        index += (n_ant - i);
    }

    char szReFits[1024], szImFits[1024];
    sprintf(szReFits, "%s/%s_re.fits", m_ImagerParameters.m_szOutputDirectory.c_str(), szBaseFitsName);
    sprintf(szImFits, "%s/%s_im.fits", m_ImagerParameters.m_szOutputDirectory.c_str(), szBaseFitsName);
    vis_re.WriteFits(szReFits);
    vis_im.WriteFits(szImFits);
    printf("DEBUG : saved %s and %s\n", szReFits, szImFits);
}

Images CPacerImager::run_imager(Visibilities &xcorr, int time_step, int fine_channel, int n_pixels, double FOV_degrees,
                              double min_uv /*=-1000*/,      // minimum UV
                              bool do_gridding /*=true*/,    // excute gridding  (?)
                              bool do_dirty_image /*=true*/, // form dirty image (?)
                              const char *weighting /*=""*/, // weighting : U for uniform (others not
                                                             // implemented)
                              const char *szBaseOutFitsName /*=NULL*/,
                              bool bCornerTurnRequired /*=true*/ // changes indexing of data
                                                                 // "corner-turn" from xGPU structure to
                                                                 // continues (FITS image-like)
)
{
    // TODO Cristian: time_step, fine_channel will be used in the to select a
    // subset of data to be imaged. ensures initalisation of object structures
    // TODO: this init function must be modified
    double initial_frequency_hz = this->get_frequency_hz(xcorr, fine_channel < 0 ? 0 : fine_channel, COTTER_COMPATIBLE);
    Initialise(initial_frequency_hz);
    int n_ant = xcorr.obsInfo.nAntennas;
    int n_pol = xcorr.obsInfo.nPolarizations;
    // Why not set in constructor?
    m_ImagerParameters.m_ImageFOV_degrees = FOV_degrees; // 2024-03-24 - reasons as above

    if (m_ImagerParameters.m_fUnixTime <= 0.0001)
    {
        m_ImagerParameters.m_fUnixTime = get_dttm_decimal();
        PRINTF_WARNING("Time of the data not specified -> setting current time %.6f\n", m_ImagerParameters.m_fUnixTime);
    }
    xcorr.to_cpu();
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/00_before_any_operation.fits");
    
    // calculate UVW (if required)
    CalculateUVW(initial_frequency_hz);

    if (m_ImagerParameters.m_bApplyGeomCorr)
        ApplyGeometricCorrections(xcorr, m_W);

    if (m_ImagerParameters.m_bApplyCableCorr)
        ApplyCableCorrections(xcorr);

    printf("DEBUG : just before run_imager(time_step=%d, fine_channel=%d )\n", time_step, fine_channel);
    fflush(stdout);
    return run_imager(xcorr, time_step, fine_channel, m_U, m_V, m_W, n_pixels, FOV_degrees, min_uv, do_gridding,
                          do_dirty_image, weighting, "", "", szBaseOutFitsName);

}

double CPacerImager::get_frequency_hz(const Visibilities &vis, int fine_channel, bool cotter_compatible)
{
    // TODO: remove hardcoded values!!
    double fine_ch_bw = 0.04, coarse_ch_bw = 1.28;
    double coarse_channel_central_freq_MHz = vis.obsInfo.coarseChannel * coarse_ch_bw;
    double channel_frequency_MHz;
    if (cotter_compatible)
    {
        // see for example awk commands in
        // experiments/image_mwa_obsid1276619416_allch.sh frequencyMHz is assumed to
        // be center frequency of coarse channel - 0.64 -> lower edge of the MWA
        // coarse channel:
        // TODO : get this information from xcorr structure or metedata, otherwise
        // it will not work for EDA2 etc:
        //        it looks to me that the required fields need to be added there
        //        first
        channel_frequency_MHz = coarse_channel_central_freq_MHz - coarse_ch_bw / 2.00 +
                                fine_ch_bw * fine_channel; // cotter has a bug ? - does not add half of fine
                                                           // channel to calculate channel frequency
    }
    else
    {
        channel_frequency_MHz =
            coarse_channel_central_freq_MHz - coarse_ch_bw / 2.00 + fine_ch_bw * fine_channel + fine_ch_bw / 2.00;
    }
    return channel_frequency_MHz * 1e6;
}