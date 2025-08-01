#include "pacer_imager.h"
#include "pacer_common.h"
#include <bg_fits.h>
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
#include "corrections.h"
#include <fstream>
#include <omp.h>


namespace {
    int calc_fft_shift(int pos, int side){
        int is_odd = side % 2;
        return (pos + side/2 + is_odd) % (side);
    }
}



// debug level : see pacer_imager_defs.h for SAVE_FILES_NONE etc
int CPacerImager::m_ImagerDebugLevel = IMAGER_INFO_LEVEL; // IMAGER_ALL_MSG_LEVEL;

// level of saving intermediate and test files , see pacer_imager_defs.h for
// defines SAVE_FILES_NONE
int CPacerImager::m_SaveFilesLevel = SAVE_FILES_ALL;


// show final image statistics :
bool CPacerImager::m_bPrintImageStatistics = false; // default disabled to make imaging as fast as possible

void CPacerImager::SetDebugLevel(int debug_level)
{
    CPacerImager::m_ImagerDebugLevel = debug_level;
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


CPacerImager::CPacerImager(const std::string metadata_file, const std::vector<int>& flagged_antennas, bool average_images) {
    this->average_images = average_images;
     // read all information from metadata
    if (metadata_file.length() > 0 && MyFile::DoesFileExist(metadata_file.c_str())) {
        PRINTF_INFO("INFO : reading meta data from file %s\n", metadata_file.c_str());
        if( !m_MetaData.ReadMetaData( metadata_file.c_str())){
            PRINTF_ERROR("ERROR : could not read meta data from file %s\n", metadata_file.c_str() );
        }
    }
    m_FlaggedAntennas = flagged_antennas;
    UpdateFlags();
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



// Based on example :
// https://github.com/AccelerateHS/accelerate-examples/blob/master/examples/fft/src-fftw/FFTW.c
void CPacerImager::dirty_image(MemoryBuffer<std::complex<float>>& grids, MemoryBuffer<float>& grids_counters,
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
    fftwf_plan pFwd = fftwf_plan_many_dft(2, n, n_images, reinterpret_cast<fftwf_complex*>(grids.data()), NULL,
        1, grid_size, reinterpret_cast<fftwf_complex*>(images_buffer.data()), NULL, 1, grid_size, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftwf_execute(pFwd);
    fftwf_destroy_plan(pFwd);
    fftwf_cleanup_threads();
    //images_buffer.dump("images_after_fft.bin");
    #pragma omp parallel for collapse(2) schedule(static)
    for (size_t time_step = 0; time_step < n_integration_intervals; time_step++)
    {
        for (size_t fine_channel = 0; fine_channel < n_frequencies; fine_channel++)
        {
            std::complex<float>* current_grid = grids.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
            std::complex<float>* current_image = images_buffer.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
            
            float* current_counter = grids_counters.data() + time_step * n_frequencies * grid_size + fine_channel * grid_size;
           
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
            PRINTF_DEBUG("DEBUG : size = %d (%d x %d), fnorm = %e (counter sum = %.8f)\n", grid_size, width, height, fnorm, counter_sum);
            for (size_t i = 0; i < grid_size; i++) current_image[i] *= fnorm;
            fft_shift(current_image, grid_side, grid_side);
        }
    }

    //images_buffer.dump("images_after_shift.bin");

    auto tstop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(tstop - tstart).count();
    std::cout << "Imaging took " << duration << "s" << std::endl;
    // TODO : re-grid to SKY COORDINATES !!!
    // convert cos(alpha) to alpha - see notes !!!
    // how to do it ???

}

bool CPacerImager::CalculateUVW(){
    m_Baselines = m_MetaData.m_AntennaPositions.CalculateUVW(u_cpu, v_cpu, w_cpu, m_bIncludeAutos);
    PRINTF_INFO("INFO : calculated UVW coordinates of %d baselines\n", m_Baselines);
    u_max = u_cpu[0];
    v_max = v_cpu[0];

    for(int i {1}; i < u_cpu.size(); i++){
    if(u_max < u_cpu[i]) u_max = u_cpu[i];
    if(v_max < v_cpu[i]) v_max = v_cpu[i];
    }
    // Bacause we are also including conjugates at (-u,-v) UV point in gridding
    // u_min = -u_max and v_min = -v_max : was -35 / +35
    u_min = -u_max;
    //  u_max = +35;
    v_min = -v_max;
    //  v_max = +35;
    return (m_Baselines > 0);
}

/* this is useful logging code, find a good place for this 
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
    */


inline int wrap_index(int i, int side){
    if(i >= 0) return i % side;
    else return (side + i);
}

void CPacerImager::gridding_fast(Visibilities &xcorr, MemoryBuffer<std::complex<float>> &grids,
                                 MemoryBuffer<float> &grids_counters, double delta_u, double delta_v, int n_pixels,
                                 double min_uv /*=-1000*/,
                                 const char *weighting) {
    PRINTF_DEBUG("DEBUG : gridding : min_uv = %.4f\n", min_uv);

    bool bStatisticsCalculated = false;
    if (m_bPrintImageStatistics)
    { // TODO : ? may require a separate flag in the
      // future, for now just using a single
      // Statistics switch ON/OFF flag
    if(v_cpu){
        u_max = u_cpu[0];
        v_max = v_cpu[0];

      for(int i {1}; i < u_cpu.size(); i++){
        if(u_max < u_cpu[i]) u_max = u_cpu[i];
        if(v_max < v_cpu[i]) v_max = v_cpu[i];
      }
      // Bacause we are also including conjugates at (-u,-v) UV point in gridding
        // u_min = -u_max and v_min = -v_max : was -35 / +35
        u_min = -u_max;
        //  u_max = +35;
        v_min = -v_max;
        //  v_max = +35;
    } 

        bStatisticsCalculated = true;
    }
    memset(grids.data(), 0, grids.size() * sizeof(std::complex<float>));
    memset(grids_counters.data(), 0, grids_counters.size() * sizeof(float));
    int grid_size = n_pixels * n_pixels;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++)
    {
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++)
        {
            std::complex<float>* current_grid = grids.data() + time_step * xcorr.nFrequencies * grid_size + fine_channel * grid_size;
            float* current_counter = grids_counters.data() + time_step * xcorr.nFrequencies * grid_size + fine_channel * grid_size;
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
                                double u = u_cpu[ant1 * n_ant + ant2] / wavelength_m;

                                // 2022-09-24 : - removed for a test on MWA data
                                // 2022-09-09 - for now sticking to - sign here to have back
                                // comatible test EDA2 data
                                //  but looks like this (-) should not be here at least does not
                                //  match UV coverage from WSCEAN (it is flipped then see :
                                //  /home/msok/Desktop/PAWSEY/PaCER/logbook/20220826_image_simulation_part3.odt
                                double v = v_cpu[ant1 * n_ant + ant2] /
                                           wavelength_m; // the - sign here fixes the Y flip, but I am
                                                         // not sure why needed ??? check RTS :
                                                         // imagefromuv.c , LM_CopyFromFFT where some
                                                         // interesting re-shuffling happens too !
                                                         // - was for EDA2 data, but + is ok for MWA
                                                         // data - it may be consistent once I start
                                                         // using TMS equation 4.1 consistently for
                                                         // both EDA2 and MWA
                               
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

Images CPacerImager::gridding_imaging(Visibilities &xcorr, double delta_u, double delta_v, int n_pixels, double min_uv, const char *weighting) {
    printf("DEBUG : gridding_imaging( Visibilities& xcorr ) in pacer_imager.cpp\n");
    // allocates data structures for gridded visibilities:
    if(xcorr.on_gpu()) xcorr.to_cpu();
    size_t n_images {xcorr.integration_intervals() * xcorr.nFrequencies};
    size_t buffer_size{n_pixels * n_pixels * n_images};
    if(!grids_counters) grids_counters.allocate(buffer_size);
    if(!grids) grids.allocate(buffer_size);
    // Should be long long int, but keeping float now for compatibility reasons
    gridding_fast(xcorr, grids, grids_counters, delta_u, delta_v, n_pixels, min_uv, weighting);

    MemoryBuffer<std::complex<float>> images_buffer_float(buffer_size, false, false);
    PRINTF_INFO("PROGRESS : executing dirty image\n");
    dirty_image(grids, grids_counters, n_pixels, xcorr.integration_intervals(), xcorr.nFrequencies, images_buffer_float);
    return {std::move(images_buffer_float), xcorr.obsInfo, xcorr.nIntegrationSteps, xcorr.nAveragedChannels, static_cast<unsigned int>(n_pixels)};
}


/** 
    @brief run the imager
    
    @param xcorr: Visibilities to be imaged.
    @param n_pixels: Image side size.
    @param min_uv: mimimum UV length (what unit??)
    @param weighting: U for uniform or N for natural.
*/
Images CPacerImager::run_imager(Visibilities &xcorr, int n_pixels, double min_uv, const char *weighting){
    // TODO: make sure the xcorr structure has the same parameters at each run call
    // This will avoid to reallocate memory
    // TODO Cristian: time_step, fine_channel will be used in the to select a
    // subset of data to be imaged. ensures initalisation of object structures
    // TODO: this init function must be modified

    m_MetaData.fix_metafits(CObsMetadata::ux2gps(xcorr.obsInfo.startTime), 1.0);
    int n_ant = xcorr.obsInfo.nAntennas;
    int n_pol = xcorr.obsInfo.nPolarizations;
    // calculate UVW (if required)
    CalculateUVW();

    // xcorr.to_gpu(); // TODO: this should be gone!!
    if(!frequencies) frequencies.allocate(xcorr.nFrequencies);
    for(size_t fine_channel {0}; fine_channel < xcorr.nFrequencies; fine_channel++)
        frequencies[fine_channel] = this->get_frequency_hz(xcorr, fine_channel, COTTER_COMPATIBLE);
     
    if (apply_geom_correction)
        ApplyGeometricCorrections(xcorr, w_cpu, frequencies);

    if (apply_cable_correction){
        if(!cable_lengths) {
            cable_lengths.allocate(xcorr.obsInfo.nAntennas);
            for(size_t a {0}; a < xcorr.obsInfo.nAntennas;  a++)
                cable_lengths[a] = m_MetaData.m_AntennaPositions[a].cableLenDelta;
        }
        ApplyCableCorrections(xcorr, cable_lengths, frequencies);
    }

    //   // based on RTS : UV pixel size as function FOVtoGridsize in
    //   /home/msok/mwa_software/RTS_128t/src/gridder.c double frequency_hz =
    //   frequency_mhz*1e6; double wavelength_m = VEL_LIGHT / frequency_hz;
    // double FoV_radians = FOV_degrees * M_PI / 180.;
    // WARNING: it actually cancels out to end up being 1/PI :
    // TODO simplity + 1/(2FoV) !!! should be
    //  double delta_u = ( (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) /
    //  wavelength_m; // in meters (NOT WAVELENGHTS) double delta_v = (
    //  (VEL_LIGHT/frequency_hz)/(FOV_degrees*M_PI/180.) ) / wavelength_m; // in
    //  meters (NOT WAVELENGHTS)
    // TODO :
    // double delta_u = 1.00 / (FoV_radians); // should be 2*FoV_radians - see TMS etc
    // double delta_v = 1.00 / (FoV_radians); // Rick Perley page 16 :
                                           // /home/msok/Desktop/PAWSEY/PaCER/doc/Imaging_basics/ATNF2014Imaging.pdf
    
        // delta_u = 2.00*(u_max)/n_pixels; // Looks like this is
    // what it should be NOT u_max/wavelength_m . So delta_u must
    // be in meters here !!! It may all depend on calculation if
    // u_index see discussion in
    // /home/msok/Desktop/PAWSEY/PaCER/logbook/20240320_gridding_delta_u_meters_vs_wavelengths.odt
    double delta_u = 2.00 * (u_max) / n_pixels; 

    double delta_v = 2.00 * (v_max) / n_pixels; // and it's not ok because then delta_u is different
                                             // for both of them, which causes exp/shrink with freq

    // automatic calculation of pixel size in radians 1/(2u_max) - see Rick
    // Perley or just Lambda/B_max divide by 2 for oversampling.
    pixsize_in_radians = 1.00 / (2.00 * u_max); // does this one need to be /wavelength or not ???

    // NEW : based on desired image resolution
    // double delta_theta = (wavelength_m/35.0)/2.00; // 2 times oversampled
    // double delta_theta = ((230.0/300.00)/(2.00*35.00)); // use maximum
    // resolution oversampled by a factor of 2., at 230 MHz -> Lambda ~1.3043m
    // double delta_theta = m_ImagerParameters.m_PixsizeInRadians;
    // MWA TEST
    //     delta_u = 1.00/(n_pixels*delta_theta);
    //     delta_v = 1.00/(n_pixels*delta_theta);
    PRINTF_DEBUG("delta_u = %.8f (u_max = %.8f), delta_v = %.8f (v_max = %.8f), "
                    "calculated as 1/FoV = 1/(%d pixels * %.5f rad), delta_theta = %.5f "
                    "[deg]\n",
                    delta_u, u_max, delta_v, v_max, n_pixels, pixsize_in_radians, pixsize_in_radians * (180.00 / M_PI));

        // PRINTF_DEBUG("delta_u = %.8f , delta_v = %.8f , calculated
        // as 2.00*u_max/n_pixels, u_max = %.8f, n_pixels =
        // %d\n",delta_u,delta_v,u_max,n_pixels);

    // virtual function calls gridding and imaging in GPU/HIP version it is
    // overwritten and both gridding and imaging are performed on GPU memory :
    auto images = gridding_imaging(xcorr, delta_u, delta_v, n_pixels, min_uv, weighting);
    
    images.ra_deg = m_MetaData.raHrs*15.00;
    images.dec_deg = m_MetaData.decDegs;
    images.pixscale.resize(images.nFrequencies);
    for(size_t f {0}; f < images.nFrequencies; f++){
        double wavelength_m = VEL_LIGHT / frequencies[f];
        images.pixscale[f] = 1.00/(2.00*(u_max/wavelength_m));
    }
    return images;
}

double CPacerImager::get_frequency_hz(const Visibilities &vis, int fine_channel, bool cotter_compatible)
{
    // TODO: remove hardcoded values!!
    double fine_ch_bw = vis.obsInfo.frequencyResolution * vis.nAveragedChannels;
    double coarse_ch_bw = vis.obsInfo.coarseChannelBandwidth;
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

void CPacerImager::ApplyGeometricCorrections( Visibilities& xcorr, MemoryBuffer<float>& w, MemoryBuffer<double>& frequencies){
    apply_geometric_corrections_cpu(xcorr, w, frequencies);
}
   
void CPacerImager::ApplyCableCorrections( Visibilities& xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies){
    apply_cable_lengths_corrections_cpu(xcorr, cable_lengths, frequencies);
}
