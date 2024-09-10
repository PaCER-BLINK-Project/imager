#include <exception>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <stdexcept>

#include <astroio.hpp>
#include <gpu_macros.hpp>
#include "common.hpp"

#include "../src/hip/gridding_gpu.h"
std::string dataRootDir;


void test_gridding_gpu(){
    ObservationInfo obs_info {VCS_OBSERVATION_INFO};
    Visibilities xcorr = Visibilities::from_fits_file(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/03_after_geo_corrections.fits", obs_info);
    xcorr.to_gpu();
    CBgFits fits_vis_u, fits_vis_v;
    fits_vis_u.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/U.fits").c_str(), 0, 1, 1 );
    fits_vis_v.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/V.fits").c_str(), 0, 1, 1 );
  
    int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    MemoryBuffer<int> antenna_flags {MemoryBuffer<int>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/antenna_flags.bin")};
    MemoryBuffer<float> antenna_weights {MemoryBuffer<float>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/antenna_weights.bin")};
    antenna_flags.to_gpu();
    antenna_weights.to_gpu();

    double delta_u = 1.221977710723877;
    double delta_v = 1.1040256023406982;
    int n_pixels = 8192;
    int min_uv = -1000;

    size_t n_images{xcorr.integration_intervals() * xcorr.nFrequencies};
    size_t buffer_size {n_pixels * n_pixels * n_images};
    MemoryBuffer<float> grids_counters_buffer(buffer_size, true);
    MemoryBuffer<std::complex<float>> grids_buffer(buffer_size,  true);
    gridding_gpu(xcorr, -1, -1, fits_vis_u, fits_vis_v, antenna_flags.data(), antenna_weights.data(), frequencies,
      delta_u, delta_v, n_pixels, min_uv, grids_counters_buffer, grids_buffer);

   grids_counters_buffer.to_cpu();
   grids_buffer.to_cpu();
   CBgFits reference_grid, reference_grid_counter;
    reference_grid.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/uv_grid_real_8192x8192.fits").c_str(), 0, 1, 1 );
    reference_grid_counter.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/uv_grid_counter_8192x8192.fits").c_str(), 0, 1, 1 );

    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(grids_counters_buffer[i] != reference_grid_counter.getXY(i % n_pixels, i / n_pixels)){
            std::cerr << "Error!! Counters are not the same at position " << i << ": " << grids_counters_buffer[i] << " != " << reference_grid_counter.getXY(i % n_pixels, i / n_pixels) << std::endl;
            break;
        }
    }
    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(std::abs(grids_buffer[i].real() - reference_grid.getXY(i % n_pixels, i / n_pixels)) > 1e-4){
            std::cerr << "Error!! Grids are not the same at position " << i << ": " << grids_buffer[i].real() << " != " << reference_grid.getXY(i % n_pixels, i / n_pixels) << std::endl;
            exit(1);
        }
    }

    std::cout << "'test_gridding_gpu' passed." << std::endl;
}




int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        #ifdef __GPU__
        test_gridding_gpu();
        #endif
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
