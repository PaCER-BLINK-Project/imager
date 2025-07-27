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
    MemoryBuffer<float> u_buff {MemoryBuffer<float>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/u_buff.bin")};
    MemoryBuffer<float> v_buff {MemoryBuffer<float>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/v_buff.bin")};    
    u_buff.to_gpu();
    v_buff.to_gpu();
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    MemoryBuffer<int> antenna_flags {MemoryBuffer<int>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/antenna_flags.bin")};
    MemoryBuffer<float> antenna_weights {MemoryBuffer<float>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/antenna_weights.bin")};
    antenna_flags.to_gpu();
    antenna_weights.to_gpu();
    double delta_u = 39.10328674, delta_v = 35.32881927;
    int n_pixels = 256;
    int min_uv = -1000;

    size_t n_images{xcorr.integration_intervals() * xcorr.nFrequencies};
    size_t buffer_size {n_pixels * n_pixels * n_images};
    MemoryBuffer<float> grids_counters(buffer_size, true);
    MemoryBuffer<std::complex<float>> grids(buffer_size,  true);
    gridding_gpu(xcorr, -1, -1, u_buff, v_buff, antenna_flags.data(), antenna_weights.data(), frequencies,
      delta_u, delta_v, n_pixels, min_uv, grids_counters, grids);

    grids_counters.to_cpu();
    grids.to_cpu();
    MemoryBuffer<std::complex<float>> reference_grid {MemoryBuffer<std::complex<float>>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/grids_buffer.bin")};
    MemoryBuffer<float> reference_grid_counter {MemoryBuffer<float>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/grids_counters_buffer.bin")};
    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(grids_counters[i] != reference_grid_counter[i]){
            std::cerr << "Error!! Counters are not the same at position " << i << ": " << grids_counters[i] << " != " << reference_grid_counter[i] << std::endl;
            throw TestFailed("'test_gridding_gpu' failed: counters are not the same.");
        }
    }
    for(size_t i {0}; i < n_pixels * n_pixels; i++){
        if(std::abs(grids[i].real() - reference_grid[i].real()) > 1e-3 || std::abs(grids[i].imag() - reference_grid[i].imag()) > 1e-3){
            std::cerr << "Error!! Grids are not the same at position " << i << ": " << grids[i] << " != " << reference_grid[i] << std::endl;
            throw TestFailed("'test_gridding_gpu' failed: grids are not the same.");
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
