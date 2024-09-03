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
#include "../src/corrections.h"
#include "../src/hip/corrections_gpu.h"
std::string dataRootDir;


void test_geometric_correction_cpu(){
    ObservationInfo obs_info {VCS_OBSERVATION_INFO};
    Visibilities xcorr = Visibilities::from_fits_file(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/01_before_geo_corrections.fits", obs_info);
    std::cout << "n_integrations = " << xcorr.integration_intervals() << ", n_frequencies = " << xcorr.nFrequencies << std::endl;
    CBgFits fits_vis_w;
    fits_vis_w.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/m_W.fits").c_str(), 0, 1, 1 );
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    apply_geometric_corrections_cpu(xcorr, fits_vis_w, frequencies);
    compare_xcorr_to_fits_file(xcorr, dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/01_after_geo_corrections.fits");
    std::cout << "'test_geometric_correction_cpu' passed." << std::endl;
}

void test_cable_lengths_correction_cpu(){

    ObservationInfo obs_info {VCS_OBSERVATION_INFO};
    Visibilities xcorr = Visibilities::from_fits_file(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/02_before_cable_corrections.fits", obs_info);
    std::cout << "n_integrations = " << xcorr.integration_intervals() << ", n_frequencies = " << xcorr.nFrequencies << std::endl;
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    MemoryBuffer<double> cable_lengths {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/cable_lengths.bin")};
    apply_cable_lengths_corrections_cpu(xcorr, cable_lengths, frequencies);
    compare_xcorr_to_fits_file(xcorr, dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/02_after_cable_corrections.fits");
    std::cout << "'test_cable_lengths_correction_cpu' passed." << std::endl;
}

#ifdef __GPU__
void test_geometric_correction_gpu(){
    ObservationInfo obs_info {VCS_OBSERVATION_INFO};
    Visibilities xcorr = Visibilities::from_fits_file(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/01_before_geo_corrections.fits", obs_info);
    std::cout << "n_integrations = " << xcorr.integration_intervals() << ", n_frequencies = " << xcorr.nFrequencies << std::endl;
    CBgFits fits_vis_w;
    fits_vis_w.ReadFits((dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/m_W.fits").c_str(), 0, 1, 1 );
    float *w_gpu;
    int xySize = xcorr.obsInfo.nAntennas * xcorr.obsInfo.nAntennas;
    gpuMalloc((void**)&w_gpu, xySize*sizeof(float));
    gpuMemcpy(w_gpu, fits_vis_w.get_data(), sizeof(float)*xySize,  gpuMemcpyHostToDevice);
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    apply_geometric_corrections_gpu(xcorr, w_gpu, frequencies);
    xcorr.to_cpu();
    gpuFree(w_gpu);
    compare_xcorr_to_fits_file(xcorr, dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/01_after_geo_corrections.fits");
    std::cout << "'test_geometric_correction_gpu' passed." << std::endl;
}

void test_cable_lengths_correction_gpu(){
    ObservationInfo obs_info {VCS_OBSERVATION_INFO};
    Visibilities xcorr = Visibilities::from_fits_file(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/02_before_cable_corrections.fits", obs_info);
    std::cout << "n_integrations = " << xcorr.integration_intervals() << ", n_frequencies = " << xcorr.nFrequencies << std::endl;
    MemoryBuffer<double> frequencies {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/frequencies.bin")};
    MemoryBuffer<double> cable_lengths {MemoryBuffer<double>::from_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/cable_lengths.bin")};
    apply_cable_lengths_corrections_gpu(xcorr, cable_lengths, frequencies);
    xcorr.to_cpu();
    compare_xcorr_to_fits_file(xcorr, dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/02_after_cable_corrections.fits");
    std::cout << "'test_cable_lengths_correction_gpu' passed." << std::endl;
}
#endif


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_geometric_correction_cpu();
        test_cable_lengths_correction_cpu();
        #ifdef __GPU__
        test_geometric_correction_gpu();
        test_cable_lengths_correction_gpu();
        #endif
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
