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
#include "../src/pacer_imager.h"

#ifdef __GPU__
#include "../src/gpu/pacer_imager_hip.h"
#endif
std::string dataRootDir;



void test_imager_common(CPacerImager& imager, bool is_cpu){

    std::string vis_file {dataRootDir + "/mwa/1192477696/1192477696_1192477703_ch109.fits"};
    std::string antennaPositionsFile {""};
    std::string output_dir { is_cpu ? 
        "test_imager_conv_cpu/" : "test_imager_conv_gpu/"};

   ObservationInfo obs_info {VCS_OBSERVATION_INFO};
   auto xcorr = Visibilities::from_fits_file(vis_file, obs_info);
   if(!is_cpu) xcorr.to_gpu();
   imager.m_nConvolvingKernelSize = 10;
   auto images = imager.run(xcorr);
   std::cout << "Saving images to disk..." << std::endl;
   images.to_fits_files(output_dir);
}

void test_imager_cpu(){
    std::string szWeighting {"N"};
    const int image_size = 512;
    double MinUV = -1000;
    std::string metadataFile {dataRootDir + "/mwa/1192477696/1192477696.metafits"};
    std::vector<int> flagged_antennas {};
    CPacerImager imager {metadataFile, image_size, flagged_antennas, true, Polarization::XX, 2.0f, MinUV, szWeighting.c_str()};
    test_imager_common(imager, true);
}

#ifdef __GPU__
void test_imager_gpu(){
    std::cout << "This is a test with Nathan !!!!!!!!." << std::endl;
    std::string szWeighting {"N"};
    const int image_size = 512;
    double MinUV = -1000;
    std::string metadataFile {dataRootDir + "/mwa/1192477696/1192477696.metafits"};
    std::vector<int> flagged_antennas {};
    CPacerImagerHip imager {metadataFile, image_size, flagged_antennas, true, Polarization::XX, 2.0f, MinUV, szWeighting.c_str()};
    test_imager_common(imager, false);
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
        test_imager_cpu();
#ifdef __GPU__
        test_imager_gpu();
#endif
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
