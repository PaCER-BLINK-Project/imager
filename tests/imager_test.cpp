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
#include "../src/hip/pacer_imager_hip.h"

std::string dataRootDir;


void test_fft_shift_simple(){
    std::complex<float> input_cpx_one[] {{0, 0}, {1, 0}, {2, 0}, {3, 0}};
    std::complex<float> ref_output_one[] {{3, 0}, {2, 0}, {1, 0}, {0, 0}};
    fft_shift(input_cpx_one, 2, 2);
    for(int x = 0; x < 4; x++)
        if(input_cpx_one[x] != ref_output_one[x])
            throw TestFailed("'test_fft_shift_simple: test 1, wrong output.");
    
    // now test a odd dimension
    std::complex<float> input_cpx_two[] {{0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}};
    std::complex<float> ref_output_two[] {{5, 0}, {4, 0}, {3, 0}, {2, 0}, {1, 0}, {0, 0}};
    fft_shift(input_cpx_two, 3, 2);
    for(int x = 0; x < 6; x++)
        if(input_cpx_two[x] != ref_output_two[x]){
            std::cout << "Elements differ at pos " << x << ": " << input_cpx_two[x] << " != " << ref_output_two[x] << std::endl;
            throw TestFailed("'test_fft_shift_simple: test 2, wrong output.");
        }
    
    std::cout << "'test_fft_shift_simple' passed." << std::endl;
}




void test_imager_common(CPacerImager& imager, bool is_cpu){

    std::string vis_file {dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/input_visibilities.fits"};
    std::string metadataFile {dataRootDir + "/mwa/1276619416/20200619163000.metafits"};
    std::string antennaPositionsFile {""};
    std::string output_dir { is_cpu ? 
        "/scratch/director2183/cdipietrantonio/test_imager_cpu" : "/scratch/director2183/cdipietrantonio/test_imager_gpu"};
    std::string szWeighting {"N"};
    const int image_size = 256;
    bool bZenithImage {false};
    double fUnixTime {1592584240};
    double MinUV = -1000;
    double FOV_degrees = 30;
    std::vector<int> flagged_antennas {21, 25, 58, 71, 80, 81, 92, 101, 108, 114, 119, 125};

    CImagerParameters::m_bApplyCableCorr = true;
    CImagerParameters::m_bApplyGeomCorr = true;
    imager.m_ImagerParameters.m_bConstantUVW = true; 
    imager.m_ImagerParameters.SetGlobalParameters(antennaPositionsFile.c_str(), bZenithImage); // Constant UVW when zenith image (-Z)
    imager.m_ImagerParameters.m_szOutputDirectory = output_dir.c_str();


    if(strlen(metadataFile.c_str())){
        imager.m_ImagerParameters.m_MetaDataFile = metadataFile.c_str();
        imager.m_ImagerParameters.m_bConstantUVW = false; // when Meta data file is provided it assumes that it will pointed observation (not all sky)
    }
    
    imager.m_ImagerParameters.m_fUnixTime = fUnixTime;
    
   imager.Initialise(0);
   // setting flagged antennas must be called / done after reading METAFITS file:
   if( flagged_antennas.size() > 0 ){
       imager.SetFlaggedAntennas( flagged_antennas );
   }
   ObservationInfo obs_info {VCS_OBSERVATION_INFO};
   auto xcorr = Visibilities::from_fits_file(vis_file, obs_info);
   auto images = imager.run_imager(xcorr, -1, -1, image_size, FOV_degrees, MinUV, true, true, szWeighting.c_str(), output_dir.c_str(), false);
   std::cout << "Saving images to disk..." << std::endl;
   images.to_fits_files(output_dir);
}

void test_imager_cpu(){
    CPacerImager imager;
    test_imager_common(imager, true);
}


void test_imager_gpu(){
    CPacerImagerHip imager;
    test_imager_common(imager, false);
}

int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_fft_shift_simple();
        test_imager_cpu();
        test_imager_gpu();
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}