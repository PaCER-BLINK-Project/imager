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

std::string dataRootDir;


template <typename T>
bool complex_vectors_equal(const std::complex<T>* a, const std::complex<T>* b, size_t length){
    double delta;
    const double TOL {1e-5};
    for(size_t i {0}; i < length; i++){
        if (std::abs(a[i]) == 0) 
            delta = std::abs(b[i]);
        else 
            delta = std::abs(a[i] - b[i]);
        
        if (delta > TOL) {
            std::cout << "Elements at position " << i << " differs (delta = " << delta <<"): " << "a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
            return false;
        }
    }
    return true;
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

void load_dump(std::string filename, char *& buffer, size_t& size){
    std::ifstream infile (filename, std::ifstream::binary);
    // get size of file
    infile.seekg (0,infile.end);
    size = infile.tellg();
    infile.seekg (0);

    // allocate memory for file content
    buffer = new char[size];

    // read content of infile
    infile.read (buffer, size);
    infile.close();
}



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




void test_fft_shift_on_reference(){

    char *input, *ref_output;
    size_t input_size, ref_output_size;
    const size_t image_side {8192ul};
    load_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/image_before_shift.bin", input, input_size);
    load_dump(dataRootDir + "/mwa/1276619416/imager_stages/1s_ch000/image_after_shift.bin", ref_output, ref_output_size);

    std::complex<float> *input_cpx {reinterpret_cast<std::complex<float>*>(input)};
    std::complex<float> *ref_output_cpx {reinterpret_cast<std::complex<float>*>(ref_output)};

    fft_shift(input_cpx, image_side, image_side);

    if(!complex_vectors_equal(input_cpx, ref_output_cpx, image_side * image_side)){
        throw TestFailed("test_fft_shift_on_reference: vectors are not equal!");
    }
    
    std::cout << "Test 'test_fft_shift_on_reference' passed." << std::endl;
}




/*

void test_imager_cpu_one_coarse_channel(){    
    auto input_visibilities = Visibilities::from_fits_file(dataRootDir + "/imager/vis_coarse_channel.fits", VCS_OBSERVATION_INFO);


    const unsigned int n_antennas {128u};
    const unsigned int n_baselines {(n_antennas + 1) * (n_antennas / 2)};
    const unsigned int n_polarisations {2u};
    const unsigned int n_fine_channels {128u};
    const unsigned int n_time_samples {1u};
    const unsigned int n_integrated_samples {52u};
    const unsigned int n_integration_intervals {n_time_samples / n_integrated_samples};
    // the following definition will make sure that the output won't be scaled by the time
    // averaging factor.
    const double time_resolution {1.0 / n_integrated_samples};
    const unsigned int n_channels_to_avg {1u};
    const unsigned int reset_visibilities {1u};
    
    size_t n_voltages {static_cast<size_t>(n_integration_intervals) * n_fine_channels * n_antennas * n_polarisations * n_integrated_samples};
    size_t n_visibilities {static_cast<size_t>(n_integration_intervals) * n_fine_channels * n_baselines * n_polarisations * n_polarisations};

    delete[] visibilities_cpu;
    delete[] outputData;
    std::cout << "'test_correlation_with_xgpu_in_mwax_data' passed." << std::endl;
}

*/




int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_fft_shift_simple();
        // the following does not work anymore the reference is in double
        //test_fft_shift_on_reference();
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
