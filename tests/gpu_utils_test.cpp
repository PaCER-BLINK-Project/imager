#include <exception>
#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <stdexcept>
#include <sstream>

#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <gpu_macros.hpp>
#include "common.hpp"
#include "../src/gpu/gpu_utils.h"
#include "../src/pacer_imager.h" // Just for the Images class. Then, we will move the class in AstroIO
std::string dataRootDir;


void test_vector_sum(){
    const int n_images {5};
    const int image_size {3000};
    float *data_cpu {new float[image_size * n_images]};
    float *sums_cpu {new float[n_images]};
    memset(sums_cpu, 0, sizeof(float) * n_images);
    
    for(int img_id {0}; img_id < n_images; img_id++){
        for(int p {0}; p < image_size; p++){
            float val = rand() % 50;
            // Initialise array
            data_cpu[img_id * image_size + p] = val;
            // and compute the result on CPU
            sums_cpu[img_id] += val;
        }
    }

    // Now the GPU implementation should return the same result
    float* data_gpu;
    gpuMalloc(&data_gpu, sizeof(float) * n_images * image_size);
    gpuMemcpy(data_gpu, data_cpu, sizeof(float) * n_images * image_size, gpuMemcpyHostToDevice);
    MemoryBuffer<float> results {n_images, MemoryType::DEVICE};
    vector_sum_gpu(data_gpu, image_size, n_images, results);
    results.to_cpu();
    // check the values
    for(int img_id {0}; img_id < n_images; img_id++){
        if(sums_cpu[img_id] != results[img_id]){
            std::stringstream ss;
            ss << "'test_vector_sum' failed: sums_cpu[" << img_id <<"] (" \
                <<  sums_cpu[img_id]  << ") != results[" << img_id << "] (" << results[img_id] << ")" << std::endl;
            throw TestFailed(ss.str().c_str());
        }
    }
    gpuFree(data_gpu);
    delete[] data_cpu;
    delete[] sums_cpu;
    std::cout << "'test_vector_sum' passed." << std::endl;
}


void test_fft_shift_and_norm(){
    std::complex<float> input_cpx_one[] {
        {0, 0}, {1, 0}, {2, 0}, {3, 0},
        {0, 0}, {1, 0}, {2, 0}, {3, 0},
        {0, 0}, {1, 0}, {2, 0}, {3, 0}
    };
    std::complex<float> ref_output_one[] {
        {3, 0}, {2, 0}, {1, 0}, {0, 0},
        {3, 0}, {2, 0}, {1, 0}, {0, 0},
        {3, 0}, {2, 0}, {1, 0}, {0, 0}
    
    };
    MemoryBuffer<std::complex<float>> input_gpu {12};
    for(int i {0}; i < 12; i++) input_gpu[i] = input_cpx_one[i];
    MemoryBuffer<float> fnorm {3};
    fnorm[0] = 1;
    fnorm[1] = 1;
    fnorm[2] = 1;
    fnorm.to_gpu();
    input_gpu.to_gpu();
    fft_shift_and_norm_gpu(reinterpret_cast<gpufftComplex*>(input_gpu.data()), 2, 2, 3, fnorm);
    input_gpu.to_cpu();
    for(int x = 0; x < 12; x++)
        if(input_gpu[x] != ref_output_one[x])
            throw TestFailed("'test_fft_shift_and_norm': wrong output.");
    
    std::cout << "'test_fft_shift_and_norm' passed." << std::endl;
}



void test_averaging_kernel(){

    const int n_images  {4};
    const int image_side {512};
    const int n_pixels {image_side * image_side};

    MemoryBuffer<std::complex<float>> input_images_data {n_images * n_pixels};
    MemoryBuffer<std::complex<float>> reference_avg_data {n_pixels};
    memset(reference_avg_data.data(), 0, sizeof(std::complex<float>) * n_pixels);

    for(int img_id {0}; img_id < n_images; img_id++){
        for(int i {0}; i < n_pixels; i++){
            float val = i % 100;
            input_images_data[img_id * n_pixels + i].real(val);
            reference_avg_data[i].real(reference_avg_data[i].real() + val);
        }
    }
    for(int i {0}; i < n_pixels; i++) reference_avg_data[i] /= n_images;
    input_images_data.to_gpu();
    Images input_images  {std::move(input_images_data), VCS_OBSERVATION_INFO, 2, 2, image_side, 0, 0, 0, 0};
    Images avg_image = image_averaging_gpu(input_images);
    std::cerr << "Invocation over." << std::endl;
    avg_image.to_cpu();
    gpuDeviceSynchronize();
    std::complex<float> *out_avg_data {avg_image.data()};
    for(int i {0}; i < n_pixels; i++){
        if(out_avg_data[i] != reference_avg_data[i]){
            std::cout << "'test_averaging_kernel' error: output differs at position " << i << \
                ": out_avg_data = " << out_avg_data[i] << ", reference = " << reference_avg_data[i] << std::endl;
            throw TestFailed("'test_averaging_kernel' failed. Output differs.");
        }
    }
    std::cout << "'test_averaging_kernel' passed." << std::endl;
}


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_vector_sum();
        test_fft_shift_and_norm();
        test_averaging_kernel();
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
