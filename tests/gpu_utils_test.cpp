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
#include <gpu_macros.hpp>
#include "common.hpp"
#include "../src/hip/gpu_utils.h"

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
    float* results_gpu = sum_gpu_atomicadd(data_gpu, image_size, n_images);
    float* results {new float[image_size * n_images]};
    gpuMemcpy(results, results_gpu, sizeof(float) * n_images, gpuMemcpyDeviceToHost);

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


int main(void){
    char *pathToData {std::getenv(ENV_DATA_ROOT_DIR)};
    if(!pathToData){
        std::cerr << "'" << ENV_DATA_ROOT_DIR << "' environment variable is not set." << std::endl;
        return -1;
    }
    dataRootDir = std::string {pathToData};

    try{
        test_vector_sum();
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
