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

std::string dataRootDir;


void test_geometric_correction_cpu(){

    throw TestFailed( "'test_cable_lengths_correction_cpu': implement me.");
    std::cout << "'test_geometric_correction_cpu' passed." << std::endl;
}

void test_cable_lengths_correction_cpu(){

    throw TestFailed( "'test_cable_lengths_correction_cpu': implement me.");
    std::cout << "'test_cable_lengths_correction_cpu' passed." << std::endl;
}



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
    } catch (std::exception& ex){
        std::cerr << ex.what() << std::endl;
        return 1;
    }
    
    std::cout << "All tests passed." << std::endl;
    return 0;
}
