# comment TODO:
# I'd like to have installation of both _cpu and _gpu executables. So, I've added compilation of _cpu and _gpu version. I tried using:
# set_target_properties(pacer_dirty_imager PROPERTIES OUTPUT_NAME pacer_dirty_imager_gpu)
# install(TARGETS pacer_dirty_imager RUNTIME DESTINATION ${RUNTIME_DIR} RENAME pacer_dirty_imager_gpu)
# but it can only create pacer_dirty_imager_gpu or pacer_dirty_imager not both so making it slower, uglier, but the only working way so far (double compilation)
cmake_minimum_required(VERSION 3.1.0)
project(imaging VERSION 0.1 LANGUAGES CXX)
include(CTest)
set(CMAKE_CXX_STANDARD 14)

option(USE_HIP "Enable the GPU acceleration for AMD GPUs." OFF)
option(USE_CUDA "Enable the GPU acceleration for NVIDIA GPUs." OFF)

if(USE_CUDA)
    message("DEBUG : USE_CUDA=ON")
    enable_language(CUDA CXX)
    # default 70 for Topaz , 61 for my old laptop (Pascal), new laptop SM=86 
    if(DEFINED SM)
       message("CUDA compute capability passed as parameter SM = ${SM}")
    else()
       set(SM 86)
       message("CUDA compute capability not set -> using default SM ${SM}")
    endif()
    message("CUDA compute capability is ${SM}")
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES ${SM})
    endif()
elseif(USE_HIP)
   message("DEBUG : USE_HIP=ON")
   set(CMAKE_CXX_COMPILER hipcc)
endif()


file(GLOB sources "src/*.cpp")
file(GLOB headers "src/*.h")
file(GLOB apps "apps/*.cpp")
file(GLOB tests "tests/*.cpp")

if(USE_HIP OR USE_CUDA)
    add_definitions("-DIMAGER_HIP")
    file(GLOB gpu_sources "src/hip/*.cpp")
    file(GLOB gpu_headers "src/hip/*.hpp")    
endif()

if(USE_CUDA)
    set_source_files_properties( ${gpu_sources} PROPERTIES LANGUAGE CUDA)
endif()

# ${gpu_sources} included into library imaging. This is required to link blink_pipeline on Setonix.
# but this may not work correctly with NVIDIA CUDA as it does not allow kernels to go to libraries
# TO BE VERIFIED and TESTED
add_library(imaging SHARED ${sources} ${gpu_sources})
set_target_properties(imaging PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
    PUBLIC_HEADER "${headers}") 

configure_file(pacerimagerlib.pc.in pacerimagerlib.pc @ONLY)


if(NOT RUNTIME_DIR)
    set(RUNTIME_DIR ${CMAKE_INSTALL_PREFIX}/bin/)
endif(NOT RUNTIME_DIR)

# Required for msfitslib includes (otherwise compilation error - I am not sure why did not investigate)
add_definitions(-D_UNIX)

# define macro for profiler:
IF(ENABLE_PROFILER)
    ADD_DEFINITIONS(-D_PACER_PROFILER_ON_)
ENDIF(ENABLE_PROFILER)

install(TARGETS imaging 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# This tries to install to /pkgconfig - ${CMAKE_INSTALL_DATAROOTDIR} not defined?
# install(FILES ${CMAKE_BINARY_DIR}/pacerimagerlib.pc DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

# External dependencies

find_package(PkgConfig REQUIRED)
find_library(CFITSIO_LIB cfitsio HINTS ENV LD_LIBRARY_PATH REQUIRED)

find_library(FFTW3_LIB fftw3 HINTS ENV LD_LIBRARY_PATH REQUIRED)

find_library(LIBNOVA_LIB nova ENV LD_LIBRARY_PATH REQUIRED)

find_library(MSFITSLIB NAMES msfitslib HINTS ENV LD_LIBRARY_PATH REQUIRED)

find_library(BLINK_ASTROIO blink_astroio HINTS ENV LD_LIBRARY_PATH)
if (NOT BLINK_ASTROIO)
   message(FATAL_ERROR "BLINK_ASTROIO library not found.")
endif()

find_library(LIBPAL_LIB pal PATHS ENV LD_LIBRARY_PATH REQUIRED)
if(NOT LIBPAL_LIB)
   message(FATAL_ERROR "LIBPAL_LIB library not found.")
endif()

# TODO add here code to find cufft or rocfft
if(USE_CUDA)
   set(GPU_LIBS "-lcufft")
elseif(USE_HIP)
   set(GPU_LIBS "-lhipfft")
else()
   set(GPU_LIBS "")
endif()

target_link_libraries(imaging ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${GPU_LIBS}  -lpthread -ldl)

# temporary changed back to NOT LINKING imaging library due to 
# see Issue #3 ( program crashes at the end when linking using imaging library )
# Very likely due to some global variable which behaves differently when linking share library.
# add_executable(pacer_dirty_imager apps/pacer_imager_main.cpp) 
# add_executable(pacer_imager_test apps/pacer_imager_test.cpp)
# add_executable(test_calsol apps/test_calsol/main.cpp)
# target_link_libraries(pacer_dirty_imager imaging)
# target_link_libraries(pacer_imager_test imaging)
# target_link_libraries(test_calsol imaging)
