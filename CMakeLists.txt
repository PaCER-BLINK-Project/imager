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

target_link_libraries(imaging ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} -lpthread -ldl)

# temporary changed back to NOT LINKING imaging library due to 
# see Issue #3 ( program crashes at the end when linking using imaging library )
# Very likely due to some global variable which behaves differently when linking share library.
# add_executable(pacer_dirty_imager apps/pacer_imager_main.cpp) 
# add_executable(pacer_imager_test apps/pacer_imager_test.cpp)
# add_executable(test_calsol apps/test_calsol/main.cpp)
# target_link_libraries(pacer_dirty_imager imaging)
# target_link_libraries(pacer_imager_test imaging)
# target_link_libraries(test_calsol imaging)
add_executable(pacer_imager_test apps/pacer_imager_test.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
add_executable(test_calsol apps/test_calsol/main.cpp src/apply_calibration.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
target_link_libraries(pacer_imager_test ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl)

# Test FFTW program :
add_executable(fftw_blocks apps/fftw_blocks.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
target_link_libraries(fftw_blocks ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl -lfftw3f_threads -lfftw3f)

target_link_libraries(test_calsol ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl)
if(USE_HIP OR USE_CUDA)
    if(USE_CUDA)
       set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_${SM},code=sm_${SM})
    endif()

    # TODO : should use something like ASTROIO_LIB path but nothing like this is defined after loading this module ...
    # include_directories("/software/projects/director2183/msok/setonix/2023.08/development/blink_astroio/master/include/")

    add_executable(pacer_dirty_imager apps/pacer_imager_main.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
    add_executable(pacer_dirty_imager_gpu apps/pacer_imager_main.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
    if(USE_CUDA)
        target_compile_options(pacer_dirty_imager PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_NVCC_FLAGS}>)
        target_compile_options(pacer_dirty_imager_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_NVCC_FLAGS}>)
    endif()
    target_link_libraries(pacer_dirty_imager ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl ${GPU_LIBS})
    target_link_libraries(pacer_dirty_imager_gpu ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl ${GPU_LIBS})

    # Gayatri's multi-block version:
    add_executable(cufft_blocks apps/cufft_blocks.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
    if(USE_CUDA)
       set_source_files_properties( apps/cufft_blocks.cpp PROPERTIES LANGUAGE CUDA)
       target_compile_options( cufft_blocks PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_NVCC_FLAGS}>)
    endif()
    target_link_libraries(cufft_blocks ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl ${GPU_LIBS})

    # gridding debugger 
    if(DEBUG_GRIDDING)
      add_executable(debug_gpu_gridding debug/debug_gpu_gridding.cu src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
      target_link_libraries(debug_gpu_gridding ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl ${GPU_LIBS})
      message("INFO : debug_gridding program will be build")
    endif()

else()
  add_executable(pacer_dirty_imager apps/pacer_imager_main.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
  target_link_libraries(pacer_dirty_imager ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl)

  add_executable(pacer_dirty_imager_cpu apps/pacer_imager_main.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
  target_link_libraries(pacer_dirty_imager_cpu ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${BLINK_ASTROIO} -lpthread -ldl)
endif()

if(BUILD_LFILE_IMAGER)
   # Multi-block imager :
   find_library(LFILEREADER NAMES liblfile2corrmatrixlib.so HINTS ENV LFILEREADER_DIR PATHS ENV LFILEREADER_DIR PATH_SUFFIXES "lib" REQUIRED)
   if(NOT LFILEREADER)
      message(FATAL_ERROR "LFILEREADER library not found")
   endif()

   add_executable(lfile_imager_test apps/lfile_imager_test.cpp src/image_analyser.cpp src/realtime_image_analyser.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
   if(USE_HIP)
       target_link_libraries(lfile_imager_test ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${LFILEREADER} -lpthread -ldl ${GPU_LIBS}) # -llfile2corrmatrixlib)

       add_executable(pacer_imager_multi apps/pacer_imager_multi.cpp src/image_analyser.cpp src/realtime_image_analyser.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp ${gpu_sources})
       target_link_libraries(pacer_imager_multi ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${LFILEREADER} -lpthread -ldl ${GPU_LIBS}) # -llfile2corrmatrixlib)
   else()
       target_link_libraries(lfile_imager_test ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${LFILEREADER} -lpthread -ldl) # -llfile2corrmatrixlib)
       add_executable(pacer_imager_multi apps/pacer_imager_multi.cpp src/image_analyser.cpp src/realtime_image_analyser.cpp src/antenna_positions.cpp src/pacer_imager.cpp src/pacer_imager_thread.cpp src/pacer_common.cpp src/apply_calibration.cpp src/pacer_imager_parameters.cpp src/observation_metadata.cpp)
       target_link_libraries(pacer_imager_multi ${MSFITSLIB} ${LIBNOVA_LIB} ${CFITSIO_LIB} ${FFTW3_LIB} ${LIBPAL_LIB} ${LFILEREADER} -lpthread -ldl) # -llfile2corrmatrixlib)
   endif()
else()
   message("L-file imager is not going to be build (use -DBUILD_LFILE_IMAGER=ON to build it too)")
endif()

# TEST
enable_testing()
find_program (BASH_PROGRAM bash)
add_test(NAME blink_test_eda2_data_20200209.sh COMMAND ../tests/blink_test_eda2_data_20200209.sh)
add_test(NAME blink_test_mwa_data_obsid1103645160.sh COMMAND ../tests/blink_test_mwa_data_obsid1103645160.sh)
add_test(NAME blink_test_mwa_simul_hydra.sh COMMAND ../tests/blink_test_mwa_simul_hydra.sh)
if(BUILD_LFILE_IMAGER)
   add_test(NAME blink_test_eda2_data_20221216.sh COMMAND ../tests/blink_test_eda2_data_20221216.sh)
   if(USE_HIP)
      add_test(NAME blink_test_eda2_data_20221216_gpu_n1.sh COMMAND ../tests/blink_test_eda2_data_20221216_gpu_n1.sh)
      add_test(NAME blink_test_eda2_data_20221216_gpu_n2.sh COMMAND ../tests/blink_test_eda2_data_20221216_gpu_n2.sh)
   endif()
   install(TARGETS pacer_imager_multi RUNTIME DESTINATION ${RUNTIME_DIR})
endif()

message("Runtime files will be installed in ${RUNTIME_DIR}")

# install(PROGRAMS "tests/blink_test_eda2_data_20200209.sh" DESTINATION ${RUNTIME_DIR})
install(FILES "tests/blink_test_eda2_data_20200209.sh" DESTINATION ${RUNTIME_DIR})
install(FILES "tests/blink_test_mwa_simul_hydra.sh" DESTINATION ${RUNTIME_DIR})
install(FILES "tests/blink_test_mwa_data_obsid1103645160.sh" DESTINATION ${RUNTIME_DIR})
install(TARGETS pacer_dirty_imager RUNTIME DESTINATION ${RUNTIME_DIR})

if(USE_HIP OR USE_CUDA)
   install(TARGETS pacer_dirty_imager_gpu RUNTIME DESTINATION ${RUNTIME_DIR})
else()
   install(TARGETS pacer_dirty_imager_cpu RUNTIME DESTINATION ${RUNTIME_DIR})
endif()
