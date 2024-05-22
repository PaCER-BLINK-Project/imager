#!/bin/bash -e

# WARNING : this script is assumed to be for a supercomputer (HPC) : Topaz, Setonix or Garrawarla etc - not for laptop or desktop as there are "module" commands
#           for laptop desktop just do : cmake .. CMAKE_OPTIONS , make , make test , sudo make install

build_type="cpu"
build_dir=build
if [[ -n "$1" && "$1" != "-" ]]; then
   build_type=$1
fi

cmake_options=""
# DEBUG : -DCMAKE_BUILD_TYPE=Debug
if [[ $build_type == "gpu" ]]; then
   build_dir=build_gpu
fi

build_opt="Release" # or  "Debug"
if [[ -n "$2" && "$2" != "-" ]]; then
   build_opt="$2"
fi


# default options decided based on build type can be over-written by the 2nd parameter:
if [[ -n "$3" && "$3" != "-" ]]; then
   cmake_options=$3
else
   if [[ $build_type == "gpu" ]]; then
      # WARNING : -DCUDA_SM=70 is specific for Topaz, the default is sm=61 for msok's laptop, but laptop/desktop is compiled using normal cmake 
      cmake_options="-DUSE_HIP=ON -DCUDA_SM=70"
      echo "INFO (build.sh) : build type ${build_type} detected -> added cmake_options = ${cmake_options}"
   fi
fi

version=""
if [[ -n "$4" && "$4" != "-" ]]; then
   version="$4"
   build_dir=${build_dir}_${version}
fi

dotests=1
if [[ -n "$5" && "$5" != "-" ]]; then
   dotests=$5
fi

# First, you need to source the bash library
module load bash-utils
echo "source ${BASH_UTILS_DIR}/build_utils.sh"
source "${BASH_UTILS_DIR}/build_utils.sh"


PROGRAM_NAME=blink-imager-${build_type}
if [[ -n "$6" && "$6" != "-" ]]; then
   PROGRAM_NAME="$6"
fi

PROGRAM_VERSION=devel
if [[ -n "$7" && "$7" != "-" ]]; then
   PROGRAM_VERSION="$7"
fi


echo "############################################"
echo "PARAMETERS (build.sh scripts) :"
echo "############################################"
echo "PROGRAM_NAME = $PROGRAM_NAME"
echo "PROGRAM_VERSION = $PROGRAM_VERSION"
echo "############################################"


 
# the following function sets up the installation path according to the
# cluster the script is running on and the first argument given. The argument
# can be:
# - "group": install the software in the group wide directory
# - "user": install the software only for the current user
# - "test": install the software in the current working directory 
echo "process_build_script_input user"
process_build_script_input user 


# load all the modules required for the program to compile and run.
# the following command also adds those module names in the modulefile
# that this script will generate.
echo "Loading required modules ..."
echo "Loading modules for PAWSEY_CLUSTER = $PAWSEY_CLUSTER"
if [ $PAWSEY_CLUSTER = "setonix" ]; then
   module reset
   # 12.1.0 -> 12.2.0
   module use /software/projects/director2183/setonix/modules/zen3/gcc/12.2.0  /software/projects/director2183/msok/setonix/modules
   # was msfitslib/master-b37lvzx -> msfitslib/devel
   
   # WARNING : module_load has to be a single command in build.sh script as it overwrite variable BU_MODULELOADED !
   count=`echo $cmake_options | grep "BUILD_LFILE_IMAGER=ON" | wc -l`
   # pal/0.9.8-pwmda33 -> pal/0.9.8-fwf66s4
   if [[ $count -gt 0 ]]; then
       echo "module_load blink_test_data/devel cfitsio/4.1.0 msfitslib/devel blink_astroio_msok/master fftw/3.3.10 pal/0.9.8-fwf66s4 libnova/0.15.0-hpvrrfr rocm/5.4.3 lfile2corrmatrix/devel"
       module_load blink_test_data/devel cfitsio/4.1.0 msfitslib/devel blink_astroio_msok/master fftw/3.3.10 pal/0.9.8-fwf66s4 libnova/0.15.0-hpvrrfr rocm/5.4.3 lfile2corrmatrix/devel
   else
       echo "module_load blink_test_data/devel cfitsio/4.1.0 msfitslib/devel blink_astroio_msok/master fftw/3.3.10 pal/0.9.8-fwf66s4 libnova/0.15.0-hpvrrfr rocm/5.4.3"
       module_load blink_test_data/devel cfitsio/4.1.0 msfitslib/devel blink_astroio_msok/master fftw/3.3.10 pal/0.9.8-fwf66s4 libnova/0.15.0-hpvrrfr rocm/5.4.3
   fi
   
   # cmake is only required at build time, so we use the normal module load
   module load cmake/3.24.3
else 
   # any HPC other than Setonix :
   module purge
   module_load gcc/8.3.0 cascadelake  blink_test_data/devel cfitsio msfitslib pal/0.9.8 blink_astroio/master cuda/11.4.2

   # cmake is only required at build time, so we use the normal module load
   module load cmake/3.18.0
   
   export CXXFLAGS=-std=c++14
fi   
# build your software..
echo "Building the software.."

[ -d ${build_dir} ] || mkdir ${build_dir}
cd ${build_dir}
pwd
echo "cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_BUILD_TYPE=${build_opt} ${cmake_options}"
cmake .. -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} -DCMAKE_BUILD_TYPE=${build_opt} ${cmake_options}
make VERBOSE=1

# Install the software
make install

# test:
if [[ $dotests -gt 0 ]]; then
   echo "make test"
   make test
else
   echo "WARNING : tests are not required"
fi   

echo "Create the modulefile in $MODULEFILE_DIR (or $INSTALL_DIR)"
export ADDITIONAL_MODULEFILE_COMMANDS="prepend_path('BLINK_IMAGER_PATH', root_dir )"
create_modulefile

echo "Done."


