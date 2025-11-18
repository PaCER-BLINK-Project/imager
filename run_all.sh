#!/bin/bash -e
# Run the following line before executing this script, to get a GPU node.
# salloc -p gpu-dev -A pawsey1154-gpu --gres=gpu:1  

./build_imager.sh gpu
cd build_gpu
rm -rf test_imager_conv* 
module load blink_test_data/devel cfitsio/4.4.0  msfitslib/master-ittkjmq blink-astroio/master fftw/3.3.10  pal/0.9.8-n3thcaw  libnova/0.15.0-iwh6cpn rocm/6.4.1
./conv_test

GPU_VERSION=`realpath test_imager_conv_gpu/start_time_1508442485_coarse_109_real.fits`
CPU_VERSION=`realpath test_imager_conv_cpu/start_time_1508442485_coarse_109_real.fits`

echo "scp setonix:$GPU_VERSION ./gpu_version.fits"
echo "scp setonix:$CPU_VERSION ./cpu_version.fits"


