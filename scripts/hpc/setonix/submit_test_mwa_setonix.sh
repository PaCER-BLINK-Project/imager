#!/bin/bash

# salloc --partition gpuq --time 1:00:00 --nodes=1

# Note : on setonix nothing is loaded by default so module purge is not needed
module use /software/projects/director2183/msok/setonix/modules/

module load cfitsio/4.0.0
module load libnova/0.15.0-pfj3ef7
module load fftw/3.3.9
module load msfitslib/devel

# 100 times (-P 100) :
/software/projects/director2183/msok/imager/pacer_dirty_image/build_release/pacer_dirty_imager 1103645160 -p _channel000_time000000 -n 2049 -f 138.88 -F 30 -m 30 -w U -a antenna_locations.txt -M 1103645160.txt -v -1000 -P 100 -V -1
