#!/bin/bash

# salloc --partition gpuq --time 1:00:00 --nodes=1

# Note : on setonix nothing is loaded by default so module purge is not needed
module use /software/projects/director2183/msok/setonix/modules/

module load cfitsio/4.0.0
module load libnova/0.15.0-pfj3ef7
module load fftw/3.3.9
module load msfitslib/devel

# pacer_dirty_imager 20191104_033537_eda2_ch1_ant256_midday_avg1 -p _channel000_time000000 -n 2049
# pacer_dirty_imager 20191104_033537_eda2_ch1_ant256_midday_avg1 -p _channel000_time000000 -n 512 -a antenna_locations.txt
/software/projects/director2183/msok/imager/pacer_dirty_image/build_release/pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 180 -w N -o miriad  -v -1 -P 100 -V -1 -Z


