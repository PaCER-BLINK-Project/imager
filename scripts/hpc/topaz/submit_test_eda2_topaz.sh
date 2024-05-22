#!/bin/bash

# salloc --partition gpuq --time 1:00:00 --nodes=1

module purge
module load cascadelake slurm/20.02.3 gcc/8.3.0 cmake/3.18.0
# module use /group/director2183/software/centos7.6/modulefiles
module use /group/director2183/msok/software/centos7.6/modulefiles
module load msfitslib/devel 
module load pacer_dirty_image/devel 

# pacer_dirty_imager 20191104_033537_eda2_ch1_ant256_midday_avg1 -p _channel000_time000000 -n 2049
# pacer_dirty_imager 20191104_033537_eda2_ch1_ant256_midday_avg1 -p _channel000_time000000 -n 512 -a antenna_locations.txt
/group/director2183/msok/imager/pacer_dirty_image/build_release/pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 180 -w N -o miriad  -v -1 -P 100 -V -1 -Z


