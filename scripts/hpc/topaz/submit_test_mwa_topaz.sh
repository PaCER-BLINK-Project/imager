#!/bin/bash

# salloc --partition gpuq --time 1:00:00 --nodes=1

module purge
module load cascadelake slurm/20.02.3 gcc/8.3.0 cmake/3.18.0
# module use /group/director2183/software/centos7.6/modulefiles
module use /group/director2183/msok/software/centos7.6/modulefiles
module load msfitslib/devel 
module load pacer_dirty_image/devel 

# 100 times (-P 100) :
/group/director2183/msok/imager/pacer_dirty_image/build_release/pacer_dirty_imager 1103645160 -p _channel000_time000000 -n 2049 -f 138.88 -F 30 -m 30 -w U -a antenna_locations.txt -M 1103645160.txt -v -1000 -P 100 -V -1
