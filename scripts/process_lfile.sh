#!/bin/bash

lfile_base=20230601_100213_100ms_ch294
if [[ -n "$1" && "$1" ]]; then
   lfile_base=$1
fi

imager_path=/software/projects/director2183/msok/imager_branches/main/imager/build_20230617/pacer_imager_multi  


echo "${imager_path} ${lfile_base} -f 229.6875 -n 180  -Z -B 1 -O fits_gpu -t 0.100 -a antenna_locations.txt"
${imager_path} ${lfile_base} -f 229.6875 -n 180  -Z -B 1 -O fits_gpu -t 0.100 -a antenna_locations.txt
