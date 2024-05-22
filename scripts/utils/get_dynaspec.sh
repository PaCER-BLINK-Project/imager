#!/bin/bash

event_path=./20230709_111828_100ms_ch294/fits_images/dynamic_spectra/0072_0122
if [[ -n "$1" && "$1" != "-" ]]; then
   event_path="$1"
fi

remote_dir=data-mover.pawsey.org.au:/scratch/director2183/msok/data/eda2/2023_07_09_eda2_10min_ch294/processing/all/
if [[ -n "$2" && "$2" != "-" ]]; then
   remote_dir="$2"
fi

local_dir="./"
if [[ -n "$3" && "$3" != "-" ]]; then
   local_dir=$3
   cd ${local_dir}
fi

subdir=`basename ${event_path}`

echo "rsync -avPL ${remote_dir}/${event_path} ."
rsync -avPL ${remote_dir}/${event_path} .
