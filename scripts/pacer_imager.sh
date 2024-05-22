#!/bin/bash

re_fits=vis_real.fits
if [[ -n "$1" && "$1" != "-" ]]; then
   re_fits=$1
fi

im_fits=vis_imag.fits
if [[ -n "$2" && "$2" != "-" ]]; then
   im_fits=$2
fi

freq_mhz=246.09375
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_mhz=$3
fi

# -o miriad
extra_options=""
if [[ -n "$4" && "$4" != "-" ]]; then
   extra_options="$4"
fi

antenna_file=antenna_locations.txt
if [[ -n "$5" && "$5" != "-" ]]; then
   antenna_file=$5
fi

echo "#######################################################"
echo "# PARAMETERS :"
echo "#######################################################"
echo "Visibility files : $re_fits / $im_fits"
echo "Freq : $freq_mhz [MHz]"
echo "extra options = $extra_options"
echo "Antenna position file = $antenna_file"
echo "#######################################################"

if [[ ! -s $antenna_file ]]; then
   echo "ERROR : antenna position file $antenna_file not found -> exiting"
   exit -1
fi


export PATH=/home/msok/github/pacer/software/imager/pacer_dirty_image/build/:$PATH

ln -s $re_fits tmp_vis_real.fits
ln -s $im_fits tmp_vis_imag.fits

echo "pacer_dirty_imager tmp -f ${freq_mhz} -a ${antenna_file} -n 256 -w U -s ${extra_options}"
pacer_dirty_imager tmp -f ${freq_mhz} -a ${antenna_file} -n 256 -w U -s ${extra_options}

