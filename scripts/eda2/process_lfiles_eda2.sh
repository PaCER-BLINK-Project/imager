#!/bin/bash

channel=175
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

lfile=20230418_130239_test_100ms_ch124
if [[ -n "$2" && "$2" != "-" ]]; then
   lfile="$2"
fi

miriad_imaging=1
if [[ -n "$3" && "$3" != "-" ]]; then
   miriad_imaging=$3
fi

cal=cal
if [[ -n "$4" && "$4" != "-" ]]; then
   cal=$4
fi

timeres=1.0001664 # in seconds
if [[ -n "$5" && "$5" != "-" ]]; then
   timeres=$5
fi

image_all_channels=1
if [[ -n "$6" && "$6" != "-" ]]; then
   image_all_channels=$6
fi


echo "##############################################"
echo "PARAMETERS:"
echo "##############################################"
echo "timeres = $timeres"
echo "image_all_channels = $image_all_channels"
echo "##############################################"

echo "~/aavs-calibration/Lfile2uvfits_eda.sh -i ${timeres} -n 1 -N 512 -C 32 -f ${channel}  -F -s eda2 ${lfile} > lfile2uvfits.out 2>&1"
~/aavs-calibration/Lfile2uvfits_eda.sh -i ${timeres} -n 1 -N 512 -C 32 -f ${channel}  -F -s eda2 ${lfile} > lfile2uvfits.out 2>&1

if [[ $miriad_imaging -gt 0 ]]; then
   ls *.uvfits > uvfits_list   
   echo "miriad_applycal_and_image_list.sh uvfits_list ${cal} 256 - - - - - 1 ${image_all_channels} > uv2image.out 2>&1"
   miriad_applycal_and_image_list.sh uvfits_list ${cal} 256 - - - - - 1 ${image_all_channels} > uv2image.out 2>&1
   xy2i.sh > xy2i.out 2>&1 

   ./remove!
else
   echo "WARNING : MIRIAD imaging is not required"
fi

# cd images_i/
