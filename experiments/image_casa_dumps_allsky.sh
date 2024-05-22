#!/bin/bash

obsID=1392048016
if [[ -n "$1" && "$1" != "-" ]]; then
   obsID=$1
fi

start_coarse_channel=120
if [[ -n "$2" && "$2" != "-" ]]; then
   start_coarse_channel=$2
fi

base=1392048016_ch120
if [[ -n "$3" && "$3" != "-" ]]; then
   base=$3
fi

max_ch=32
if [[ -n "$4" && "$4" != "-" ]]; then
   max_ch=$4
fi

# size=480
size=1024
if [[ -n "$5" && "$5" != "-" ]]; then
   size=$5
fi

fov=180
if [[ -n "$6" && "$6" != "-" ]]; then
   fov=$6
fi

options=""
if [[ -n "$7" && "$7" != "-" ]]; then
   options="$7"
fi


# cd /media/msok/5508b34c-040a-4dce-a8ff-2c4510a5d1a3/mwa/data/1276619416/20240120/single_script_execution_24ch_SCRATCH/000000/

fine_ch=0
while [[ $fine_ch -lt $max_ch ]]; 
do
   echo 
   echo "------------------------------ channel = $fine_ch ------------------------------"
   fine_ch_str=`echo $fine_ch | awk '{printf("%03d\n",$1);}'`
   freq_mhz=`echo $fine_ch | awk -v start_coarse_channel=${start_coarse_channel} '{fine_ch=$1;freq_mhz=(start_coarse_channel*1.28 -0.64) + 0.04*fine_ch + 0.04/2.00;printf("%.8f\n",freq_mhz);}'`
   
   mkdir -p ${fine_ch_str}
#   echo "/usr/local/bin/pacer_dirty_imager_cpu ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -F ${fov} -w N -a antenna_locations.txt -M ${obsID}.metafits -O ${fine_ch_str}/ ${options} > ${fine_ch_str}/out 2>&1"
#   /usr/local/bin/pacer_dirty_imager_cpu ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -F ${fov} -w N -a antenna_locations.txt -M ${obsID}.metafits -O ${fine_ch_str}/ ${options} > ${fine_ch_str}/out 2>&1
   echo "/usr/local/bin/pacer_dirty_imager_cpu ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -a antenna_locations.txt  -w N -O ${fine_ch_str}/ -Z -v 100 -V 100 ${options} > ${fine_ch_str}/out 2>&1"
   /usr/local/bin/pacer_dirty_imager_cpu ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -a antenna_locations.txt  -w N -O ${fine_ch_str}/ -Z -v 100 -V 100 ${options} > ${fine_ch_str}/out 2>&1
   
   fine_ch=$(($fine_ch+1))
done

ls ???/dirty*_real.fits > fits_list

echo "avg_images fits_list avg.fits rms.fits -r 100000.000"
avg_images fits_list avg.fits rms.fits -r 100000.000 
