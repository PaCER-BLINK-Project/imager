#!/bin/bash

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --time=23:59:00
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --output=./image_casa_dumps_perfreq_gpu.o%j
#SBATCH --error=./image_casa_dumps_perfreq_gpu.e%j
#SBATCH --export=NONE

echo "Load modules:"
module load blink-imager-gpu/gpu
module load msfitslib/devel
echo "Modules loaded"


obsID=1276619416
if [[ -n "$1" && "$1" != "-" ]]; then
   obsID=$1
fi

start_coarse_channel=133
if [[ -n "$2" && "$2" != "-" ]]; then
   start_coarse_channel=$2
fi

base=1276619416_20200619163000
if [[ -n "$3" && "$3" != "-" ]]; then
   base=$3
fi

max_ch=768
if [[ -n "$4" && "$4" != "-" ]]; then
   max_ch=$4
fi

size=4096
if [[ -n "$5" && "$5" != "-" ]]; then
   size=$5
fi

fov=16.72
if [[ -n "$6" && "$6" != "-" ]]; then
   fov=$6
fi

options=""
if [[ -n "$7" && "$7" != "-" ]]; then
   options="$7"
fi

postfix=""
if [[ -n "$8" && "$8" != "-" ]]; then
   postfix="$8"
fi

imager_path=/software/projects/director2183/msok/setonix/2023.08/development/blink-imager-gpu/gpu/bin/pacer_dirty_imager_gpu
if [[ -n "$9" && "$9" != "-" ]]; then
   imager_path="$9"
fi

force=0
fine_ch=0
while [[ $fine_ch -lt $max_ch ]]; 
do
   # freq_mhz=`echo $center_coarse_channel | awk '{print $1*1.28;}'`
   freq_mhz=`echo $fine_ch | awk -v start_coarse_channel=${start_coarse_channel} '{fine_ch=$1;freq_mhz=(start_coarse_channel*1.28 -0.64) + 0.04*fine_ch + 0.04/2.00;printf("%.8f\n",freq_mhz);}'`
   fine_ch_str=`echo $fine_ch | awk '{printf("%03d\n",$1);}'`
   echo 
   echo "------------------------------ channel = $fine_ch (freq_mhz = $freq_mhz - using centre freq. ) ------------------------------"
   
   count=`ls ${fine_ch_str}${postfix}/dirty_image*_real.fits | wc -l`      
   if [[ $count -le 0 || $force -gt 0 ]]; then     
      mkdir -p ${fine_ch_str}${postfix}
      echo "${imager_path} ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -F ${fov} -w N -M ${obsID}.metafits -O ${fine_ch_str}${postfix}/ ${options} > ${fine_ch_str}${postfix}/out 2>&1"
      ${imager_path} ${base} -p _channel${fine_ch_str}_time000000_pol0 -n ${size} -f ${freq_mhz} -F ${fov} -w N -M ${obsID}.metafits -O ${fine_ch_str}${postfix}/ ${options} > ${fine_ch_str}${postfix}/out 2>&1
   else
      echo "WARNING : already imaged -> skipped"
   fi
   
   fine_ch=$(($fine_ch+1))
done

ls ???${postfix}/dirty*_real.fits > fits_list

echo "avg_images fits_list avg.fits rms.fits -r 100000.000"
avg_images fits_list avg.fits rms.fits -r 100000.000 
