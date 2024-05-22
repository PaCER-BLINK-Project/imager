#!/bin/bash

corr_basename=corr_matrix
if [[ -n "$1" && "$1" != "-" ]]; then
   corr_basename=$1
fi

freq_ch=315
if [[ -n "$2" && "$2" != "-" ]]; then
   freq_ch=$2
fi

# -o miriad
extra_options=""
if [[ -n "$3" && "$3" != "-" ]]; then
   extra_options="$3"
fi

antenna_file=~/aavs-calibration/config/eda2/antenna_locations.txt
if [[ -n "$4" && "$4" != "-" ]]; then
   antenna_file=$4
fi

fits2jpg=0
if [[ -n "$5" && "$5" != "-" ]]; then
   fits2jpg=$5
fi


echo "#######################################################"
echo "# PARAMETERS :"
echo "#######################################################"
echo "Visibility files basename : $corr_basename"
echo "Freq channel : $freq_ch"
echo "extra options = $extra_options"
echo "Antenna position file = $antenna_file"
echo "Convert FITS to JPG : $fits2jpg"
echo "#######################################################"

if [[ ! -s $antenna_file ]]; then
   echo "ERROR : antenna position file $antenna_file not found -> exiting"
   exit -1
fi


export PATH=/home/msok/github/pacer/software/imager/pacer_dirty_image/build/:/home/msok/github/pacer/software/imager/pacer_dirty_image/scripts/:$PATH

for timestamp in `ls -d ????????`
do
   cd ${timestamp}
   echo
   echo
   echo "------------------------------ Time = $timestamp ------------------------------"
   for ch in `ls -d ???`
   do
      cd $ch
      pwd
   
      freq_mhz=`echo "$freq_ch $ch" | awk '{freq_ch=$1;ch=$2;fine_freq_mhz=freq_ch*(400.00/512.00)-((400.00/512.00)*(32.00/27.00))/2 + (((400.00/512.00)*(32.00/27.00))/32.00)/2.00;printf("%.8f\n",fine_freq_mhz);}'`
      echo "------------------------------ channel $ch / frequency $freq_mhz [MHz] ------------------------------"
      date
      echo "ln -s ${antenna_file}"
      ln -s ${antenna_file}
   
      file_base=${corr_basename}_ch${ch}
   
      echo "pacer_dirty_imager $file_base -f ${freq_mhz} -a ${antenna_file} -n 256 -w U -s ${extra_options}"
      pacer_dirty_imager $file_base -f ${freq_mhz} -a ${antenna_file} -n 256 -w U -s ${extra_options}

      if [[ $fits2jpg -gt 0 ]]; then   
         path=`which imagefits.py`
         echo "python $path dirty_test_real_fftshift_256x256.fits --ext=png"
         python $path dirty_test_real_fftshift_256x256.fits --ext=png
      fi
 
# WARNING : make sure the link is removed not the target file !!!     
#      echo "rm -f ${antenna_file}"
#      rm -f ${antenna_file}
   
      cd ../
      pwd
      date
      echo "-----------------------------------------------------------------------------------------------------"
   done
   
   cd ..
   
done
   

