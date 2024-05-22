#!/bin/bash

rm -f dirty*.fits 

# was -f 138.88 -> 138.90 
# /home/msok/github/pacer/software/imager/build/pacer_dirty_imager 1103645160 -p _channel000_time000000_pol0 -n 456 -f 138.90 -F 16.72 -w N -M 1103645160.metafits 

ch=0
while [[ $ch -lt 768 ]];
do
  freq=`echo $ch | awk '{ch=$1;freq=1.28*109-1.28/2.00+(ch+0.5)*(1.28/32.00);printf("%.8f\n",freq);}'`
  ch_str=`echo $ch | awk '{printf("%03d",$1);}'`
  echo 
  echo "ch = $ch -> freq = $freq"
  echo "/home/msok/github/pacer/software/imager/build/pacer_dirty_imager 1103645160 -p _channel${ch_str}_time000000_pol0 -n 456 -f ${freq} -F 16.72 -w U -M 1103645160.metafits"
  /home/msok/github/pacer/software/imager/build/pacer_dirty_imager 1103645160 -p _channel${ch_str}_time000000_pol0 -n 456 -f ${freq} -F 16.72 -w U -M 1103645160.metafits
  
  ch=$(($ch+1))
done


ls dirty*_real.fits > fits_list 
echo "avg_images fits_list"
avg_images fits_list
