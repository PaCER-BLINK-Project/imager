#!/bin/bash

uvfits=chan_315_20220319T04240130.uvfits
if [[ -n "$1" && "$1" != "-" ]]; then
   uvfits=$1
fi
base=${uvfits%%.uvfits}

n_chan=32
if [[ -n "$2" && "$2" != "-" ]]; then
   n_chan=$2
fi

freq_channel=315
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_channel=$3
fi

freq_mhz=`echo $freq_channel | awk '{printf("%.8f\n",($1*(400.00/512.00)));}'`
bw=`echo 1 | awk '{printf("%.8f\n",(400.00/512.00)*(32.00/27.00));}'`
channel_bw=`echo "$bw $n_chan" | awk '{printf("%.8f\n",($1/$2));}'`
half_channel_bw=`echo $channel_bw | awk '{printf("%.8f\n",($1/2.00));}'`
freq_start=`echo "$freq_mhz $bw" | awk '{printf("%.8f\n",($1-$2/2.00));}'`
freq_end=`echo "$freq_mhz $bw" | awk '{printf("%.8f\n",($1+$2/2.00));}'`

echo "Frequency channel $freq_channel -> $freq_mhz [MHz] ( $freq_start - $freq_end [MHz] )"
echo "Bandwidth : $bw [MHz] , fine channel bw = $channel_bw [MHz] (half channel bw = $half_channel_bw ) number of channels = $n_chan"

ch=0
while [[ $ch -lt $n_chan ]]; 
do
   ch_str=`echo $ch | awk '{printf("%03d",$1);}'`
   mkdir -p ${ch_str}
   cd ${ch_str}
   ln -s ../${uvfits}
   ln -s ../antenna_locations.txt .

   echo "python ${BLINK_DIR}/imager/scripts/eda2_calibration.py $uvfits - --channel=${ch} --dump"
   python ${BLINK_DIR}/imager/scripts/eda2_calibration.py $uvfits - --channel=${ch} --dump

   echo "ln -sf ${base}_CorrMatrix_RE.fits ${base}_vis_real.fits"
   ln -sf ${base}_CorrMatrix_RE.fits ${base}_vis_real.fits   
   
   echo "ln -sf ${base}_CorrMatrix_IM.fits ${base}_vis_imag.fits"
   ln -sf ${base}_CorrMatrix_IM.fits ${base}_vis_imag.fits
   
   ch_freq_mhz=`echo "$freq_start $channel_bw $ch" | awk '{printf("%.8f\n",$1+$2*$3+$2/2.00);}'`
   
   # 246.09375
   echo "/home/msok/github/pacer/software/imager/pacer_dirty_image/build/pacer_dirty_imager $base -f ${ch_freq_mhz} -a antenna_locations.txt -n 256 -w N -s -o miriad"
   /home/msok/github/pacer/software/imager/pacer_dirty_image/build/pacer_dirty_imager $base -f ${ch_freq_mhz} -a antenna_locations.txt -n 256 -w N -s -o miriad
   
   cd ..
   
   ch=$(($ch+1))
done

