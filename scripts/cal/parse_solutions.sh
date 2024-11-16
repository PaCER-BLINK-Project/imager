#!/bin/bash

first_coarse_channel=133
if [[ -n "$1" && "$1" != "-" ]]; then
   first_coarse_channel=$1
fi

fine_bw=0.04

ch=0
while [[ $ch -lt 768 ]];
do
   ch_str=`echo $ch | awk '{printf("%03d\n",$1);}'`   
   ch_str_float=`echo $ch | awk '{printf("%.4f\n",$1);}'`

   freq_mhz=`echo $ch | awk -v cc=${first_coarse_channel} -v fine_bw=${fine_bw} '{freq_mhz=cc*1.28-0.64+fine_bw/2.00+fine_bw*$1;printf("%.4f\n",freq_mhz);}'`
   coarse_channel=`echo $ch | awk -v cc=${first_coarse_channel} '{printf("%d\n",cc+($1/32));}'`
   
   echo "# frequency channel $coarse_channel = ${freq_mhz} MHz" 
   echo "# frequency channel $coarse_channel = ${freq_mhz} MHz"  > calsolutions_chan${ch_str}_xx.txt

   echo "awk -v ch_str_float=${ch_str_float} ..."
   awk -v ch_str_float=${ch_str_float} '{if($1==ch_str_float){print $4" "$2" "$3" "$2" "$3;}}' calsolutions_xx.txt  >> calsolutions_chan${ch_str}_xx.txt
   
   ch=$(($ch+1))
done
