#!/bin/bash

casa_ms=1276619416_20200619163000.ms
if [[ -n "$1" && "$1" != "-" ]]; then
   casa_ms="$1"
fi

max_ch=768
if [[ -n "$2" && "$2" != "-" ]]; then
   max_ch=$2
fi

max_time=100
if [[ -n "$3" && "$3" != "-" ]]; then
   max_time=$3
fi

column=CORRECTED_DATA
if [[ -n "$4" && "$4" != "-" ]]; then
   column=$4
fi


base_casa_ms=${casa_ms%%.ms}

pol=0
timeidx=0

while [[ $timeidx -lt $max_time ]];
do
   timeidx_str=`echo $timeidx | awk '{printf("%06d",$1);}'`

   mkdir -p ${timeidx_str}
   cd ${timeidx_str}
   ln -s ../${casa_ms}
   ch=0
   while [[ $ch -lt $max_ch ]];
   do
      ch_str=`echo $ch | awk '{printf("%03d",$1);}'`
      corr_re="${base_casa_ms}_vis_real_channel${ch_str}_time${timeidx_str}_pol${pol}.fits"       
       
      if [[ -s ${corr_re} ]]; then
         echo "INFO : FITS file ${corr_re} already exists -> skipped"
      else
         # was ~/github/pacer/software/imager/
         echo "casapy --nologger -c $BLINK_DIR/scripts/casa/get_corrmatrix_from_casa.py ${casa_ms} --data_column=${column} --save_full_matrix --channel=${ch} --time_index=${timeidx} --pol=${pol}"
         casapy --nologger -c $BLINK_DIR/scripts/casa/get_corrmatrix_from_casa.py ${casa_ms} --data_column=${column} --save_full_matrix --channel=${ch} --time_index=${timeidx} --pol=${pol}
      fi
#      sleep 1

      ch=$(($ch+1))
   done
   cd ..

   timeidx=$(($timeidx+1))
done
