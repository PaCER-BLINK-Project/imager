#!/bin/bash

channel=175
if [[ -n "$1" && "$1" != "-" ]]; then
   channel=$1
fi

miriad_imaging=1
if [[ -n "$2" && "$2" != "-" ]]; then
   miriad_imaging=$2
fi

cal=cal
if [[ -n "$3" && "$3" != "-" ]]; then
   cal=$3
fi

for lfile_c in `ls *.LCCSPC`
do
   lfile=${lfile_c%%.LCCSPC}
   
   echo 
   echo "Processing lfile = $lfile started at:"
   date   
   echo "process_lfiles_eda2.sh ${channel} ${lfile} ${miriad_imaging} ${cal}"
   process_lfiles_eda2.sh ${channel} ${lfile} ${miriad_imaging} ${cal}
done

   
