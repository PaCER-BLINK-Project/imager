#!/bin/bash

datadir=/data_archive/2023_05_28_ch94_1second/
if [[ -n "$1" && "$1" != "-" ]]; then
   datadir="$1"
fi

for lfile_c in `ls ${datadir}/*.LCCSPC`
do
   lfile=${lfile_c%%.LCCSPC}
   lfile_a=${lfile_c%%LCCSPC}LACSPC
   
   echo "ln -s ${lfile_c}"
   ln -s ${lfile_c}

   echo "ln -s ${lfile_a}"
   ln -s ${lfile_a}
done
