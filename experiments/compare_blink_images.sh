#!/bin/bash

reference_path=/media/msok/5508b34c-040a-4dce-a8ff-2c4510a5d1a3/mwa/data/1276619416/FINAL/setonix/gpu/casa_dump/nimbus4_casa_dump/
if [[ -n "$1" && "$1" != "-" ]]; then
   reference_path="$1"
fi

for chdir in `ls -d ???perfreq8192x8192`
do
   echo "Comparing data in ${chdir} vs. ${reference_path}/${chdir}/ :"
   cd ${chdir}/
   for fits in `ls *.fits`; do    
#       echo "calcfits_bg $fits = ${reference_path}/${fits} | grep \"COMPARISON RESULT\""
       # dirty_image_20240404T081559000_real.fits
       fits_ref="${reference_path}/${chdir}/${fits}"
       fits_ref_prefix=`echo $fits | cut -b 1-11`
       if [[ $fits_ref_prefix == "dirty_image" ]]; then
          fits_ref_postfix=`echo $fits | cut -b 31-`
          echo "ls ${reference_path}/${chdir}/${fits_ref_prefix}*${fits_ref_postfix}"
          fits_ref=`ls ${reference_path}/${chdir}/${fits_ref_prefix}*${fits_ref_postfix} | head -1`          
       fi

       echo "calcfits_bg $fits = ${fits_ref}"
#       calcfits_bg $fits = ${fits_ref} | grep "COMPARISON RESULT for file $fits"
       calcfits_bg $fits = ${fits_ref} | grep "COMPARISON RESULT" 
       echo
       sleep 5
   done
   echo "----------------------------------------"
   echo
   echo
   sleep 2;
   cd -
done
