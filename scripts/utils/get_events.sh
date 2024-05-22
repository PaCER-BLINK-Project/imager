#!/bin/bash

export PATH=~/github/mwafrb/scripts:$PATH

remote_dir=data-mover.pawsey.org.au:/scratch/director2183/msok/data/eda2/2023_06_01_eda2_10min_ch294/processing/cpu_all/flipped_as_on_laptop_REPEAT_NEW_CODE
# remote_dir=setonix.pawsey.org.au:/scratch/director2183/msok/data/eda2/2023_06_01_eda2_10min_ch294/processing/cpu_all/flipped_as_on_laptop_REPEAT_NEW_CODE
if [[ -n "$1" && "$1" != "-" ]]; then
   remote_dir="$1"
fi

cand_file=candidates_threshold_15sigma_nonzeroDM.txt
if [[ -n "$2" && "$2" != "-" ]]; then
   cand_file="$2"
fi

view=1
if [[ -n "$3" && "$3" != "-" ]]; then
   view=$3
fi

copy=1
if [[ -n "$4" && "$4" != "-" ]]; then
   copy=$4
fi

copy_images=0
if [[ -n "$5" && "$5" != "-" ]]; then
   copy_images=$5
fi

evt_name=
if [[ -n "$6" && "$6" != "-" ]]; then
   evt_name="$6"
fi

show_images=0
if [[ -n "$7" && "$7" != "-" ]]; then
   show_images=$7
fi


echo "##############################################"
echo "PARAMETERS:"
echo "##############################################"
echo "cand_file  = $cand_file"
echo "remote_dir = $remote_dir"
echo "view       = $view"
echo "copy       = $copy"
echo "copy_images = $copy_images"
echo "evt_name   = $evt_name"
echo "show_images = $show_images"
echo "##############################################"

while read line # example 
do
   path=`echo $line | awk '{print $1;}'`
   b=`basename $path`
   d=`dirname $path`
   evt=`echo $b | cut -b 1-9`
   dtm=`echo $d | awk -F "/" '{print $2;}'`
   echo "DEBUG : $d -> dtm = $dtm"
   
   if [[ -n "$evt_name" ]]; then
      if [[ $evt_name != $evt ]]; then
         echo "Event $evt != required $evt_name -> skipped"
         continue;
      fi
   fi

   echo "$path -> $b -> $evt"
   
   mkdir -p ${dtm}/${evt}
   cd ${dtm}/
   
   candfile=${evt}_series.cand
   fitsfile=${candfile%%cand}fits # use this line to show dynamic spectrum
   fitsfile=${evt}.fits
 
   if [[ $copy -gt 0 ]]; then
      if [[ ! -s ${evt}/${evt}.fits ]]; then      
         echo "rsync -avPL ${remote_dir}/${d}/ ${evt}/"
         rsync -avPL ${remote_dir}/${d}/* ${evt}/
      else
         echo "Event ${evt} already copied"
      fi

      if [[ $copy_images -gt 0 ]]; then
         # mkdir -p ${evt}/fits_images
         # cd ${evt}/fits_images
         mkdir -p fits_images
         cd fits_images

         for ch in {0..31}
         do
            ch_str=`echo $ch | awk '{printf("%05d",$1)}'`
            mkdir -p ${ch_str}
            cd ${ch_str}
            if [[ ! -s fits_list ]]; then
               echo "rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/fits_list ."
               rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/fits_list .
            fi         
         
            pwd
            for evt_t in `cat ../../${evt}/*cand | awk '{if($1!="#"){print $2;}}'`
            do
               list=""
               for fits in `cat fits_list | awk -v evt_t=${evt_t} '{if(NR>=(evt_t+1-1) && NR<=(evt_t+1+2)){print $1;}}'`
               do            
                  # echo "rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/${fits} ."
                  # rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/${fits} .               
                  list="${fits},${list}"
                  # sleep 1
               done
               # list="${list}${fits}"
               list2=`echo $list | awk '{l=length($1);print substr($1,1,l-1);}'`
               echo "rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/{${list2}} ."
               echo "rsync -avP ${remote_dir}/${dtm}/fits_images/${ch_str}/{${list2}} ." > scp!
               chmod +x scp!
               ./scp!
               sleep 1 # otherwise "ssh connection closed"
            done
            cd ../
         done
         cd ..
         
         # create symbolic links
         echo "INFO : creating symbolic links"
         mkdir -p ${evt}/fits_images
         cd ${evt}/fits_images
         for ch in {0..31}
         do
            ch_str=`echo $ch | awk '{printf("%05d",$1)}'`
            cd ${ch_str}
            for fits_file in `ls ../../fits_images/${ch_str}/*.fits`
            do
               echo "ln -s ${fits_file}" 
               ln -s ${fits_file}
            done
            cd ..
         done
         cd ../../
      else
         echo "INFO : copy of images is not required"
      fi
   else
      echo "WARNING : copying is not required"
   fi

   if [[ $view -gt 0 ]]; then
      cd ${evt}
      echo "~/github/mwafrb/scripts/ds9_fredda_all! 6 100000 1 ${candfile} ${fitsfile} -1 6 10"
      ~/github/mwafrb/scripts/ds9_fredda_all! 6 100000 1 ${candfile} ${fitsfile} -1 6 10

      if [[ $show_images -gt 0 ]]; then
         cd fits_images      
         for ch in {0..31}
         do
            ch_str=`echo $ch | awk '{printf("%05d",$1)}'`
            for evt_t in `cat ../*cand | awk '{if($1!="#"){print $2;}}'`
            do
               fits=`cat ${ch_str}/fits_list | awk -v evt_t=${evt_t} '{if(NR==(evt_t+1)){print $1;}}'`
          
               echo "ds9 -scale zscale ${ch_str}$fits &"
               ds9 -scale zscale ${ch_str}/$fits &
               sleep 5
            done
         done
         cd ..
      fi
      cd ..
   fi      
   cd ..
done < $cand_file
