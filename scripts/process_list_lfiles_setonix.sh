#!/bin/bash

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
# -account=pawsey0809
#SBATCH --account=director2183
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem=64gb
#SBATCH --output=./process_all_lfiles_setonix.o%j
#SBATCH --error=./process_all_lfiles_setonix.e%j
#SBATCH --export=NONE

lafile_list=lafiles.txt
if [[ -n "$1" && "$1" != "-" ]]; then
   lafile_list=$1
fi

datadir=/scratch/director2183/msok/data/eda2/2023_06_01_eda2_10min_ch294/
if [[ -n "$2" && "$2" != "-" ]]; then
   datadir="$2"
fi

compress=1
if [[ -n "$3" && "$3" != "-" ]]; then
   compress=$3
fi

force=0
if [[ -n "$4" && "$4" != "-" ]]; then
   force=$4
fi

# center_freq_mhz=229.6875 # ch*(400/512)
freq_channel=294
if [[ -n "$5" && "$5" != "-" ]]; then
   # center_freq_mhz=$4
   freq_channel=$5
fi
center_freq_mhz=`echo $freq_channel | awk '{printf("%.8f\n",$1*(400.00/512.00));}'`

run_frb_search=0
if [[ -n "$6" && "$6" != "-" ]]; then
   run_frb_search=$6
fi

pacer_imager_path=/software/projects/director2183/msok/imager_branches/main/imager/build/pacer_imager_multi
if [[ -n "$7" && "$7" != "-" ]]; then
   pacer_imager_path="$7"
fi

modules="blink-imager-cpu/devel"
if [[ -n "$8" && "$8" != "-" ]]; then
   modules="$8"
fi

imager_options=""
if [[ -n "$9" && "$9" != "-" ]]; then
   imager_options="$9"
fi

echo "#####################################"
echo "PARAMETERS:"
echo "#####################################"
echo "datadir  = $datadir"
echo "lafile_list = $lafile_list"
echo "compress = $compress"
echo "force    = $force"
echo "freq_ch  = $freq_channel -> center_freq_mhz = $center_freq_mhz"
echo "run_frb_search = $run_frb_search"
echo "pacer_imager_path = $pacer_imager_path"
echo "modules = $modules"
echo "imager_options = $imager_options"
echo "#####################################"

# 2024-01-09 : may have to be removed:
module reset
module use /software/projects/director2183/setonix/modules/zen3/gcc/12.1.0  /software/projects/director2183/msok/setonix/modules

echo "module load $modules"
module load $modules

# /software/projects/director2183/msok/imager_branches/main/imager/build_20230617/pacer_imager_multi  20230601_100213_100ms_ch294 -f 229.6875 -n 180  -Z -B 1 -O fits_gpu -t 0.100 -C 32 -j -1

for lfile_a in `cat ${lafile_list}`
do
   lfile_c=${lfile_a%%LACSPC}LCCSPC
   lfile_base_name=`basename $lfile_a`
   l_file_base=${lfile_base_name%%.LACSPC}   
   
   tar_file=${l_file_base}.tar.gz
   
   if [[ ! -s ${tar_file} || $force -gt 0 ]]; then   
     echo 
     echo "-------------------------- Processing L-file $lfile_base --------------------------"
     mkdir -p ${l_file_base}
     cd ${l_file_base}
     pwd
   
     echo "ln -s ${datadir}/${lfile_a}"
     ln -s ${datadir}/${lfile_a}
   
     echo "ln -s ${datadir}/${lfile_c}"
     ln -s ${datadir}/${lfile_c}

     if [[ -s antenna_locations.txt ]]; then   
        echo "INFO : file antenna_locations.txt already exists"
     else
        echo "cp $BLINK_IMAGER_PATH/config/eda2/antenna_locations.txt ."
        cp $BLINK_IMAGER_PATH/config/eda2/antenna_locations.txt .
     fi

     if [[ -s instr_config.txt ]]; then   
        echo "INFO : file instr_config.txt already exists"
     else
        echo "cp $BLINK_IMAGER_PATH/config/eda2/instr_config.txt ."
        cp $BLINK_IMAGER_PATH/config/eda2/instr_config.txt .
     fi
   
     echo "time ${pacer_imager_path} ${l_file_base} -f ${center_freq_mhz} -n 180  -Z -B 1 -O fits_gpu -t 0.100 -C 32 -j -1 $imager_options > out 2>&1"
     time ${pacer_imager_path} ${l_file_base} -f ${center_freq_mhz} -n 180  -Z -B 1 -O fits_gpu -t 0.100 -C 32 -j -1 $imager_options > out 2>&1 
   
     cd ..      
     
     if [[ $run_frb_search -gt 0 ]]; then
        echo "INFO : frb_search is required executing ..."
        
        # TODO : create frb_search project/module and use it here, rather than script from by home directory:
        # no sbatch here as we want to do it here and now, before compression takes place!
        echo "~/bin/process_frb_search.sh $l_file_base 1 ${freq_channel}"
        ~/bin/process_frb_search.sh $l_file_base 1 ${freq_channel}
     else
        echo "WARNING : frb_search is not required"
     fi
     
     if [[ $compress -gt 0 ]]; then
        echo "time tar zcvf ${tar_file} ${l_file_base}"
        time tar zcvf ${tar_file} ${l_file_base}
        
        echo "time rm -fr ${l_file_base}/"
        time rm -fr ${l_file_base}/
     fi
  else
     echo "INFO : L-files $lfile_a and $lfile_c already processed ( file $tar_file found )"  
  fi
done
