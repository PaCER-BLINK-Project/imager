#!/bin/bash

lfiles_per_job=10
if [[ -n "$1" && "$1" != "-" ]]; then
   lfiles_per_job=$1
fi

curr_path=`pwd`
options=" 1 - 294 1"
if [[ -n "$2" && "$2" != "-" ]]; then
   options="$2"
fi

echo "################################"
echo "PARAMETERS:"
echo "################################"
echo "lfiles_per_job = $lfiles_per_job"
echo "options        = $options"
echo "################################"


ls *.LACSPC > all_lafiles_list.txt
split --additional-suffix=.txt --lines=10 all_lafiles_list.txt lafiles_list

echo "Commands to run:"
for lafile in `ls lafiles_list*.txt`
do
   echo "sbatch /software/projects/director2183/msok/imager_branches/main/imager/scripts/process_list_lfiles_setonix.sh ${lafile} ${curr_path} ${options}"
done
