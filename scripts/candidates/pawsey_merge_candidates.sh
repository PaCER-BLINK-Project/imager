#!/bin/bash

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
#SBATCH --account=pawsey0809
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=32gb
#SBATCH --output=./pawsey_merge_candidates.o%j
#SBATCH --error=./pawsey_merge_candidates.e%j
#SBATCH --export=NONE

# find . -name "*.cand" -exec cat {} \; | awk '{if($1!="#" && $1>=10){print $0;}}' > candidates_snr10_or_higher.txt

threshold=10
if [[ -n "$1" && "$1" != "-" ]]; then
   threshold=$1
fi

outfile=candidates_threshold_${threshold}sigma.txt
if [[ -n "$2" && "$2" != "-" ]]; then
   outfile="$2"
fi

echo "###################################"
echo "PARAMETERS:"
echo "###################################"
echo "threshold = $threshold"
echo "outfile   = $outfile"
echo "###################################"


find . -name "*.cand" > cand_files.txt

first_file=`head -1 cand_files.txt`
head -1 ${first_file} > ${outfile}

for candfile in `cat cand_files.txt`
do
   echo "cat ${candfile} | awk -v file=${candfile} -v t=${threshold} '{if(\$1!=\"#\" && \$1>=t){print file\" : \"$0;}}' >> ${outfile}"
   cat ${candfile} | awk -v file=${candfile} -v t=${threshold} '{if($1!="#" && $1>=t){print file" : "$0;}}' >> ${outfile}   
done
