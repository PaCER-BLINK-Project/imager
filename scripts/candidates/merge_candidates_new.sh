#!/bin/bash

# INFO : this is the proper script for merging as it also adds DATE/TIME as the last column

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
#SBATCH --account=pawsey0809
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16gb
#SBATCH --output=./merge_candidates.o%j
#SBATCH --error=./merge_candidates.e%j
#SBATCH --export=NONE

# find . -name "*.cand" -exec cat {} \; | awk '{if($1!="#" && $1>=10){print $0;}}' > candidates_snr10_or_higher.txt
template=20??????_??????_100ms_ch???
if [[ -n "$1" && "$1" != "-" ]]; then
   template="$1"   
fi

merged_file=all_merged_candidates.txt
if [[ -n "$2" && "$2" != "-" ]]; then
   merged_file="$2"
fi


first_dtm=`ls -d ${template} | head -1`
cat ${first_dtm}/fits_images/dynamic_spectra/????_????/*.cand | head -1 | awk '{print $0" DTM";}' > ${merged_file}

for dtm in `ls -d ${template}`
do
   cat ${dtm}/fits_images/dynamic_spectra/????_????/*.cand | awk -v dtm=${dtm} '{if($1!="#"){print $0" "dtm;}}' >> ${merged_file}
done
