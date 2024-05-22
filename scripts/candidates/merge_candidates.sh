#!/bin/bash

# BETTER TO USE : merge_candidates_new.sh
# INFO : this is an obsolate version of the script which makes a simple merge and selects only SNR>=10 candidates (does not include DATE/TIME column) 
#  

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
#SBATCH --account=pawsey0809
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=16gb
#SBATCH --output=./merge_candidates.o%j
#SBATCH --error=./merge_candidates.e%j
#SBATCH --export=NONE

find . -name "*.cand" -exec cat {} \; | awk '{if($1!="#" && $1>=10){print $0;}}' > candidates_snr10_or_higher.txt


