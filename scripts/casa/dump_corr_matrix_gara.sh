#!/bin/bash -l

#SBATCH --account=mwavcs
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64gb
#SBATCH --output=./dump_corr_matrix_gara.o%j
#SBATCH --error=./dump_corr_matrix_gara.e%j
#SBATCH --export=NONE

module load casa/5.6.1-8
export BLINK_DIR=/software/projects/mwavcs/msok/imager/

casa_ms=1276619416_20200619163000.ms
if [[ -n "$1" && "$1" != "-" ]]; then
   casa_ms=$1
fi

# When dumping data from SMART pipeline there is no CORRECTED_DATA as calibration is applied on-the-fly by cotter and saved to column DATA :
echo "/software/projects/mwavcs/msok/imager/scripts/casa/get_corr_matrix.sh $casa_ms 768 1 DATA"
/software/projects/mwavcs/msok/imager/scripts/casa/get_corr_matrix.sh $casa_ms 768 1 DATA
