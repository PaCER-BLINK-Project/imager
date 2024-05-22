#!/bin/bash -l

#SBATCH --account=director2183-gpu # your account pawsey0809-gpu or director2183-gpu
#SBATCH --partition=gpu            # Using the gpu partition
#SBATCH --ntasks=1                 # Total number of tasks
#SBATCH --ntasks-per-node=1        # Set this for 1 mpi task per compute device
#SBATCH --cpus-per-task=1          # How many OpenMP threads per MPI task
#SBATCH --threads-per-core=1       # How many omp threads per core
#SBATCH --gres=gpu:1
#SBATCH --gpu-bind=closest         # Bind each MPI taks to the nearest GPU
#SBATCH --mem=64gb                 #Indicate the amount of memory per node when asking for share resources
#SBATCH --time=00:15:00
#SBATCH --output=./lfile_imager.o%j
#SBATCH --error=./lfile_imager.e%j

lfile=20230505_151304_test_100ms_ch124
if [[ -n "$1" && "$1" != "-" ]]; then
   lfile="$1"
fi

# load modules :
export IMAGER_DIR=/software/projects/director2183/msok/imager_branches/main_setonix_20230515/imager/
source ${IMAGER_DIR}/runtime.sh

echo "Started at :"
date

echo "${IMAGER_DIR}/build_gpu/pacer_imager_multi ${lfile} -f 94.531250 -n 180  -G -Z -B 1 -O fits_gpu -t 1.000000 > ${lfile}.out 2>&1"
${IMAGER_DIR}/build_gpu/pacer_imager_multi ${lfile} -f 94.531250 -n 180  -G -Z -B 1 -O fits_gpu -t 1.000000 > ${lfile}.out 2>&1

echo "Finished at:"
date

