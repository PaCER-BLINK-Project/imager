#!/bin/bash
#SBATCH --account=pawsey1154-gpu
#SBATCH --partition=gpu-dev
#SBATCH --exclusive
#SBATCH --time=01:00:00


module use  /software/projects/pawsey1154/setonix/2025.08/modules/zen3/gcc/14.2.0;
module use /software/projects/pawsey1154/cdipietrantonio/setonix/2025.08/modules/zen3/gcc/14.2.0;
module use /software/projects/pawsey1154/njayamanne/setonix/2025.08/modules/zen3/gcc/14.2.0;

module load blink-pipeline-gpu/main;
ldd `which blink_pipeline`
rm /scratch/pawsey1154/njayamanne/output_images/start_time_1592584200_coarse_133_real.fits
blink_pipeline -R 128 -c 4 -t 1s -L -o /scratch/pawsey1154/njayamanne/output_images -n 8192  -O 2 -M /scratch/pawsey1154/cdipietrantonio/1276619416/1276619416.metafits -r -s /scratch/pawsey1154/cdipietrantonio/1276619416/1276625432.bin -b 0 -I /scratch/pawsey1154/cdipietrantonio/1276619416/combined -X 0 -f -1 -P 276.6859632019899,-5.686834509583719 -A 71,21,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127 -Q 10 ;

echo "scp setonix:$(ls -1t /scratch/pawsey1154/njayamanne/1225713560_output_nathan/*_real.fits 2>/dev/null | head -n 1) ."

