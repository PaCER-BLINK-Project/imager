#!/bin/bash
#SBATCH --account=pawsey1154-gpu
#SBATCH --partition=gpu-dev
#SBATCH --exclusive
#SBATCH --time=01:00:00

module use  /software/projects/pawsey1154/setonix/2025.08/modules/zen3/gcc/14.2.0;
module use /software/projects/pawsey1154/cdipietrantonio/setonix/2025.08/modules/zen3/gcc/14.2.0;
module use /software/projects/pawsey1154/njayamanne/setonix/2025.08/modules/zen3/gcc/14.2.0;

module load blink-pipeline-gpu/main;
if [ -e /scratch/pawsey1154/njayamanne/1225713560_output_nathan ]; then rm -r /scratch/pawsey1154/njayamanne/1225713560_output_nathan ; fi
blink_pipeline -R 128 -c 4 -t 1s -L -o /scratch/pawsey1154/njayamanne/1225713560_output_nathan -n 512 -O 2 -M /scratch/pawsey1154/cdipietrantonio/1225713560/1225713560.metafits -r -s /scratch/pawsey1154/cdipietrantonio/1225713560/1225736152.bin -b 0 -I /scratch/pawsey1154/cdipietrantonio/1225713560/combined -X 0 -f -1 -Q 5 -P 353.834580131841,-40.52218992025596 -A 48,66,87,25,108,46 ;


echo "scp setonix:$(ls -1t /scratch/pawsey1154/njayamanne/1225713560_output_nathan/*_real.fits 2>/dev/null | head -n 1) ~/Downloads/"





