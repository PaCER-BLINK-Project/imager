#!/bin/bash

# SETONIX : --account=director2183 - use explicit option of sbatch vs. 
#SBATCH --account=pawsey0809
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --mem=8gb
#SBATCH --output=./pawsey_rsync.o%j
#SBATCH --error=./pawsey_rsync.e%j
#SBATCH --export=NONE

# module load msfitslib/devel

# load python virtual env :
echo "source /scratch/director2183/msok/data/eda2/scripts/py-virtual-env/bin/activate"
source /scratch/director2183/msok/data/eda2/scripts/py-virtual-env/bin/activate 

# 20230601_100213_100ms_ch294/fits_images/00000/
template=20230601_??????_100ms_ch294
if [[ -n "$1" && "$1" != "-" ]]; then
   template="$1"
fi

echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "template = $template"
echo "##########################################"


date
pwd

for dtm in `ls -d ${template}`
do
   if [[ -s ${dtm}/fits_images/ ]]; then
      cd ${dtm}/fits_images/
      for ch in `ls -d ?????`
      do
         cd ${ch}
         for dirty_image in `ls dirty*_real.fits`
         do
            echo "python /software/projects/director2183/msok/imager_branches/main/imager/scripts/utils/wcs/setkey.py ${dirty_image}"
            python /software/projects/director2183/msok/imager_branches/main/imager/scripts/utils/wcs/setkey.py ${dirty_image}
         done
         cd ..    
      done
      cd ../..
   else
      echo "WARNING : ${dtm}/fits_images/ does not exist -> skipped"
   fi
done
