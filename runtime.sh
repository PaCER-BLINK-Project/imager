#!/bin/bash

if [[ -n $PAWSEY_CLUSTER && $PAWSEY_CLUSTER = "setonix" ]]; then
   echo "Setonix CPU environment"
   
   module reset
   module use /software/projects/director2183/setonix/modules/zen3/gcc/12.1.0  /software/projects/director2183/msok/setonix/modules
   module load blink_test_data/devel cfitsio/4.1.0 msfitslib/devel blink_astroio/master fftw/3.3.9 pal/0.9.8-pwmda33 libnova/0.15.0-l354muq rocm/5.4.3
   module load lfile2corrmatrix/devel
else  
   echo "Non-setonix HPC environment"

   module purge
   module load gcc/8.3.0 cascadelake  blink_test_data/devel cfitsio msfitslib pal/0.9.8 blink_astroio/master cuda/11.4.2
fi
