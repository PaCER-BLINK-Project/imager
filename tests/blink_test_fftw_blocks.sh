#!/bin/bash

# PARAMETERS :
n_images=1   # blocks for CUDA plan many 
if [[ -n "$1" && "$1" != "-" ]]; then
   n_images=$1
fi

n_pixels=180 # size of images  
if [[ -n "$2" && "$2" != "-" ]]; then
   n_pixels=$2
fi

n_iterations=100 # number of iterations (to make it run longer - maybe not the best way)
if [[ -n "$3" && "$3" != "-" ]]; then
   n_iterations=$3
fi

threads=1
if [[ -n "$4" && "$4" != "-" ]]; then
   threads=$4
fi


echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "n_images = $n_images"
echo "n_pixels = $n_pixels"
echo "n_iterations = $n_iterations"
echo "threads  = $threads"
echo "##########################################"

if [[ $n_pixels != 180 && $n_pixels != 1024 && $n_pixels != 4096 ]]; then
   echo "ERROR : there are no template images for image size $n_pixels x $n_pixels -> cannot contine, please ask MSOK to generate template images of the required size"
   exit -1
else
   echo "OK : image size ($n_pixels x $n_pixels) accepeted (template images exist)"
fi


# load test data module to set ENV variables :
module load blink_test_data/devel

# ln -s $BLINK_TEST_DATADIR/
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/u.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/v.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/w.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_real.fits vis_real.fits
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_imag.fits vis_imag.fits

# cleaning old FITS files first :
echo "rm -f re_??.fits im_??.fits"
rm -f re_??.fits im_??.fits

# TODO : change to -n -p etc :
echo "time ./fftw_blocks -n $n_images -p $n_pixels -f $n_iterations -F 1 -s ${threads}"
time ./fftw_blocks -n $n_images -p $n_pixels -f $n_iterations -F 1 -s ${threads}

exit_code=0

# for fits_real in `ls re_??.fits`
# do
#   echo "calcfits_bg ${fits_real} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_real.fits"
#   calcfits_bg ${fits_real} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_real.fits
#   
#   # repeat to get exit code :
#   equal=`calcfits_bg ${fits_real} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_real.fits | grep "Images are EQUAL" | wc -l`
#   if [[ $equal -le 0 ]]; then
#      exit_code=1
#      # to exit on the first error :
#      # break
#   fi
# done


#if [[ $exit_code -le 0 ]]; then
#   for fits_imag in `ls im_??.fits`
#   do
#      echo "calcfits_bg ${fits_imag} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_imag.fits"
#      calcfits_bg ${fits_imag} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_imag.fits
#   
#      # repeat to get exit code :
#      equal=`calcfits_bg ${fits_imag} = ${BLINK_TEST_DATADIR}/eda2/20200209/gpu/cufft_blocks/template_${n_pixels}_imag.fits | grep "Images are EQUAL" | wc -l`
#      if [[ $equal -le 0 ]]; then
#         exit_code=1
#         # to exit on the first error :
#         # break
#      fi
#   done
#fi   

# echo "Exiting script blink_test_cufft_blocks.sh with exit code = $exit_code"

# if [[ $exit_code == 0 ]]; then
#    echo "OK : test passed"
# else
#    echo "ERROR : test failed"
# fi

exit $exit_code
