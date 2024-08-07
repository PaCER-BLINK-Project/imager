#!/bin/bash

# PARAMETERS :
n_images=10   # blocks for CUDA plan many 
if [[ -n "$1" && "$1" != "-" ]]; then
   n_images=$1
fi

n_pixels=1024 # size of images  
if [[ -n "$2" && "$2" != "-" ]]; then
   n_pixels=$2
fi

n_iterations=1000 # number of iterations (to make it run longer - maybe not the best way)
if [[ -n "$3" && "$3" != "-" ]]; then
   n_iterations=$3
fi

n_streams=1
if [[ -n "$4" && "$4" != "-" ]]; then
   n_streams=$4
fi


echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "n_images = $n_images"
echo "n_pixels = $n_pixels"
echo "n_iterations = $n_iterations"
echo "n_streams = $n_streams"
echo "##########################################"



# load test data module to set ENV variables :
module load blink_test_data/devel

# ln -s $BLINK_TEST_DATADIR/
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/u.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/v.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/images/w.fits 
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_real.fits vis_real.fits
ln -s ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_imag.fits vis_imag.fits

# TODO : change to -n -p etc :
echo "time ./fftw_blocks -n $n_images -p $n_pixels -f $n_iterations -s $n_streams"
time ./fftw_blocks -n $n_images -p $n_pixels -f $n_iterations -s $n_streams

