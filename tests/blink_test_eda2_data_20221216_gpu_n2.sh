#!/bin/bash

module load blink_test_data/devel

# extra_options="- -p 0.01" # temporary until GPU vs. CPU resoloved 
extra_options=""

ln -s ${BLINK_TEST_DATADIR}/eda2/20221216_101030/20221216_101030_test_100ms_ch124.LCCSPC
ln -s ${BLINK_TEST_DATADIR}/eda2/20221216_101030/20221216_101030_test_100ms_ch124.LACSPC
cp ${BLINK_TEST_DATADIR}/eda2/20221216_101030/instr_config.txt .
cp ${BLINK_TEST_DATADIR}/eda2/20221216_101030/antenna_locations.txt .
# 
pwd
which pacer_imager_multi
echo "./pacer_imager_multi 20221216_101030_test_100ms_ch124 -f 96.875 -n 180 -N 2 -G -Z -v 200 -B 1 -o save_imag"
./pacer_imager_multi 20221216_101030_test_100ms_ch124 -f 96.875 -n 180 -N 2 -G -Z -v 200 -B 1 -o save_imag

# first_sky_re=gpumulti_dirty_image_ch000_real.fits
# first_sky_im=gpumulti_dirty_image_ch000_imag.fits
first_sky_re="ch00000/gpumulti_dirty_image_ch000_20221216T101030100_real.fits"
first_sky_im="ch00000/gpumulti_dirty_image_ch000_20221216T101030100_imag.fits"

if [[ ! -s ${first_sky_re} ]]; then
   echo "ERROR : real sky image (file $first_sky_re) not created"
   exit 1
fi

if [[ ! -s ${first_sky_im} ]]; then
   echo "ERROR : imaginary sky image (file $first_sky_im) not created"
   exit 1
fi


if [[ ! -s ${first_sky_re} ]]; then
   echo "ERROR : real sky image (file $first_sky_re) not created"
   exit 1
fi

if [[ ! -s ${first_sky_im} ]]; then
   echo "ERROR : imaginary sky image (file $first_sky_im) not created"
   exit 1
fi

# execute to have it in output :
echo "calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_real.fits ${extra_options}"
calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_real.fits ${extra_options}

# repeat to get exit code :
exit_code=0
equal=`calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_real.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# execute to have it in output :
echo "calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_imag.fits ${extra_options}"
calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_imag.fits ${extra_options}

# repeat to get exit code :
equal=`calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_20221216T101030100_imag.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# Cannot compare UV real for n=2 because it's only saved for n=1 in the CPU version ???
# echo "calcfits_bg gpu_uvgrid_real_ch000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_real_180x180.fits - - -p 0.001"
# calcfits_bg gpu_uvgrid_real_ch000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_real_180x180.fits - - -p 0.001

# echo "calcfits_bg gpu_uvgrid_imag_ch000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_imag_180x180.fits"
# calcfits_bg gpu_uvgrid_imag_ch000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_imag_180x180.fits

echo "calcfits_bg counter_channel000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_counter_180x180.fits"
calcfits_bg counter_channel000.fits = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/uv_grid_counter_180x180.fits


echo "Exiting script blink_test_eda2_data_20200209.sh with exit code = $exit_code"
exit $exit_code

