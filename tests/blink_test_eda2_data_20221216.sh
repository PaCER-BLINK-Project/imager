#!/bin/bash

# module load blink_test_data/devel

# extra_options="- -p 0.01" # temporary until GPU vs. CPU resoloved 
extra_options=""

ln -s ${BLINK_TEST_DATADIR}/eda2/20221216_101030/20221216_101030_test_100ms_ch124.LCCSPC
ln -s ${BLINK_TEST_DATADIR}/eda2/20221216_101030/20221216_101030_test_100ms_ch124.LACSPC
cp ${BLINK_TEST_DATADIR}/eda2/20221216_101030/instr_config.txt .
cp ${BLINK_TEST_DATADIR}/eda2/20221216_101030/antenna_locations.txt .
# 
pwd
which pacer_imager_multi
echo "./pacer_imager_multi 20221216_101030_test_100ms_ch124 -f 96.875 -n 180 -N 1 -v 200 -B 1 -o save_imag"
./pacer_imager_multi 20221216_101030_test_100ms_ch124 -f 96.875 -n 180 -N 1 -v 200 -B 1 -o save_imag

cd fits_images/

# different 
first_sky_re=`ls -tr 00000/dirty_image_*_real.fits | tail -1`
first_sky_im=`ls -tr 00000/dirty_image_*_imag.fits | tail -1`


# execute to have it in output :
echo "calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_real.fits ${extra_options}"
calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_real.fits ${extra_options}

# repeat to get exit code :
exit_code=0
equal=`calcfits_bg ${first_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_real.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# execute to have it in output :
echo "calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_imag.fits ${extra_options}"
calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_imag.fits ${extra_options}

# repeat to get exit code :
equal=`calcfits_bg ${first_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20221216_101030/dirty_image_imag.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "Exiting script blink_test_eda2_data_20200209.sh with exit code = $exit_code"
exit $exit_code

