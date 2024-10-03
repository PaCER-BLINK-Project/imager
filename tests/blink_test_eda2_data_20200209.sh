#!/bin/bash

extra_options="- -p 0.01 -d -d -d" # temporary until GPU vs. CPU resoloved 
if [[ -n "$PAWSEY_CLUSTER" ]]; then
   echo "module load blink_test_data/devel"
   module load blink_test_data/devel
   
   if [[ -n "$PAWSEY_CLUSTER" && "$PAWSEY_CLUSTER" == "setonix" ]]; then
      extra_options="- -p 0.2 -d -d -d" # the differences are slightly higher on setonix 
   fi
else
   echo "INFO : variable PAWSEY_CLUSTER is not defined not using -> modules (desktop/server like environment)"
fi

# overwrite BLINK_TEST_DATADIR with ../data/ to use local test data from ../data/eda2/20200209 :
BLINK_TEST_DATADIR=../data/

cp ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_????.fits .
cp ${BLINK_TEST_DATADIR}/eda2/20200209/antenna_locations.txt .
cp ${BLINK_TEST_DATADIR}/eda2/20200209/calsol_merged.txt .

# 
pwd
which pacer_dirty_imager
echo "./pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 180 -w N -o miriad -c calsol_merged.txt -Z -v 100 -O 20200209_eda2/"
./pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 180 -w N -o miriad -c calsol_merged.txt -Z -v 100 -O 20200209_eda2/

last_sky_re=`ls -tr 20200209_eda2/dirty_image_*_real.fits | tail -1`
last_sky_im=`ls -tr 20200209_eda2/dirty_image_*_imag.fits | tail -1`

# execute to have it in output :
echo "calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_real.fits ${extra_options}"
calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_real.fits ${extra_options}

# repeat to get exit code :
exit_code=0
equal=`calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_real.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# execute to have it in output :
echo "calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_imag.fits ${extra_options}"
calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_imag.fits ${extra_options}

# repeat to get exit code :
equal=`calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/eda2/20200209/dirty_image_imag.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "Exiting script blink_test_eda2_data_20200209.sh with exit code = $exit_code"
exit $exit_code
