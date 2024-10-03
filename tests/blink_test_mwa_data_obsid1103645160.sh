#!/bin/bash

extra_options="- -p 0.01" # temporary until GPU vs. CPU resoloved 
if [[ -n "$PAWSEY_CLUSTER" ]]; then
   # echo "module load blink_test_data/devel"
   # module load blink_test_data/devel
   
   if [[ -n "$PAWSEY_CLUSTER" && "$PAWSEY_CLUSTER" == "setonix" ]]; then
      extra_options="- -p 5 -d -d -d" # the differences are slightly higher on setonix 
   fi
else
   echo "INFO : variable PAWSEY_CLUSTER is not defined not using -> modules (desktop/server like environment)"
fi

# overwrite BLINK_TEST_DATADIR with ../data/ to use local test data from ../data/mwa/1103645160/ :
BLINK_TEST_DATADIR=../data/

ln -s $BLINK_TEST_DATADIR/mwa/1103645160/1103645160_vis_real_channel000_time000000.fits
ln -s $BLINK_TEST_DATADIR/mwa/1103645160/1103645160_vis_imag_channel000_time000000.fits

# 
pwd
which pacer_dirty_imager
echo "./pacer_dirty_imager 1103645160 -p _channel000_time000000 -n 456 -f 138.88 -F 16.72 -w N -a $BLINK_TEST_DATADIR/mwa/1103645160/antenna_locations.txt -M $BLINK_TEST_DATADIR/mwa/1103645160/1103645160.metafits -v 100 -U 1419609943 -O mwa_obsid1103645160"
./pacer_dirty_imager 1103645160 -p _channel000_time000000 -n 456 -f 138.88 -F 16.72 -w N -a $BLINK_TEST_DATADIR/mwa/1103645160/antenna_locations.txt -M $BLINK_TEST_DATADIR/mwa/1103645160/1103645160.metafits -v 100 -U 1419609943 -O mwa_obsid1103645160

last_sky_re=`ls -tr mwa_obsid1103645160/dirty_image_*_real.fits | tail -1`
last_sky_im=`ls -tr mwa_obsid1103645160/dirty_image_*_imag.fits | tail -1`

# execute to have it in output :
echo "calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_real.fits ${extra_options}"
calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_real.fits ${extra_options}

# repeat to get exit code :
exit_code=0
equal=`calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_real.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# execute to have it in output :
echo "calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_imag.fits ${extra_options}"
calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_imag.fits ${extra_options}

# repeat to get exit code :
equal=`calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/1103645160/images/dirty_image_imag.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "Exiting script blink_test_mwa_data_obsid1103645160.sh with exit code = $exit_code"
exit $exit_code
