#!/bin/bash

extra_options="" # temporary until GPU vs. CPU resoloved 
if [[ -n "$PAWSEY_CLUSTER" ]]; then
   # echo "module load blink_test_data/devel"
   # module load blink_test_data/devel
   
   if [[ -n "$PAWSEY_CLUSTER" && "$PAWSEY_CLUSTER" == "setonix" ]]; then
      extra_options="- -p 4 -d -d -d" # the differences are slightly higher on setonix 
   fi
else
   echo "INFO : variable PAWSEY_CLUSTER is not defined not using -> modules (desktop/server like environment)"
fi


ln -s $BLINK_TEST_DATADIR/mwa/simul/hydra/HydA_154MHz.mwa_ant_vis_real_channel000_time000000_pol0.fits
ln -s $BLINK_TEST_DATADIR/mwa/simul/hydra/HydA_154MHz.mwa_ant_vis_imag_channel000_time000000_pol0.fits

# 
pwd
which pacer_dirty_imager
echo "./pacer_dirty_imager HydA_154MHz.mwa_ant -n 378 -f 154.88 -F 13.86 -p  _channel000_time000000_pol0 -M $BLINK_TEST_DATADIR/mwa/simul/hydra/metadata.txt -V 100 -v 100 -w N -a $BLINK_TEST_DATADIR/mwa/simul/hydra/mwa_antenna_locations.txt -E -J dirty_image_hydraa_simul -O hydra_simul_test/"
./pacer_dirty_imager HydA_154MHz.mwa_ant -n 378 -f 154.88 -F 13.86 -p  _channel000_time000000_pol0 -M $BLINK_TEST_DATADIR/mwa/simul/hydra/metadata.txt -V 100 -v 100 -w N -a $BLINK_TEST_DATADIR/mwa/simul/hydra/mwa_antenna_locations.txt -E -J dirty_image_hydraa_simul -O hydra_simul_test/

last_sky_re=`ls -tr hydra_simul_test/dirty_image_hydraa_simul_real.fits | tail -1`
last_sky_im=`ls -tr hydra_simul_test/dirty_image_hydraa_simul_imag.fits | tail -1`

# execute to have it in output :
echo "calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra//images/dirty_image_simul_hydra_real.fits ${extra_options}"
calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra/images/dirty_image_simul_hydra_real.fits ${extra_options}

# repeat to get exit code :
exit_code=0
equal=`calcfits_bg ${last_sky_re} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra/images/dirty_image_simul_hydra_real.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

# execute to have it in output :
echo "calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra/images/dirty_image_simul_hydra_imag.fits ${extra_options}"
calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra/images/dirty_image_simul_hydra_imag.fits ${extra_options}

# repeat to get exit code :
equal=`calcfits_bg ${last_sky_im} = ${BLINK_TEST_DATADIR}/mwa/simul/hydra/images/dirty_image_simul_hydra_imag.fits ${extra_options} | grep "Images are EQUAL" | wc -l`
if [[ $equal -le 0 ]]; then
   exit_code=1
fi

echo "Exiting script blink_test_mwa_simul_hydra.sh with exit code = $exit_code"
exit $exit_code
