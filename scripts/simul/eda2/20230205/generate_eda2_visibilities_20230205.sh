#!/bin/bash

# based on data /data/2023_02_05_ch204_corr/chan_204_20230205T034618.uv .
# date2date -ut2ux=20230205_034618
# unixtime = 1675568778

export AAVSCAL=/home/msok/aavs-calibration/

src=chan_204_20230205T034618
if [[ -n "$1" && "$1" != "-" ]]; then
   src="$1"
fi

do_cross_check=0
if [[ -n "$2" && "$2" != "-" ]]; then
   do_cross_check=$2
fi

freq_mhz=159.375
if [[ -n "$3" && "$3" != "-" ]]; then
   freq_mhz=$3
fi
# freq_mhz_int=`echo $freq_mhz | awk '{printf("%d\n",$1);}'`
freq_mhz_int=`echo ${freq_mhz} | awk '{printf "%03.0f",$1 }'`

uxtime=1675568778
if [[ -n "$4" && "$4" != "-" ]]; then
   uxtime=$4
fi


template_path=$BLINK_DIR/imager/scripts/simul/eda2/20230205/templates/
if [[ -n "$5" && "$5" != "-" ]]; then
   template_path="$5"
fi


# clean output data
echo "rm -fr *.xy *.uv *.fits *.map *.beam *.uvfits blink/"
rm -fr *.xy *.uv *.fits *.map *.beam *.uvfits blink/


echo "cp -a ${template_path}/${src}.uv ."
cp -a ${template_path}/${src}.uv .

# generate sky model using HASLAM map :
echo "$AAVSCAL/SkyModel/make_skymodel_msok.sh -f ${freq_mhz} -T ${uxtime} -B EDA2 -I 512"
$AAVSCAL/SkyModel/make_skymodel_msok.sh -f ${freq_mhz} -T ${uxtime} -B EDA2 -I 512

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_Xsky.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_Xsky.fits

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_Ysky.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_Ysky.fits


# generate visibilities using images in X and Y polarisations :
echo "uvcat vis=${src}.uv stokes=xx out=${src}_XX.uv options=nopass,nocal"
uvcat vis=${src}.uv stokes=xx out=${src}_XX.uv options=nopass,nocal

echo "uvcat vis=${src}.uv stokes=yy out=${src}_YY.uv options=nopass,nocal"
uvcat vis=${src}.uv stokes=yy out=${src}_YY.uv options=nopass,nocal

echo "uvmodel vis=${src}_XX.uv model=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv options=replace"
uvmodel vis=${src}_XX.uv model=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv options=replace

echo "uvmodel vis=${src}_YY.uv model=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv options=replace"
uvmodel vis=${src}_YY.uv model=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv options=replace

# create test images from simulated data using MIRIAD :
imsize=180
robust=-0.5

echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX.beam robust=$robust options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize}"
invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX.beam robust=$robust options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize}

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX.fits

echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY.beam robust=$robust options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize}"
invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY.beam robust=$robust options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize}

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY.fits

# naturally weighted MIRIAD images :
echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.beam options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize} sup=0"
invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.beam options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize} sup=0

echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.beam options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize} sup=0"
invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.beam options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize} sup=0

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX_Natural.fits

echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.fits"
fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY_Natural.fits

# convert .uv to .uvfits to use my code to dump correlation matrix :
echo "fits op=uvout in=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv out=${src}_simul_vis_XX.uvfits"
fits op=uvout in=${uxtime}_${freq_mhz_int}_simul_vis_XX.uv out=${src}_simul_vis_XX.uvfits

echo "fits op=uvout in=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv out=${src}_simul_vis_YY.uvfits"
fits op=uvout in=${uxtime}_${freq_mhz_int}_simul_vis_YY.uv out=${src}_simul_vis_YY.uvfits

# dump correlation matrix to be used by blink imager :
mkdir -p blink/
cd blink/
ln -s ../${src}_simul_vis_XX.uvfits ${src}_simul_XX_vis.uvfits
ln -s ../${src}_simul_vis_YY.uvfits ${src}_simul_YY_vis.uvfits

echo "python $BLINK_DIR/imager/scripts/eda2_calibration.py ${src}_simul_XX_vis.uvfits"
python $BLINK_DIR/imager/scripts/eda2_calibration.py ${src}_simul_XX_vis.uvfits

ln -s ${src}_simul_XX_vis_CorrMatrix_IM.fits ${src}_simul_XX_vis_imag.fits
ln -s ${src}_simul_XX_vis_CorrMatrix_RE.fits ${src}_simul_XX_vis_real.fits

echo "python $BLINK_DIR/imager/scripts/eda2_calibration.py ${src}_simul_YY_vis.uvfits"
python $BLINK_DIR/imager/scripts/eda2_calibration.py ${src}_simul_YY_vis.uvfits

ln -s ${src}_simul_YY_vis_CorrMatrix_IM.fits ${src}_simul_YY_vis_imag.fits
ln -s ${src}_simul_YY_vis_CorrMatrix_RE.fits ${src}_simul_YY_vis_real.fits

# create images using BLINK imager :
echo "cp ${template_path}/antenna_locations.txt ."
cp ${template_path}/antenna_locations.txt .

echo "pacer_dirty_imager ${src}_simul_XX -f ${freq_mhz} -a antenna_locations.txt -n 180 -w N -o miriad -O XX/ -Z > blink_XX.out 2>&1"
pacer_dirty_imager ${src}_simul_XX -f ${freq_mhz} -a antenna_locations.txt -n 180 -w N -o miriad -O XX/ -Z > blink_XX.out 2>&1 

echo "pacer_dirty_imager ${src}_simul_YY -f ${freq_mhz} -a antenna_locations.txt -n 180 -w N -o miriad -O YY/ -Z > blink_YY.out 2>&1"
pacer_dirty_imager ${src}_simul_YY -f ${freq_mhz} -a antenna_locations.txt -n 180 -w N -o miriad -O YY/ -Z > blink_YY.out 2>&1 
cd ../

if [[ $do_cross_check -gt 0 ]]; then
   echo "INFO : generating visibilities and simulated images using a different approach is request. WARNING : this may produce twice as much data products !"
   
   # test2 - using visibitiles generated in slightly different way:
   # test2 : generate in a slightly different way:
   echo "uvmodel vis=${src}_XX.uv model=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.uv options=replace select= \"polarization(xx)\""
   uvmodel vis=${src}_XX.uv model=${uxtime}_${freq_mhz_int}_Xsky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.uv options=replace select="polarization(xx)"

   echo "uvmodel vis=${src}_YY.uv model=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.uv options=replace select=\"polarization(yy)\""
   uvmodel vis=${src}_YY.uv model=${uxtime}_${freq_mhz_int}_Ysky.xy out=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.uv options=replace select="polarization(yy)"

   echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.beam robust=$robust options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize}"
   invert vis=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.uv map=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.map beam=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.beam robust=$robust options=double,mfs stokes=XX select='uvrange(0.0,100000)' imsize=${imsize},${imsize}

   echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.fits"
   fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.map out=${uxtime}_${freq_mhz_int}_simul_vis_XX_test2.fits

   echo "invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.beam robust=$robust options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize}"
   invert vis=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.uv map=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.map beam=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.beam robust=$robust options=double,mfs stokes=YY select='uvrange(0.0,100000)' imsize=${imsize},${imsize}

   echo "fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.fits"
   fits op=xyout in=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.map out=${uxtime}_${freq_mhz_int}_simul_vis_YY_test2.fits
else
   echo "WARNING : generating visibilities and simulated images using a different approach is not required"
fi   



