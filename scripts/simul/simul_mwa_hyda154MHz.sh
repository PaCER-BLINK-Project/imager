#!/bin/bash

# WARNING : this only works correctly with OLD CASA 
# export PATH=/home/msok/.casa/bin:/opt/caastro/ext/dload/casa/current/bin/:$PATH
# export CASAPATH=/opt/caastro/ext/dload/casa/current/
# export PAGER=less
# alias casapy='casa'
# NEW CASA CRASHES -> GREAT !!!!


# simulate image using source model in image HydA_154.88MHz.im
# # for 154 MHz , 3000 km baseline
casapy -c "simobserve( project='HydA_154MHz',skymodel='HydA_154.88MHz.im',totaltime='1s',integration='1s', antennalist='mwa_ant.cfg', thermalnoise='tsys-manual', t_sky=200, tau0=0.00, obsmode='int', incell='2.2arcmin', mapsize=['3deg','3deg'] )"

# create images from the simulated visibilities :
casapy -c "simanalyze('HydA_154MHz',niter=0,weighting = 'natural',showuv=True,showdifference=True)"

# dump visibilities to FITS files suitable for the PACER imager :
cd HydA_154MHz/
echo "casapy -c ../../get_corrmatrix_from_casa.py --save_full_matrix HydA_154MHz.mwa_ant.ms"
casapy -c ../../get_corrmatrix_from_casa.py --save_full_matrix HydA_154MHz.mwa_ant.ms

echo "wsclean -name wsclean378x378 -j 1 -size 378 378 -pol XX -absmem 16 -weight natural -scale 2.2arcmin -niter 0  -nwlayers 1 -save-uv  -grid-mode nn -save-weights HydA_154MHz.mwa_ant.ms"
wsclean -name wsclean378x378 -j 1 -size 378 378 -pol XX -absmem 16 -weight natural -scale 2.2arcmin -niter 0  -nwlayers 1 -save-uv  -grid-mode nn -save-weights HydA_154MHz.mwa_ant.ms

echo "ds9 -scale zscale wsclean378x378-dirty.fits &"
ds9 -scale zscale wsclean378x378-dirty.fits &

# BLINK IMAGER AS INSTALLED ON THE SYSTEM OR MODULE LOADED :
echo "pacer_dirty_imager HydA_154MHz.mwa_ant -n 378 -f 154.88 -F 13.86 -p  _channel000_time000000_pol0 -M metadata.txt -V 100 -w N -a mwa_antenna_locations.txt -E -U 1400665108"
pacer_dirty_imager HydA_154MHz.mwa_ant -n 378 -f 154.88 -F 13.86 -p  _channel000_time000000_pol0 -M metadata.txt -V 100 -w N -a mwa_antenna_locations.txt -E -U 1400665108

last_fits=`ls dirty*_real.fits | head -1`
echo "ds9 -scale zscale $last_fits &"
ds9 -scale zscale $last_fits &
