Real EDA2 visibilities as a template (antenna positions etc) : 
   chan_204_20230205T034618.uv

   antenna_locations.txt - EDA2 antenna locations 

Sky model based on HASLAM map in orthographic projection
   1675568778_159_Xsky.fits
   1675568778_159_Ysky.fits

Dumped visibilities to FITS files :
   chan_204_20230205T034618_simul_vis_XX_CorrMatrix*fits - correlation matrix in XX polarisation
   chan_204_20230205T034618_simul_vis_YY_CorrMatrix*fits - correlation matrix in YY polarisation

MIRIAD images :
  Robust = -0.5 weighting (NOT TO BE COMPARED WITH BLINK IMAGER) :
     1675568778_159_simul_vis_XX.fits
     1675568778_159_simul_vis_YY.fits

     1675568778_159_simul_vis_XX_test2.fits - slightly differnt way of generating .uv file with visibilities
     1675568778_159_simul_vis_YY_test2.fits - slightly differnt way of generating .uv file with visibilities

  Natural weighting - these ones should be compared with BLINK imager :
     1675568778_159_simul_vis_XX_Natural.fits 
     1675568778_159_simul_vis_XX_Natural.fits


Blink images generated using command : pacer_dirty_imager ${src}_simul_XX -f 159.375 -a antenna_locations.txt -n 180 -w N -o miriad -O XX/ -Z
  1/ GPU version :
     dirty_image_20230527T062200428_XX_real_GPU.fits
     dirty_image_20230527T062200428_YY_real_GPU.fits 


  2/ CPU version :
     dirty_image_20230527T063039862_XX_real_CPU.fits
     dirty_image_20230527T063039904_YY_real_CPU.fits
