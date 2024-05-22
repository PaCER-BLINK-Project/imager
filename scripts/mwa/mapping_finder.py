#!/opt/caastro/ext/anaconda/bin/python


import astropy.io.fits as pyfits
# import pylab
import math
# from array import *
# import matplotlib.pyplot as plt
import numpy
# import string
import sys
import os
import errno
import getopt

def find_ant_index( value, array ) :
   l = len(array)
   for i in range(0,l):
      if math.fabs(value-array[i]) < 1e-20 :
         return i
         
   return -1
   
def find_value( value, data, x_size, y_size ) :
   # filling correlation matrix :
   for x in range(0,x_size) : # step of 2 floats to save RE/IM :
      for y in range(0,y_size ):
         val = data[y,x]
         if math.fabs(value-val) < 0.0000001 : 
            # print("FOUND %.20f at (%d,%d) = %.20f (next imaginary value is zero %s)" % (value,x,y,gpu_val,im_ok))
            return (True,x,y,val)

   return (False,-1,-1,numpy.nan)

def find_value_reim( value_re, value_im, data_re, data_im, x_size, y_size, allow_conjugate=True, limit=0.0000001 ) : # 0.0000001
   if allow_conjugate : 
      value_im = math.fabs( value_im )
   
   for x in range(0,x_size) : # step of 2 floats to save RE/IM :
      for y in range(0,y_size ):
         val_re = data_re[y,x]
         if allow_conjugate :
            val_im = math.fabs(data_im[y,x]) # allow conjugates for now 
         else :
            val_im = data_im[y,x] 
         
         if math.fabs(value_re-val_re) < limit and math.fabs(value_im-val_im) < limit : 
            # print("FOUND %.20f at (%d,%d) = %.20f (next imaginary value is zero %s)" % (value,x,y,gpu_val,im_ok))
            return (True,x,y,val_re,val_im)

   return (False,-1,-1,numpy.nan,numpy.nan)



# global parameters :
debug=0
# fitsname="re.fits"
fitsname_x_re = "1276619416_20200619163000_vis_real_channel000_time000000_pol00.fits"
fitsname_x_im = "1276619416_20200619163000_vis_imag_channel000_time000000_pol00.fits"

fitsname_y_re = "1276619416_20200619163000_vis_real_channel000_time000000_pol03.fits"
fitsname_y_im = "1276619416_20200619163000_vis_imag_channel000_time000000_pol03.fits"

fits_x_re = pyfits.open(fitsname_x_re)
fits_x_im = pyfits.open(fitsname_x_im)

fits_y_re = pyfits.open(fitsname_y_re)
fits_y_im = pyfits.open(fitsname_y_im)


x_size=fits_x_re[0].header['NAXIS1']
y_size=fits_x_re[0].header['NAXIS2']
print("Read fits files : %s %s %s %s" % (fitsname_x_re,fitsname_x_im,fitsname_y_re,fitsname_y_im))
print 'FITS size = %d x %d' % (x_size,y_size)

data_x_re=fits_x_re[0].data
data_x_im=fits_x_im[0].data
data_y_re=fits_y_re[0].data
data_y_im=fits_y_im[0].data


# data2=fits2[0].data
# for x in range(x_size) :
#   print("%.8f %.8f" % (data1[x][x],data2[x][x]))   
   

# read GPU file:
fitsname_gpu = "1276619416_20200619163000_gpubox24_00.fits"
fits_gpu = pyfits.open(fitsname_gpu)
x_size_gpu=fits_gpu[0].header['NAXIS1']
y_size_gpu=fits_gpu[0].header['NAXIS2']
print 'Read fits file %s' % fitsname_gpu
print 'FITS size = %d x %d' % (x_size_gpu,y_size_gpu)
gpu_data = fits_gpu[0].data
ch = 0
ant1=0 # 0-based index 
ant2=0 # 0-based index
n_vis = int( (ant1+1)*(ant1+1+1)/2 ) # ant1+1 is ant count at this stage
imzero_counter=0 # calculator of IM=0

n_ant=128
hdu_re = pyfits.PrimaryHDU()
hdu_re.data = numpy.zeros( (n_ant,n_ant))
hdu_im = pyfits.PrimaryHDU()
hdu_im.data = numpy.zeros( (n_ant,n_ant))

# filling correlation matrix :
for x in range (0,x_size_gpu,8) : # step of 2 floats to save RE/IM :
  xx_re = gpu_data[ch,x]
  xx_im = gpu_data[ch,x+1]
  xy_re = gpu_data[ch,x+2]
  xy_im = gpu_data[ch,x+3]
  yx_re = gpu_data[ch,x+4]
  yx_im = gpu_data[ch,x+5]
  yy_re = gpu_data[ch,x+6]
  yy_im = gpu_data[ch,x+7]
  
  if math.fabs(xx_im) < 1e-20 and math.fabs(xx_re)>1e-20 :
     imzero_counter += 1

#  if (x+1) >= n_vis :  
  if ant2 > ant1 : 
     ant1 += 1
     ant2 = 0
     n_vis = int( (ant1+1)*((ant1+1)+1)/2 ) # ant1+1 is ant count at this stage
  else :
     if ant2 == ant1 : 
        if math.fabs(xx_im) > 1e-20 :
           print("ERROR : in code auto-correlation has to have IM = 0 , ant1,ant2 = (%d,%d)" % (ant1,ant2))
           sys.exit(-1)

        # ant_idx = find_ant_index( xx_re , casa_autos )   
        #if ant_idx >= 0 :
        #   print("MAPPING : GPUBOX (visidx,ant1,ant2) = (%d,%d,%d) %.20f -> antenna %d" % (x,ant1,ant2,xx_re,ant_idx))
        #else :
        #   print("ERROR in code : auto-correlation value %.20f not found in CASA for antennas (%d,%d) !!!???" % (xx_re,ant1,ant2))

  print("%d ( (%d,%d), vis %d out of %d) ( %.8f + j%.8f )  ( %.8f + j%.8f )  ( %.8f + j%.8f )  ( %.8f + j%.8f ) " % (int(x/8),ant1,ant2,x,x_size_gpu,xx_re,xx_im,xy_re,xy_im,yx_re,yx_im,yy_re,yy_im))

  hdu_re.data[ant2,ant1] = xx_re
  hdu_re.data[ant1,ant2] = xx_re

  hdu_im.data[ant2,ant1] = xx_im
  hdu_im.data[ant1,ant2] = -xx_im

  if ant1 == ant2 :
     print("AUTOS : XX = %.8f , YY = %.8f" % (xx_re,yy_re))

  ant2 += 1

print("Number of vis with IM=0 is %d" % (imzero_counter))

# save :     
hdulist_re = pyfits.HDUList([hdu_re])
hdulist_re.writeto('gpu_re.fits',clobber=True)
hdulist_im = pyfits.HDUList([hdu_im])
hdulist_im.writeto('gpu_im.fits',clobber=True)

mapping_array=numpy.zeros( (x_size_gpu,4,2) )

# go through each value in GPU file (ch=y=0) and find it in CASA dump corr-matrix:
for x in range (0,x_size_gpu,8) : # step of 2 floats to save RE/IM :
  xx_re = gpu_data[ch,x]
  xx_im = gpu_data[ch,x+1]
  xy_re = gpu_data[ch,x+2]
  xy_im = gpu_data[ch,x+3]
  yx_re = gpu_data[ch,x+4]
  yx_im = gpu_data[ch,x+5]
  yy_re = gpu_data[ch,x+6]
  yy_im = gpu_data[ch,x+7]
  
  gpu_index = int( x / 8 )
  gpu_ant1 = ( -1 + math.sqrt(1+8*(gpu_index+1)) )/2  

  # def find_value_reim( value_re, value_im, data_re, data_im, x_size, y_size ) :
  (is_xx,ant1,ant2,val_re,val_im) = find_value_reim( xx_re, xx_im, data_x_re, data_x_im, x_size, y_size )
  
  if not is_xx :
     (is_yy,ant1,ant2,val_re,val_im) = find_value_reim( xx_re, xx_im, data_y_re, data_y_im, x_size, y_size )

  if is_xx :
     print("MAPPING : x = %d (gpu_index = %d, gpu_ant = %d) vis=(%.8f + j%.8f) found in XX matrix at (ant1,ant2) = (%d,%d)" % (x,gpu_index,gpu_ant1,xx_re,xx_im,ant1,ant2))

  if is_yy :
     print("MAPPING : x = %d (gpu_index = %d, gpu_ant = %d) vis=(%.8f + j%.8f) found in YY matrix at (ant1,ant2) = (%d,%d)" % (x,gpu_index,gpu_ant1,xx_re,xx_im,ant1,ant2))

  if not is_xx and not is_yy :
     print("ERROR : x = %d (gpu_index = %d, gpu_ant = %d) vis=(%.8f + j%.8f) not found in either XX nor YY !!!???" % (x,gpu_index,gpu_ant1,xx_re,xx_im))  
     
     # (is_xx,ant1,ant2,val_re,val_im) = find_value_reim( xx_re, -xx_im, data_x_re, data_x_im, x_size, y_size )
     # (is_yy,ant1,ant2,val_re,val_im) = find_value_reim( xx_re, -xx_im, data_y_re, data_y_im, x_size, y_size )     
     # if is_xx or is_yy :
     #   print("DEBUG : conjugate found !")
        
     # (is_xx,ant1_x,ant2_x,val_re) = find_value( xx_re, data_x_re, x_size, y_size )   
     # (is_yy,ant1_y,ant2_y,val_re) = find_value( xx_re, data_y_re, x_size, y_size )                  
     # if is_xx or is_yy :
     #   print("DEBUG : real found at (%d,%d) is_xx=%s and (%d,%d) is_yy=%s" % (ant1_x,ant2_x,is_xx,ant1_y,ant2_y,is_yy))