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

casa_autos=[119448.0000,133580.7500,116702.7500,160176.0000,142252.5000,167756.2500,145155.2500,136278.2500,151409.5000,167463.0000,169468.2500,242343.5000,141075.5000,165487.2500,122375.7500,196583.5000,183161.2500,171682.5000,196191.2500,174761.0000,203969.5000,195011.5000,175193.7500,125399.7500,179297.7500,172837.5000,76575.5000,137492.7500,125328.2500,190963.5000,142012.5000,148432.0000,163298.2500,187270.5000,179644.2500,183631.5000,154322.7500,190315.0000,156500.0000,154259.0000,166652.2500,152703.2500,163290.0000,164028.5000,156783.5000,188054.0000,146688.5000,129199.0000,173291.0000,159845.7500,20271.7500,131479.7500,172248.2500,194296.5000,166811.7500,167659.2500,136442.0000,112583.5000,115301.2500,170547.7500,101825.2500,152018.2500,144496.7500,146991.7500,210285.0000,183970.5000,167831.0000,147015.7500,157989.0000,183903.7500,180366.0000,19442.7500,125185.2500,94135.7500,77534.7500,105600.0000,186188.2500,57734.2500,172250.2500,145955.2500,34846.0000,56502.7500,70244.0000,127714.0000,139072.2500,123637.0000,68420.0000,151034.0000,96232.5000,57615.5000,189113.5000,134545.7500,150582.7500,102077.0000,193315.2500,121740.0000,155028.0000,46626.2500,151264.5000,124094.7500,104574.0000,33048.7500,124783.0000,142705.7500,115125.0000,17780.5000,74063.5000,130080.7500,0.0000,84115.2500,126535.5000,50550.2500,74370.5000,98443.2500,0.0000,144062.5000,70233.5000,72266.5000,139187.5000,12139.2500,189862.2500,128239.0000,72394.7500,32906.0000,33746.5000,91778.2500,208084.0000,121911.0000]
casa_autos=numpy.array(casa_autos)

def find_ant_index( value, array ) :
   l = len(array)
   for i in range(0,l):
      if math.fabs(value-array[i]) < 1e-20 :
         return i
         
   return -1
   
def find_value( value, data ) :
   # filling correlation matrix :
   for x in range(0,x_size) : # step of 2 floats to save RE/IM :
      for y in range(0,y_size ):
         gpu_val = gpu_data[y,x]
         if math.fabs(value-gpu_val) < 0.00001 : 
            gpu_val_im = gpu_data[y,x+1]
            im_ok=True
            if math.fabs(gpu_val_im) > 1e-20 :
               im_ok=False

            # print("FOUND %.20f at (%d,%d) = %.20f (next imaginary value is zero %s)" % (value,x,y,gpu_val,im_ok))
            return (True,x,y,gpu_val,gpu_val_im)

   return False

value=0.00
if len(sys.argv) > 1:
   value = float( sys.argv[1] )

# global parameters :
debug=0
# fitsname="re.fits"
fitsname = "1276619416_20200619163000_gpubox24_00.fits"

fits = pyfits.open(fitsname)
x_size=fits[0].header['NAXIS1']
# channels=100
y_size=fits[0].header['NAXIS2']
print 'Read fits file %s' % fitsname
print 'FITS size = %d x %d' % (x_size,y_size)

gpu_data = fits[0].data


# find_value( value, gpu_data )
ant_idx=0
#for value in casa_autos :
for ant_index in range(0,len(casa_autos)):
   value = casa_autos[ant_index]
   (found,x,y,gpu_val,gpu_val_im) = find_value( value, gpu_data )
   gpu_index = int( x / 8 )
# same equation in : /home/msok/mwa_software/MWA_Tools/build_lfiles/build_lfiles.c
#         /* update globals based on in-file metadata */
#        /* The number of correlation products is slightly different to the standard n(n+1)/2
#           formula since there is a small fraction of redundant info in the file. If n is the
#           number of stations (not stations*pols), and N is the dimension of the data in the file
#           then n = (-1 + sqrt(1+N))/2
#        */
# here 8 due to 8 floats for XX XY YX YY in re/im (4*2 = 8 floats) : 
   gpu_ant1 = ( -1 + math.sqrt(1+8*(gpu_index+1)) )/2
   im_zero=True
   
   if found :
      if math.fabs(gpu_val_im) > 1e-20 :
         im_zero = False
      print("FOUND %.20f at (%d,%d) = %.20f (next imaginary value is zero %.20f -> %s) : casa_ant_index = %d, gpu_index = %d -> gpu_ant1=%d" % (value,x,y,gpu_val,gpu_val_im,im_zero,ant_index,gpu_index,gpu_ant1))
   