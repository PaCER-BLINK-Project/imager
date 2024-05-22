#!/opt/caastro/ext/anaconda/bin/python


import astropy.io.fits as pyfits
# import pylab
# import math as m
# from array import *
# import matplotlib.pyplot as plt
import numpy as np
# import string
import sys
import os
import errno
import getopt

# global parameters :
debug=0
# fitsname="re.fits"
fitsname = "1276619416_20200619163000_vis_real_channel000_time000000.fits"
fitsname2 = "1276619416_20200619163000_vis_imag_channel000_time000000.fits"

fits = pyfits.open(fitsname)
fits2 = pyfits.open(fitsname2)
x_size=fits[0].header['NAXIS1']
# channels=100
y_size=fits[0].header['NAXIS2']
print 'Read fits file %s' % fitsname
print 'FITS size = %d x %d' % (x_size,y_size)

# CRPIX2  =                    1
# CTYPE2  = 'Frequency'
# CUNIT2  = 'MHz     '
# CRVAL2  =                 200.
# CDELT2  =                  0.1
# CRPIX1  =                    1
# CTYPE1  = 'Time    '
# CUNIT1  = 'sec     '
# CDELT1  =                 0.01

data1=fits[0].data
data2=fits2[0].data
for x in range(x_size) :
   print("%.8f %.8f" % (data1[0][x],data2[0][x]))   
   
   