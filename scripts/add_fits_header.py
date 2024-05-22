

import astropy.io.fits as pyfits
import pylab
import math 
from array import *
import matplotlib.pyplot as plt
import numpy as np
import string
import sys
import os
import errno
import getopt
import optparse


# global parameters :
debug=0
fitsname="file.fits"
out_fitsname="fits_keys.fits"
do_show_plots=0
do_gif=0

center_x=1025
center_y=1025
radius=600

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise
                                            
def usage():
   print "add_fits_header.py FITS_FILE" % out_fitsname
   print "\n"
   print "-d : increases verbose level"
   print "-h : prints help and exists"
   print "-g : produce gif of (channel-avg) for all integrations"

# functions :
def parse_command_line():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvdg", ["help", "verb", "debug", "gif"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    for o, a in opts:
        if o in ("-d","--debug"):
            debug += 1
        if o in ("-v","--verb"):
            debug += 1
        if o in ("-g","--gif"):
            do_gif = 1
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    # ...

# 
if len(sys.argv) > 1:
   fitsname = sys.argv[1]

pixscale=0.52941176
if len(sys.argv) > 2:
   pixscale = float( sys.argv[2] )


parser=optparse.OptionParser()
parser.set_usage("""add_fits_header.py""")
# parser.add_option('-i','--int','--integer',action="store_true",dest="integer",default=False, help="Integer ?")
# parser.add_option('-f','--float','--float',action="store_true",dest="float",default=False, help="Float ?")
# parser.add_option("--ap_radius","-a","--aperture","--aperture_radius",dest="aperture_radius",default=0,help="Sum pixels in aperture radius [default: %default]",type="int")
# parser.add_option("--verbose","-v","--verb",dest="verbose",default=0,help="Verbosity level [default: %default]",type="int")
# parser.add_option("--outfile","-o",dest="outfile",default=None,help="Output file name [default:]",type="string")
# parser.add_option('--use_max_flux','--use_max_peak_flux','--max_flux',dest="use_max_peak_flux",action="store_true",default=False, help="Use maximum flux value around the source center [default %]")
(options,args)=parser.parse_args(sys.argv[4:])


print "####################################################"
print "PARAMTERS :"
print "####################################################"
print "fitsname       = %s"   % fitsname
print "pixscale       = %.6f [deg]" % (pixscale) 
print "####################################################"

fits = pyfits.open(fitsname)
x_size = fits[0].header['NAXIS1']
y_size = fits[0].header['NAXIS2']

# TELESCOP= 'eda2    '  /
# CRPIX1  =    9.10000000000E+01  /
# CDELT1  =   -7.71037590484E-01  /
# CRVAL1  =    2.13673600000E+02  /
# CTYPE1  = 'RA---SIN'  /
# CRPIX2  =    9.10000000000E+01  /
# CDELT2  =    7.71037590484E-01  /
# CRVAL2  =   -2.67033000000E+01  /
# CTYPE2  = 'DEC--SIN'  /

sign_delta = 1

fits[0].header['TELESCOP'] = 'eda2'
fits[0].header['CRPIX1'] = int(x_size/2) + 1
fits[0].header['CDELT1'] = sign_delta*pixscale # -7.71037590484E-01
fits[0].header['CRVAL1'] = 2.13673600000E+02
fits[0].header['CTYPE1'] = 'RA---SIN'

fits[0].header['CRPIX2'] = int(y_size/2) + 1
fits[0].header['CDELT2'] = sign_delta*pixscale # -7.71037590484E-01
fits[0].header['CRVAL2'] = -2.67033000000E+01
fits[0].header['CTYPE2'] = 'DEC---SIN'

fits.writeto( out_fitsname, clobber=True ) 
