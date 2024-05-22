#!/opt/caastro/ext/anaconda/bin/python

import astropy.io.fits as pyfits
# import pylab
import math 
from array import *
# import matplotlib.pyplot as plt
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
out_fitsname="scaled.fits"
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
   print("set_keyword.py FITS_FILE FREQ" % (out_fitsname))
   print("\n")
   print("-d : increases verbose level")
   print("-h : prints help and exists")
   print("-g : produce gif of (channel-avg) for all integrations")

# functions :
def parse_command_line():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvdg", ["help", "verb", "debug", "gif"])
    except : # getopt.GetoptError, err:
        # print help information and exit:
        # print str(err) # will print something like "option -a not recognized"
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


def parse_date( fitsname ):
   dtm_utc=fitsname[12:30]
   year=int(dtm_utc[0:4])
   month=int(dtm_utc[4:6])
   day=int(dtm_utc[6:8])
   hour=int(dtm_utc[9:11])
   min=int(dtm_utc[11:13])
   sec=int(dtm_utc[13:15])
   msec=int(dtm_utc[15:18])

   # 20230601T104250099 -> 2023-08-22T17:18:41.1
   # fits[0].header['DATE-OBS'] = '2023-06-01T10:42:50.1'
   dtm_out=('%04d-%02d-%02dT%02d:%02d:%02d.%03d' % (year,month,day,hour,min,sec,msec))
   print("%s -> %s" % (dtm_utc,dtm_out))

   return dtm_out

# 
if len(sys.argv) > 1:
   fitsname = sys.argv[1]

keyword='TEST'
if len(sys.argv) > 2:
   keyword = sys.argv[2]

value='VALUE'
if len(sys.argv) > 3:
   value = sys.argv[3]

parser=optparse.OptionParser()
parser.set_usage("""setkey.py""")
parser.add_option('-i','--int','--integer',action="store_true",dest="integer",default=False, help="Integer ?")
parser.add_option('-f','--float','--float',action="store_true",dest="float",default=False, help="Float ?")
# parser.add_option("--ap_radius","-a","--aperture","--aperture_radius",dest="aperture_radius",default=0,help="Sum pixels in aperture radius [default: %default]",type="int")
# parser.add_option("--verbose","-v","--verb",dest="verbose",default=0,help="Verbosity level [default: %default]",type="int")
# parser.add_option("--outfile","-o",dest="outfile",default=None,help="Output file name [default:]",type="string")
# parser.add_option('--use_max_flux','--use_max_peak_flux','--max_flux',dest="use_max_peak_flux",action="store_true",default=False, help="Use maximum flux value around the source center [default %]")
(options,args)=parser.parse_args(sys.argv[4:])


print("####################################################")
print("PARAMTERS :")
print("####################################################")
print("fitsname       = %s"   % fitsname)
print("SET %s := %s" % (keyword,value))
print("integer = %s" % (options.integer))
print("float   = %s" % (options.float))
print("####################################################")

fits = pyfits.open(fitsname)
x_size = fits[0].header['NAXIS1']
y_size = fits[0].header['NAXIS2']

# fits[0].header[keyword] = value
# if options.integer :
#    fits[0].header[keyword] = int(value)
#else :
#    if options.float :
#        fits[0].header[keyword] = float(value)
#    else :
#        fits[0].header[keyword] = value    

# 20230601T104250099 -> 2023-08-22T17:18:41.1
# HARDCODED : fits[0].header['DATE-OBS'] = '2023-06-01T10:42:50.1'
fits[0].header['DATE-OBS'] = parse_date(fitsname)

# fits[0].header['CELLSCAL'] = 'CONSTANT'
# fits[0].header['SPECSYS'] = 'TOPOCENT'
# fits[0].header['EPOCH'] = 2000.00
# fits[0].header['PV2_1'] = -0.00197647886293958
# fits[0].header['PV2_2'] = -1.3141281461018E-06
fits[0].header['CTYPE2']  = 'DEC--SIN'

# CELLSCAL= 'CONSTANT'  /
# SPECSYS = 'TOPOCENT'  /
# LTYPE   = 'channel '  /
# LSTART  =    1.00000000000E+00  /
# LWIDTH  =    1.00000000000E+00  /
# LSTEP   =    1.00000000000E+00  /
# OBJECT  = 'eda2cal '  /
# RMS     =    9.04570400715E-01  /
# EPOCH   =    2.00000000000E+03  /


print("Writing fits with %s = %s" % (keyword,value))
fits.writeto( fitsname, overwrite=True ) # , clobber=True ) 


