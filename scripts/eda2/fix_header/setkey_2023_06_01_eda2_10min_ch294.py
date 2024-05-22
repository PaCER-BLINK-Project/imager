#!/opt/caastro/ext/anaconda/bin/python

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

from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, FK5, FK4
MWA_POS=EarthLocation.from_geodetic(lon="116:40:14.93",lat="-26:42:11.95",height=377.8)

DEG2RAD = (math.pi/180.00)
RAD2DEG = (180.00/math.pi)


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
   print("set_keyword.py FITS_FILE FREQ" % out_fitsname)
   print("\n")
   print("-d : increases verbose level")
   print("-h : prints help and exists")
   print("-g : produce gif of (channel-avg) for all integrations")

# https://python.hotexamples.com/examples/astropy.coordinates/AltAz/-/python-altaz-class-examples.html
# def azel2radec(az_deg, el_deg, lat_deg, lon_deg, dtime):
#    obs = EarthLocation(lat=lat_deg * u.deg, lon=lon_deg * u.deg)
#    direc = AltAz(location=obs, obstime=Time(dtime),
#                  az=az_deg * u.deg, alt=el_deg * u.deg)
#    sky = SkyCoord(direc.transform_to(ICRS()))
#    return sky.ra.deg, sky.dec.deg
# t0 = Time("%d-%02d-%02d 00:00:00" % (year, month, day), format = 'iso', scale = 'utc')
def azim2radec( az_deg, el_deg, dtime, debug=True, astro_azim=True ) : 
   # obstime=Time(dtime)
   fk5_2000 = FK5(equinox=Time(2000, format='jyear'))
   fk4_2000 = FK4(equinox=Time(2000, format='jyear'))
   direc = AltAz(location=MWA_POS, obstime=dtime, az=az_deg * u.deg, alt=el_deg * u.deg)
   sky = SkyCoord(direc.transform_to(ICRS()))
#   sky = SkyCoord(direc.transform_to(fk4_2000))

   print("DEBUG : ra, dec = %.6f, %.6f [deg]" % (sky.ra.deg, sky.dec.deg))
  
   return ( sky.ra.deg, sky.dec.deg, az_deg, el_deg, dtime  )   



# functions :
def parse_command_line():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hvdg", ["help", "verb", "debug", "gif"])
#    except getopt.GetoptError, err:
    except:
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

# 
if len(sys.argv) > 1:
   fitsname = sys.argv[1]

parser=optparse.OptionParser()
parser.set_usage("""setkey.py""")
parser.add_option('-i','--int','--integer',action="store_true",dest="integer",default=False, help="Integer ?")
parser.add_option('-f','--float','--float',action="store_true",dest="float",default=False, help="Float ?")
parser.add_option('-r','--ra','--ra_deg',dest="ra_deg",default=165.41518729, help="RA [deg]",type="float")
parser.add_option('-d','--dec','--dec_deg',dest="dec_deg",default=-26.703319, help="DEC [deg]",type="float")

parser.add_option('-a','--azim','--az_deg','--azim_deg',dest="az_deg",default=0.00, help="Azim [deg]",type="float")
parser.add_option('-e','--elev','--el_deg',dest="el_deg",default=90.00, help="Elev [deg]",type="float")

parser.add_option('--use_filename','--dtmfile',action="store_true",dest="use_filename",default=False, help="Get UTC from file name [default %default]")
(options,args)=parser.parse_args(sys.argv[2:])


print("####################################################")
print("PARAMTERS :")
print("####################################################")
print("fitsname       = %s"   % fitsname)
print("integer = %s" % (options.integer))
print("float   = %s" % (options.float))
print("(AZ,EL)  = (%.6f,%.6f) [deg]" % (options.az_deg,options.el_deg))
print("(RA,DEC) = (%.6f,%.6f) [deg]" % (options.ra_deg,options.dec_deg))
print("use_filename = %d" % (options.use_filename))
print("####################################################")

fits = pyfits.open(fitsname)
x_size = fits[0].header['NAXIS1']
y_size = fits[0].header['NAXIS2']

# dirty_image_20230601T103621900_real.fits
if options.use_filename :
   if fitsname[0:11] == "dirty_image" :
      dtime=fitsname[12:30]
      if dtime[8] == "T" :
         # 20230601T103621900
         year=dtime[0:4]    
         mon=dtime[4:6]
         day=dtime[6:8]
         h=dtime[9:11]
         m=dtime[11:13]
         s=dtime[13:15] + "." + dtime[15:]
         dtime_new=("%s-%s-%s %s:%s:%s" % (year,mon,day,h,m,s))
         print("DEBUG : %s -> %s" % (dtime,dtime_new))
         dtime = dtime_new
         
         print("DEBUG : before call of azim2radec ?")   
         ( options.ra_deg, options.dec_deg, azim_deg, elev_deg, dtm_utc ) = azim2radec( options.az_deg, options.el_deg, dtime  )   
         print("DEBUG : new radec = (%.6f,%.6f) [deg]" % (options.ra_deg,options.dec_deg))



fits[0].header["CDELT1"] = 1.0 # float(180.00/120.00)
fits[0].header["CDELT2"] = 1.0 # float(180.00/120.00)
fits[0].header["CRVAL1"] = options.ra_deg
fits[0].header["CRVAL2"] = options.dec_deg


fits.writeto( fitsname, clobber=True ) 


