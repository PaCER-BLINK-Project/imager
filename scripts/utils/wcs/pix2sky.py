#!/opt/caastro/ext/anaconda/bin/python

from __future__ import print_function
import sys,os,logging,shutil,datetime,re,subprocess,math,tempfile,string,glob,platform
import os.path
from optparse import OptionParser,OptionGroup
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
    _useastropy=True
except ImportError:
    import pywcs,pyfits
    _useastropy=False
import numpy
import math
import os

from astropy import units as u
from astropy.coordinates import SkyCoord

def pix2sky_all( fits, fitsname ):
    x_size=fits[0].header['NAXIS1']
    # channels=100
    y_size=fits[0].header['NAXIS2']

    center_x = x_size/2.00
    center_y = y_size/2.00



    print('Read fits file %s' % fitsname)
    print('FITS size = %d x %d -> center at (%d,%d)' % (x_size,y_size,center_x,center_y))

    ext=0

    h=fits[ext].header
    wcs=pywcs.WCS(h)
    naxes=h['NAXIS']

    x=numpy.arange(1,h['NAXIS1']+1)
    y=numpy.arange(1,h['NAXIS2']+1)

    ff=1
    Y,X=numpy.meshgrid(y,x)

    Xflat=X.flatten()
    Yflat=Y.flatten()
    FF=ff*numpy.ones(Xflat.shape)

    print("naxes = %d" % naxes)
    if naxes >= 4 :     
       Tostack=[Xflat,Yflat,FF] 
       for i in xrange(3,naxes):
          Tostack.append(numpy.ones(Xflat.shape))
    else :
       Tostack=[Xflat,Yflat]
    pixcrd=numpy.vstack(Tostack).transpose()
                                        
    try:
       # Convert pixel coordinates to world coordinates
       # The second argument is "origin" -- in this case we're declaring we
       # have 1-based (Fortran-like) coordinates.
       if _useastropy:
          print("_useastropy = True ?")
          sky = wcs.wcs_pix2world(pixcrd, 1)  
       else:
          print("_useastropy = False ?")
          sky = wcs.wcs_pix2sky(pixcrd, 1)

    except Exception as e:
       print('Problem converting to WCS: %s' % e)

    # extract the important pieces
    ra=sky[:,0]
    dec=sky[:,1]
           
    print("sky.shape = %d x %d -> ra.shape = %d" % (sky.shape[0],sky.shape[1],ra.shape[0]))   

    # and make them back into arrays
    RA=ra.reshape(X.shape)
    Dec=dec.reshape(Y.shape)
    
    return (RA,Dec)

def pix2sky( fits, x, y, fitsname=None, debug=False, debug_level=0 ):
    x_size=fits[0].header['NAXIS1']
    # channels=100
    y_size=fits[0].header['NAXIS2']

    center_x = x_size/2.00
    center_y = y_size/2.00

    if fitsname is not None and debug :
       print('\tRead fits file %s' % fitsname)
    if debug :
       print('\tFITS size = %d x %d -> center at (%d,%d)' % (x_size,y_size,center_x,center_y))

    ext=0
    h=fits[ext].header
    wcs=pywcs.WCS(h)
#    naxes=h['NAXIS']
#    dim = fits[0].data.ndim
    dim=wcs.naxis
#    if dim <= 0 :
#       dim = fits[0].data.shape
    naxes = dim
    if debug_level > 0 :
       print("\tFITS ndim = %d" % (dim))

    if naxes == 2 :
       (ra, dec ) = wcs.all_pix2world( x, y, 1 )
    elif naxes == 3 :
       (ra, dec, rubbish1 ) = wcs.all_pix2world( x, y, numpy.array([0]), 1 )
    else :
       (ra, dec, rubbish1, rubbish2 ) = wcs.all_pix2world( x, y, numpy.array([0]), numpy.array([0]), 1 )

    return (ra,dec)



def sky2pix( fits, ra, dec, fitsname=None ):
    x_size=fits[0].header['NAXIS1']
    # channels=100
    y_size=fits[0].header['NAXIS2']

    center_x = x_size/2.00
    center_y = y_size/2.00

    if fitsname is not None :
       print('Read fits file %s' % fitsname)
    print('FITS size = %d x %d -> center at (%d,%d)' % (x_size,y_size,center_x,center_y))

    ext=0
    h=fits[ext].header
    wcs=pywcs.WCS(h)
#    naxes=h['NAXIS']
#    dim = fits[0].data.ndim
    dim=wcs.naxis
#    if dim <= 0 :
#       dim = fits[0].data.shape
    naxes = dim
    print("\tFITS ndim = %d" % (dim))

    if naxes == 2 :
       (x, y ) = wcs.all_world2pix( ra, dec, 1 )
    elif naxes == 3 :
       (x, y, rubbish1 ) = wcs.all_world2pix( ra, dec, numpy.array([0]), 1 )
    else :
       (x, y, rubbish1, rubbish2 ) = wcs.all_world2pix( ra, dec, numpy.array([0]), numpy.array([0]), 1 )

    return (x,y)



if __name__ == '__main__':

    if len(sys.argv) > 1:
       fitsname = sys.argv[1]

    x=0
    if len(sys.argv) > 2:
       x = float(sys.argv[2])

    y=0
    if len(sys.argv) > 3:
       y = float(sys.argv[3])

   
    print("####################################################")
    print("PARAMTERS :")
    print("####################################################")
    print("fitsname       = %s"   % fitsname)
    print("(x,y)       = (%.4f,%.4f)" % (x,y))
    print("####################################################")

    fits = pyfits.open(fitsname)
    x = numpy.array( [x] )
    y = numpy.array( [y] )
    (ra,dec) = pix2sky( fits, x, y  )
    print("(x,y) = (%.2f,%.2f) -> (RA,DEC) = (%.4f,%.4f) [deg]" % (x[0],y[0],ra[0],dec[0]))


#    print "RA.shape = %d x %d" % (RA.shape[0],RA.shape[1])
#    print "Dec.shape = %d x %d" % (Dec.shape[0],Dec.shape[1])

#    ra_xc = RA[xc,yc]
#    if ra_xc < 0 :
#      ra_xc = ra_xc + 360.0
#    dec_xc = Dec[xc,yc]

#    c = SkyCoord(ra=ra_xc*u.degree, dec=dec_xc*u.degree, frame='icrs')
#    dec_dms=c.to_string('dms')
#    radec_str=c.to_string('hmsdms')

#    print "RADEC of (%d,%d) = (%.4f,%.4f) = %s" % (xc,yc,ra_xc,dec_xc,radec_str)

#    hdu_out = pyfits.PrimaryHDU()
#    hdu_out.header = fits[0].header
#    RA[RA<0] += 360
#    hdu_out.data=RA.transpose().copy()
#    hdulist = pyfits.HDUList([hdu_out])
#    hdulist.writeto("ra.fits",clobber=True)

#    hdu_out = pyfits.PrimaryHDU()
#    hdu_out.header = fits[0].header
#    hdu_out.data=Dec.transpose().copy()
#    hdulist = pyfits.HDUList([hdu_out])
#    hdulist.writeto("dec.fits",clobber=True)





        

       
