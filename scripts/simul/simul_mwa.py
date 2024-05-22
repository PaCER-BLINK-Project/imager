import subprocess, re
import numpy as n,os,sys,shutil
from pylab import *
import optparse

# create command line arguments
parser=optparse.OptionParser()
parser.set_usage("Usage: casapy [--nologger] -c simul_mwa.py [options]")
# parse command line arguments
# parser.add_option("--exclude_AAVS",dest="exclude_AAVS",default=False,action="store_true",
#                          help="If True, exclude AAVS (Tile079) from calibration. [default: %default]")
# parser.add_option("--uvdist_low",dest="uvdist_low",default=30,
#                          help="Minimum baseline UV distance (wavelengths) to be included. [default: %default]")
# parser.add_option("--uvdist_high",dest="uvdist_high",default=1e6,
#                          help="Maximum baseline UV distance (wavelengths) to be included. [default: %default]")
# casa_index=sys.argv.index('-c')+2  # remove CASA options
# (options,args)=parser.parse_args(sys.argv[casa_index:])

# Overwrite files
clobber=True

# Get the frequency information of the measurement set
# ms.open(vis)
# rec = ms.getdata(['axis_info'])
# df,f0 = (rec['axis_info']['freq_axis']['resolution'][len(rec['axis_info']['freq_axis']['resolution'])/2],rec['axis_info']['freq_axis']['chan_freq'][len(rec['axis_info']['freq_axis']['resolution'])/2])
# F =rec['axis_info']['freq_axis']['chan_freq'].squeeze()/1e6
# df=df[0]*len(rec['axis_info']['freq_axis']['resolution'])
# f0=f0[0]
# rec_time=ms.getdata(['time'])
# sectime=qa.quantity(rec_time['time'][0],unitname='s')
sectime=1661505463
time_utc="2022-08-26 09:17:43"

f0=121*1.28*1000000.00 # Hz 
df=24*1.28*1000000.00 # Hz 
midfreq=str(f0/1.e6)+'MHz'
bandwidth=str(df)+'Hz'
freq_array=[midfreq]
spectral_beam=False

calibrator="HydA"
modeldir="./"

# Start models are imagefreq Jy/pixel fits files in a known directory
model=modeldir+calibrator+'.fits'
# With a corresponding spectral index map
spec_index=modeldir+calibrator+'_spec_index.fits'
if not os.path.exists(spec_index):
    print 'Could not find spectral index map: %s'%spec_index
    print 'Using -0.83 for Hydra A'
    #Create uniform spectral map same size as HydrA image
    immath(imagename=[model],mode='evalexpr',outfile=spec_index[:-5]+'.im',
           expr='(IM0*0-0.83)')
    exportfits(fitsimage=spec_index,imagename=spec_index[:-5]+'.im',overwrite=clobber)
    print 'Exported to %s'%spec_index
    spec_index=modeldir+calibrator+'_spec_index.fits'

print 'Using spectral index map: %s'%spec_index 

if not os.path.exists(model):
  print 'Could not find calibrator model %s' % model
  raise KeyboardInterrupt

#Get centre frequency of model 
#importfits(fitsimage=model,imagename='test.im')

try:
    # See if there is a corresponding im file
    model_im=modeldir+calibrator+'.im'
    imagefreq=imhead(model_im,mode='get',hdkey='crval3')['value'] #Freq in Hz
except:
    print "Could not find image frequency"
    raise KeyboardInterrupt
print 'Image frequency is %s'%imagefreq
#  imhead(imagename=outname,mode='put',hdkey='crval3',hdvalue='150MHz')

#===============================================================================
# # Generate the primary beam
# delays=info.delays
# str_delays=','.join(map(str,delays))
# print 'Delays are: %s' % str_delays
#===============================================================================

# do this for the start, middle, and end frequencies
for freq in freq_array:
# We'll generate images in the local directory at the right frequency for this ms
  outname=calibrator+'_'+freq+'.im' #This is the name for the model image
  outnt2=calibrator+'_'+freq+'_nt2.im'
  print 'Generating image %s for calibration %s' % (outname,"")
  print 'Generating 2nd-order term image %s for calibration %s' % (outnt2,"")
# import model, edit header so make_beam generates the right beam in the right place
  if os.path.exists(outname) and clobber:
    print 'Overwriting %s' % (outname)
    rmtables(outname)
    #Why don't we remove outnt2?

  importfits(fitsimage=model,imagename=outname)
  imhead(outname,mode='put',hdkey='crval3',hdvalue=freq)
  imhead(outname,mode='put',hdkey='cdelt3',hdvalue=bandwidth)
  #TMC: the qa.time output is broken and causes an error 
  #imhead(outname,mode='put',hdkey='date-obs',hdvalue=qa.time(sectime,form=["ymd"])) 
   
  #=============================================================================
  # exportfits(fitsimage=outname+'.fits',imagename=outname,overwrite=clobber)
  # print 'Creating primary beam models...'
  # subprocess.call(['python',bindir+'make_beam.py','-f',outname+'.fits','-d',str_delays])
  #=============================================================================
  
# delete the temporary model
#  rmtables(outname)

#===============================================================================
#  beamimage={}
#  fitsimage={}
#  for stokes in ['XX','YY']:
# # import the beams from make_beam.py into beamimage dictionary (for XX and YY)
#        fitsimage[stokes]=calibrator+'_'+freq+'.im_beam'+stokes+'.fits'
#        beamimage[stokes]=calibrator+'_'+freq+'_beam'+stokes+'.im'
#        if os.path.exists(beamimage[stokes]) and clobber:
#          rmtables(beamimage[stokes])
#        importfits(fitsimage=fitsimage[stokes],imagename=beamimage[stokes],
#                   overwrite=clobber)
#===============================================================================

# scale by the primary beam
# Correct way of doing this is to generate separate models for XX and YY
# Unfortunately, ft doesn't really understand cubes
# So instead we just use the XX model, and then scale the YY solution later

#freq=midfreq
#outname=calibrator+'_'+freq+'.im'
#outnt2=calibrator+'_'+freq+'_nt2.im'

#===============================================================================
# beamarray=[calibrator+'_'+freq+'_beamXX.im',calibrator+'_'+freq+'_beamYY.im']
# 
# # Hardcoded to use the XX beam in the model
# beam=calibrator+'_'+freq+'_beamXX.im'
# ratio=calibrator+'_'+freq+'_beam_ratio.im'
# # divide to make a ratio beam, so we know how to scale the YY solution later
# if os.path.exists(ratio) and clobber:
#	rmtables(ratio)
# immath(imagename=beamarray,mode='evalexpr',expr='(IM0/IM1)',outfile=ratio)
# ratio=imstat(ratio)['mean'][0]
#===============================================================================

if os.path.exists(outname) and clobber:
    rmtables(outname)

# Models are at 150MHz
# Generate scaled image at correct frequency
#exp='IM0/((150000000/'+str(f0)+')^(IM1))'
exp='IM0/(('+str(imagefreq)+'/'+str(f0)+')^(IM1))'

#Do maths on image
print 'Evaluating expression %s'%exp
print 'IM0 is %s'%model
print 'IM1 is %s'%spec_index 
immath(imagename=[model,spec_index],mode='evalexpr',expr=exp,
       outfile=outname)
#===============================================================================
# exp='IM2*IM0/((150000000/'+str(f0)+')^(IM1))' #beam*model/(150MHz/freq)^spec_index
# if os.path.exists(outname) and clobber:
#	rmtables(outname)
# immath(imagename=[model,spec_index,beam],mode='evalexpr',expr=exp,
#       outfile=outname)
#===============================================================================
print 'sectime2', sectime
imhead(outname,mode='put',hdkey='crval3',hdvalue=freq)
imhead(outname,mode='put',hdkey='cdelt3',hdvalue=bandwidth)
imhead(outname,mode='put',hdkey='date-obs',hdvalue=time_utc) # qa.time(sectime,form=["fits"])[0]) # use ymd or fits format
imhead(outname,mode='put',hdkey='crval4',hdvalue='I')

# Generate 2nd Taylor term
if os.path.exists(outnt2) and clobber:
    rmtables(outnt2)

if spectral_beam:
# Generate spectral image of the beam
  exp='log(IM0/IM1)/log('+str(f0-df/2)+'/'+str(f0+df/2)+')'
  beam_spec=calibrator+'_'+startfreq+'--'+endfreq+'_beamXX.im'
  immath(imagename=[calibrator+'_'+startfreq+'_beamXX.im',calibrator+'_'+endfreq+'_beamXX.im'],
        mode='evalexpr',expr=exp,outfile=beam_spec)

  immath(imagename=[outname,beam,spec_index,beam_spec],mode='evalexpr',outfile=outnt2, 
       expr='(IM0*IM1*(IM2+IM3))')
else:
  #=============================================================================
  # immath(imagename=[outname,beam,spec_index],mode='evalexpr',outfile=outnt2, 
  #     expr='(IM0*IM1*IM2)')
  #=============================================================================
    immath(imagename=[outname,spec_index],mode='evalexpr',outfile=outnt2, 
       expr='(IM0*IM1)')

imhead(outnt2,mode='put',hdkey='crval3',hdvalue=freq)
imhead(outnt2,mode='put',hdkey='cdelt3',hdvalue=bandwidth)
imhead(outnt2,mode='put',hdkey='date-obs',hdvalue=time_utc) # qa.time(sectime,form=["ymd"])[0])
imhead(outnt2,mode='put',hdkey='crval4',hdvalue='I')

#Export so we can easily see what our models were
exportfits(fitsimage=outname+'.fits',imagename=outname,overwrite=clobber)
exportfits(fitsimage=outnt2+'.fits',imagename=outnt2,overwrite=clobber)

# print 'Fourier transforming model...' #This applies it to the model column of the ms
# ft(vis=vis,model=[outname,outnt2],nterms=2,usescratch=True)

simobserve( project='HydA_154MHz',
            skymodel='HydA_154.88MHz.im',
            totaltime='1s',integration='1s',
            antennalist='mwa_ant.cfg',
            thermalnoise='tsys-manual', t_sky=200, tau0=0.00,
            obsmode='int',
            incell='2.2arcmin', # for 154 MHz , 3000 km baseline 
            mapsize=['3deg','3deg']
           )