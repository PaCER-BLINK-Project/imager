# test :

# import pdb

import pyuvdata
import numpy
import pickle
import sys
import os
import errno
import math
# import matplotlib.pyplot as plt
import time
import datetime
import h5py
from scipy import stats
import astropy.io.fits as pyfits
from scipy.io import loadmat

import os
import mwaconfig
from optparse import OptionParser,OptionGroup,Option

# to live with old values :
# astropy.utils.iers.iers.IERSRangeError: (some) times are outside of range covered by IERS table.
from astropy.utils.iers import conf
conf.auto_max_age = None


g_debug_level = 0
b_save_with_respect_to_ant0 = True # True is default  , REALLY after addition of --reference_antenna option it can be more than just ant0

def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass


def parse_options(idx):
   usage="Usage: %prog UVFITS SIMUL_FILE(.hdf5 or .mat) [options]\n"
#   usage+='\tASE_FILE.txt format - 2 columns Freq[MHz] and Phase[deg]\n'
   parser = OptionParser(usage=usage,version=1.00)
   parser.add_option('--transpose',action="store_true",dest="transpose",default=False, help="Transpose correlation matrix before saving ? [default %]")
   parser.add_option('--ignore_missing',action="store_true",dest="ignore_missing",default=False, help="Ignore missing UVFITS file [default %]")
   parser.add_option('--outdir','-o',dest="outdir",default="./",help="Output directory [default %default]")
   parser.add_option('--refant','--ref_ant','--reference_antenna','--reference_ant','-r',dest="reference_antenna",default=0,help="Output directory [default %default]",type="int")
   parser.add_option('--config','--config_file','--cfg','--instr_file','--instr_config',dest="config_file",default="instr_config_eda2.txt",help="Config file (with flags) [default %default]")
   parser.add_option('--channel','--ch','--dump_channel','--dump_ch','-c',dest="dump_channel",default=-1,help="Specify channel to dump, <0 -> average visibility over all channels [default %default]",type="int")
   parser.add_option('--dump_only','--dump',action="store_true",dest="dump_only",default=False, help="Only dump data and exit (no calibration) [default %]")
   
   (options, args) = parser.parse_args(sys.argv[idx:])


   return (options, args)


def calc_corr_matrix( uvfitsfile = "chan_204_20190818T120912.uvfits" , single_channel=-1 , n_ants=256, n_chan=32, n_pols=2, skip_ch=4, ignore_missing=False , outdir="./", transpose=False ) :
   global g_debug_level

   phase_vis_fits = os.path.basename(uvfitsfile).replace(".uvfits","_corr_phase.fits" )

#   n_ants=256
   n_inputs = n_ants*2
   n_baselines=((n_inputs*(n_inputs-1))/2)
#   n_chan=32
#   n_pols=2
#   skip_ch=4

   pol=0

   UV = pyuvdata.UVData() 
   corr_matrix=numpy.zeros( (n_ants,n_ants) , dtype=numpy.complex64 )
   
   if os.path.exists( uvfitsfile ) :
      UV.read( uvfitsfile, file_type='uvfits') 

      hdu_phases = pyfits.PrimaryHDU()
      hdu_phases.data = numpy.zeros( (n_ants,n_ants))


      for ant1 in range(0,n_ants):
         for ant2 in range(ant1,n_ants):
             bl = UV.antnums_to_baseline(ant1, ant2) 
             bl_ind = numpy.where(UV.baseline_array == bl)[0]
       
             mean_vis =  0
             count    = 0
             
             start_range = skip_ch
             stop_range  = (n_chan-skip_ch)
             if single_channel >= 0 :
                start_range = single_channel
                stop_range  = start_range + 1
                print("DEBUG : single channel requested -> forced channel range to be %d - %d" % (start_range,stop_range))


             # if there is only one channel just use this channel
             if UV.data_array.shape[2] == 1 :
                start_range = 0
                stop_range  = 1
             
             for ch in range(start_range,stop_range) :
                 if single_channel<0 or ch == single_channel :
                     l = len(UV.data_array[bl_ind, 0, ch, pol])
                     if l > 0 :            
                        mean_vis += UV.data_array[bl_ind, 0, ch, pol][0]
                        count += 1
                     else :
                        print "WARNING : 0 values for antenna1 = %d, antenna2 = %d , ch = %d , pol = %d (bl_ind=%s)" % (ant1,ant2,ch,pol,bl_ind)

             if count > 0 :       
                mean_vis = mean_vis / count
                if g_debug_level > 0 :
                    print "Baseline %d-%d visibility calculated from %d channels" % (ant1,ant2,count)
             else :
                print "WARNING : 0 values for baseline %d-%d" % (ant1,ant2)
       
             corr_matrix[ant1,ant2] = mean_vis
             hdu_phases.data[ant1,ant2] = numpy.angle(corr_matrix[ant1,ant2])*(180.00/math.pi)
             if ant1 != ant2 :
                corr_matrix[ant2,ant1] = mean_vis.conjugate()             
                hdu_phases.data[ant2,ant1] = numpy.angle(corr_matrix[ant2,ant1])*(180.00/math.pi) 

      hdulist = pyfits.HDUList([hdu_phases])
      full_path=outdir + "/" + phase_vis_fits
      print "Saving file %s" % (full_path)
      hdulist.writeto( full_path,overwrite=True)

             
      if g_debug_level > 1 : 
         for ant1 in range(0,n_ants):
            line = ""
            for ant2 in range(0,n_ants):
               if ant1<16 and ant2<16 :            
                  re = corr_matrix[ant1,ant2].real
                  im = corr_matrix[ant1,ant2].imag
           
                  line +=  ("%06.2f+j%06.2f " % (re,im))
            
            print line      
   else :
      if ignore_missing :
         print "WARNING : uvfits file %s does not exist, trying to continue ..." % (uvfitsfile)
         corr_matrix=numpy.ones( (n_ants,n_ants) , dtype=numpy.complex64 )
      else :
         print "ERROR : uvfits file %s does not exist !!! -> cannot continue !" % (uvfitsfile)
         os.sys.exit(-1)
         
             
   return corr_matrix             
   
def save_object( obj , name , outdir="./" ):
    pickle_file = name
    if name.find(".pkl") < 0 :
        pickle_file = name + '.pkl'

    pickle_file = outdir + pickle_file
    print "Saving file %s" % (pickle_file)


    with open( pickle_file, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        
    print "Saved object to pickle file %s" % (pickle_file)

def load_object( pickle_file, show=True ) :
   obj = None
 
   with open( pickle_file, 'rb') as f:
       obj = pickle.load(f)
       
   if show :
      print_coeff( obj, pickle_file )

   return obj


def save_matrix( matrix , uvfitsfile, output_postfix, outdir="./" , transpose=False) :
   # save all matrices for inspection or further optimisation and solving of gain parameters :
   postfix_re = "_" + output_postfix + "_RE"
   postfix_im = "_" + output_postfix + "_IM"   
   postfix_phase = "_" + output_postfix + "_PHASE"
   postfix_mag   = "_" + output_postfix + "_MAG"

   outfile_txt_re = uvfitsfile.replace(".uvfits", postfix_re + ".txt" )
   outfile_txt_im = uvfitsfile.replace(".uvfits", postfix_im + ".txt" )
   outfile_txt_phase = uvfitsfile.replace(".uvfits", postfix_phase + ".txt" )
   
   outfile_fits_re = uvfitsfile.replace(".uvfits", postfix_re + ".fits" )
   outfile_fits_im = uvfitsfile.replace(".uvfits", postfix_im + ".fits" )
   outfile_fits_phase = uvfitsfile.replace(".uvfits", postfix_phase + ".fits" )
   outfile_fits_mag = uvfitsfile.replace(".uvfits", postfix_mag + ".fits" )
   
   hdu_outfits_re = pyfits.PrimaryHDU()
   hdu_outfits_re.data = numpy.zeros( (n_ants,n_ants) )
   hdu_outfits_im = pyfits.PrimaryHDU()
   hdu_outfits_im.data = numpy.zeros( (n_ants,n_ants) )
   hdu_outfits_phase = pyfits.PrimaryHDU()
   hdu_outfits_phase.data = numpy.zeros( (n_ants,n_ants) )
   hdu_outfits_mag = pyfits.PrimaryHDU()
   hdu_outfits_mag.data = numpy.zeros( (n_ants,n_ants) )
   
   full_path=outdir + "/" + outfile_txt_re
   print "Saving file %s" % (full_path)
   out_f_txt_re = open( full_path , "w" )
   
   full_path=outdir + "/" + outfile_txt_im
   print "Saving file %s" % (full_path)
   out_f_txt_im = open( full_path , "w" )

   full_path=outdir + "/" + outfile_txt_im
   print "Saving file %s" % (full_path)   
   out_f_txt_phase = open( full_path , "w" )
   
   # TODO : loop to save 
   for ant1 in range(0,n_ants) :
      line_matrix_txt_re = ""   # corr_matrix RE
      line_matrix_txt_im = ""   # corr_matrix IM 
      line_matrix_txt_phase = "" 
   
      for ant2 in range(0,n_ants) :
          vis_obs   = matrix[ant1,ant2]
          vis_phase_deg = numpy.angle(vis_obs)*(180.00/math.pi)

          line_matrix_txt_re += ("%+07.2f " % vis_obs.real)
          line_matrix_txt_im += ("%+07.2f " % vis_obs.imag)
          line_matrix_txt_phase += ("%+07.2f " % vis_phase_deg )
          
          
          hdu_outfits_re.data[ant1,ant2] = vis_obs.real
          hdu_outfits_im.data[ant1,ant2] = vis_obs.imag
          hdu_outfits_phase.data[ant1,ant2] = vis_phase_deg
          hdu_outfits_mag.data[ant1,ant2] = vis_obs.real*vis_obs.real + vis_obs.imag*vis_obs.imag 
      
      out_f_txt_re.write( line_matrix_txt_re + "\n" )
      out_f_txt_im.write( line_matrix_txt_im + "\n" )
      out_f_txt_phase.write( line_matrix_txt_phase + "\n" )

   
   out_f_txt_re.close()
   out_f_txt_im.close()
   out_f_txt_phase.close()
   
   if transpose :
      hdu_outfits_re.data = hdu_outfits_re.data.transpose()
      hdu_outfits_im.data = hdu_outfits_im.data.transpose()
      hdu_outfits_phase.data = hdu_outfits_phase.data.transpose()
      hdu_outfits_mag.data = hdu_outfits_mag.data.transpose()
   

   hdulist_re = pyfits.HDUList([hdu_outfits_re])
   full_path=outdir + "/" + outfile_fits_re
   print "Saving file %s" % (full_path)   
   hdulist_re.writeto( full_path ,overwrite=True)

   full_path=outdir + "/" + outfile_fits_im
   print "Saving file %s" % (full_path)   
   hdulist_im = pyfits.HDUList([hdu_outfits_im])
   hdulist_im.writeto( full_path ,overwrite=True)

   hdulist_phase = pyfits.HDUList([hdu_outfits_phase])
   full_path=outdir + "/" + outfile_fits_phase
   print "Saving file %s" % (full_path)   
   hdulist_phase.writeto( full_path ,overwrite=True)

   hdulist_mag = pyfits.HDUList([hdu_outfits_mag])
   full_path=outdir + "/" + outfile_fits_mag
   print "Saving file %s" % (full_path)   
   hdulist_mag.writeto( full_path ,overwrite=True)


if __name__ == "__main__":
   n_ants=256

   uvfitsfile = "chan_204_20190818T120912.uvfits"
   if len(sys.argv) > 1:   
       uvfitsfile = sys.argv[1]
   base_uvfitsfile = os.path.basename(uvfitsfile)
   utc = base_uvfitsfile[9:24]       
   uxtime = time.mktime(datetime.datetime.strptime(utc, "%Y%m%dT%H%M%S").timetuple()) + 8*3600 # just for Perth !!!

   simul_hdf5_file = None # "UTC20190809_125000EDA_X.hdf5"
   if len(sys.argv) > 2 and sys.argv[2] != "-" :
       simul_hdf5_file = sys.argv[2]   

   (options, args) = parse_options(2)
   
   print "#########################################################"
   print "PARAMETERS :"
   print "#########################################################"
   print "UV-fits = %s -> utc = %s -> %.2f (base_uvfits = %s)" % (uvfitsfile,utc,uxtime,base_uvfitsfile)
   print "Outdir = %s" % (options.outdir)   
   print "Reference antenna = %d" % (options.reference_antenna)
   print "Config file = %s" % (options.config_file)
   print("Dump single channel = %d" % (options.dump_channel))
   print("Dump data only      = %s" % (options.dump_only))
   print "#########################################################"
   
   mkdir_p( options.outdir )
       
   hdf5_file=True
   if simul_hdf5_file is not None :
       if simul_hdf5_file.find(".mat") >= 0 :
           hdf5_file=False

   config_file = options.config_file # "instr_config_eda2.txt"
   config = None
   flagged_antenna_list = []
   if os.path.exists( config_file ) :
       config = mwaconfig.MWAConfig( config_file )       
       config_inputs = config.read_mapping_file( config_file )
       print "INFO : read %d inputs from config file %s OK" % (len(config_inputs),config_file)
       
       flagged_antenna_list = config.GetFlaggedAntList() 
       if flagged_antenna_list is None :
          print "\t !!! WARNING : no flagged antennas specified in the config file %s - really, please verify !!!" % (config_file)
       else :
          print "\tThere are %d flagged antennas = %s" % (len(flagged_antenna_list),flagged_antenna_list)
       
#       time.sleep(5)
   # read config file if there is one :
   
    
   
   outfile    = base_uvfitsfile.replace(".uvfits",".txt") 
   out_eigen_vector1 =  base_uvfitsfile.replace(".uvfits","_EigenVector1.txt")
   corr_matrix_pickle = base_uvfitsfile.replace(".uvfits","_CorrMatrix.pkl")
   

   corr_matrix = calc_corr_matrix( uvfitsfile=uvfitsfile , single_channel=options.dump_channel , n_ants=n_ants, ignore_missing=options.ignore_missing, outdir=options.outdir )           
   (lamb,eigen_vectors) = numpy.linalg.eig( corr_matrix )   
   # lamb - eigen values
   # eigen_vectors
   save_object( corr_matrix, corr_matrix_pickle, options.outdir )
   
   # moved from the very end here:
   save_matrix( corr_matrix             , base_uvfitsfile, "CorrMatrix" , options.outdir, transpose=options.transpose )
   
   if options.dump_only :
      print("INFO : only data dump is required -> exiting now")
      sys.exit(0)

   corr_matrix_simul = numpy.ones( (n_ants,n_ants) , dtype=numpy.complex64 )
   lamb_simul        = numpy.zeros( n_ants )
   eigen_vectors_simul = [] 
   for i in range(0,n_ants) :
       eigen_vectors_simul.append( numpy.zeros( n_ants, dtype=numpy.complex64 ) )

   # Initialise gains with ones :
   g_abs_sq            = numpy.ones( n_ants, dtype=numpy.complex64 ) # None # |g|^2
   g_abs               = numpy.ones( n_ants, dtype=numpy.complex64 ) # None # |g|    
   g_abs_r             = numpy.ones( n_ants, dtype=numpy.complex64 ) # Normalised gains = gi / |g1|
   r_phase_ev          = numpy.zeros( n_ants )
   if simul_hdf5_file is not None :
       corr_matrix_simul = None
       corr_matrix_simul_phase = None
       if hdf5_file :               
           DATA = h5py.File( simul_hdf5_file , 'r' )
           corr_matrix_simul = numpy.array( DATA["C"] )
           print "Read HDF5 file %s correctly !" % (simul_hdf5_file)
       else :
           matlab_simul_file = loadmat( simul_hdf5_file )               
#           corr_matrix_simul_phase = matlab_simul_file['POI_new']
           # corr_matrix_simul = numpy.cos(corr_matrix_simul_phase)+1j*numpy.sin(corr_matrix_simul_phase)
           corr_matrix_simul = matlab_simul_file['POI'] # was POI_new
           
           # test :
           tmp = corr_matrix_simul.conjugate()
           tmp = corr_matrix_simul.transpose()
           corr_matrix_simul = tmp.copy()
           print "Transpose test"
#           time.sleep(10)
           
       (lamb_simul,eigen_vectors_simul) = numpy.linalg.eig( corr_matrix_simul ) 

       # if simulation file provided -> initialise gains with zeros :              
       g_abs_sq = numpy.zeros( n_ants, dtype=numpy.complex64 )
       g_abs    = numpy.zeros( n_ants, dtype=numpy.complex64 )
       g_abs_r  = numpy.zeros( n_ants, dtype=numpy.complex64 )
   else :
       print "WARNING : no hdf5 simulation file provided - using just ONEs (Galactic Transit phase calibration - no gain calibration)"              
#       time.sleep(1)

   # gains from auto-correlations 
#   (g_abs_sq,g_abs) = calc_gains( corr_matrix, corr_matrix_simul, n_ants, g_abs_sq, g_abs )
   full_path=options.outdir + "/" + outfile
   print "Saving file %s" % (full_path)   
   out_f = open( full_path,"w")         
   header = "# Index Corr_matrix_diagonal_DATA  EigenValues_DATA Corr_matrix_diagonal_SIMUL EigenValues_SIMUL GAIN^2 GAIN GAIN_r"
   print "%s" % (header)
   out_f.write( header + "\n" )   
   n_lamb = len(lamb)
   for l in range(0,n_lamb) :
      g_abs_sq[l] = corr_matrix[l,l] / corr_matrix_simul[l,l]
      g_abs[l]    = math.sqrt( g_abs_sq[l] )
      g_abs_r[l]  = g_abs[l] / g_abs[0]
   
      if lamb_simul is not None and eigen_vectors_simul is not None :
          line = ( "%d %s %s %s %s %s %s %s" % (l,corr_matrix[l,l],lamb[l],corr_matrix_simul[l,l],lamb_simul[l], g_abs_sq[l], g_abs[l] , g_abs_r[l]  ) )
      else :
          line = ( "%d %s %s" % (l,corr_matrix[l,l],lamb[l]) )
      print "%s" % (line)
      out_f.write( line + "\n" )               
   out_f.close()
   
   # save eigen vectos 1 
   full_path = options.outdir + "/" + out_eigen_vector1
   print "Saving file %s" % (full_path)   
   out_f = open( full_path , "w" )
   header =  "# ANTENNA EigenVector_DATA EigenVector_SIMUL DATA/SIMUL MAG PHASE[deg]"
   out_f.write( header + "\n" )
   for l in range(0,n_lamb) :
      if lamb_simul is not None and eigen_vectors_simul is not None :      
          evec_data  = eigen_vectors[0]
          evec_simul = eigen_vectors_simul[0]
          
          r = evec_data[l] / evec_simul[l]
          r_mag = abs(r)**2
          r_phase = numpy.angle(r)*(180.00/math.pi) 
          r_phase_ev[l] = r_phase
          
          line = "%d %s %s %s %.8f %.2f" % (l,evec_data[l],evec_simul[l],r,r_mag,r_phase)
      else :
          print "NOT IMPLEMENTED !!!\n"
      
      out_f.write( line + "\n" )          
   out_f.close()

   # solve for phases using gains calculated from AUTOs :
   hdu_simul_phases = pyfits.PrimaryHDU()
   hdu_simul_phases.data = numpy.zeros( (n_ants,n_ants))

   g_reference_antenna = options.reference_antenna # global variable for a reference antenna 
   corr_matrix_simul_solve = corr_matrix_simul.copy()
   for r in range(0,n_ants) :
      for c in range(0,n_ants) :
         hdu_simul_phases.data[r,c] = numpy.angle( corr_matrix_simul[r,c] )*(180.00/math.pi)  
      
         if r == c :
            # AUTO-CORR :
            corr_matrix_simul_solve[r,c] = corr_matrix_simul[r,c] * g_abs_sq[r]
         else :
            # CROSS-CORR :
            corr_matrix_simul_solve[r,c] = g_abs[r] * corr_matrix_simul[r,c] * g_abs[c] # conjugate omitted as phase not yet included 
            
   hdulist = pyfits.HDUList([hdu_simul_phases])
   corr_matrix_simul_fits = "simul_phase.fits"
   if simul_hdf5_file is not None :
       if hdf5_file :
          corr_matrix_simul_fits = os.path.basename( simul_hdf5_file ).replace( ".hdf5", "_phase.fits" )
       else :
          corr_matrix_simul_fits = os.path.basename( simul_hdf5_file ).replace( ".mat", "_phase.fits" )
          
   full_path = options.outdir + "/" + corr_matrix_simul_fits
   print "Saving file %s" % (full_path)   
          
   hdulist.writeto( full_path ,overwrite=True)


# TODO : use each antenna as reference, subtract phase of antenna =0 to get it back to referenced to zero - I will have N-1 phases for a given antenna and 
#        I can later average or find median to get it back to have one value !


   # First just print equations with reference antenna included :
   eq = 1
   eq_filename = "eq_%04d.txt" % (eq)

   full_path = options.outdir + "/" + eq_filename
   print "Saving file %s" % (full_path)   
   out_eq_f = open( full_path , "w" )
   for r in range(0,n_ants) :
      for c in range(0,n_ants) :
         # autos skipped (not phase info only gain already used : 
         if r == g_reference_antenna or c == g_reference_antenna :
             ratio = ( corr_matrix[r,c] / corr_matrix_simul_solve[r,c] )
             r_phase = numpy.angle(ratio)*(180.00/math.pi)
             vis = corr_matrix[r,c]                 
             vis_simul = corr_matrix_simul_solve[r,c]

#             if c == g_reference_antenna :                    
#                line = ( "%d : (%06.2f,%06.2f) = (%06.2f,%06.2f) * exp( - i Phi_%03d ) -> phase = %.2f [deg] vs. %.2f [deg]" % (r,vis.real,vis.imag,vis_simul.real,vis_simul.imag,r,r_phase,r_phase_ev[r]) )
#                out_eq_f.write( line + "\n" )

# Do not save complex conjugate equation :
             if r == g_reference_antenna :
                 line = ( "%d : %s = %s * exp( + i Phi_%03d ) -> phase = %.2f [deg] vs. %.2f [deg]" % (c,corr_matrix[r,c],corr_matrix_simul_solve[r,c],c,r_phase,r_phase_ev[c]) )
                 out_eq_f.write( line + "\n" )
                 
   out_eq_f.close()
   
   # full loop over reference antennas :
   # solve for phases :   
   cal_sol_per_ref_ant = {}
   for ref_antenna in range(0,n_ants) : # loop over reference antennas to solve equations (note ref_antenna != g_reference_antenna) 
       cal_sol_per_ref_ant[ref_antenna] = numpy.zeros(n_ants)
# TODO : use each antenna as reference, subtract phase of antenna =0 to get it back to referenced to zero - I will have N-1 phases for a given antenna and 
#        I can later average or find median to get it back to have one value !
# print equations with reference antenna included :
       out_filename = "eq_refant%04d.txt" % (ref_antenna)
       full_path = options.outdir + "/" + out_filename
       print "Saving file %s" % (full_path)   

       out_f = open(  full_path , "w" )
       for r in range(0,n_ants) :
           for c in range(0,n_ants) :
               # autos skipped (not phase info only gain already used : 
               if r == ref_antenna or c == ref_antenna :
                   ratio = ( corr_matrix[r,c] / corr_matrix_simul_solve[r,c] )
                   r_phase = numpy.angle(ratio)*(180.00/math.pi)
                   vis = corr_matrix[r,c] 
                   vis_simul = corr_matrix_simul_solve[r,c]

# DO NOT SAVE COMPLEX CONJUGATE EQUATION :
#                   if c == ref_antenna :                    
#                       line = ( "%d : (%06.2f,%06.2f) = (%06.2f,%06.2f) * exp( - i Phi_%03d ) -> phase = %.2f [deg] vs. %.2f [deg]" % (r,vis.real,vis.imag,vis_simul.real,vis_simul.imag,r,r_phase,r_phase_ev[r]) )
#                       out_f.write( line + "\n" )
                           
                       # fill output numpy array                            
#                       cal_sol_per_ref_ant[ref_antenna][r] = r_phase

# DO NOT SAVE COMPLEX CONJUGATE EQUATION :                           
                   if r == ref_antenna :
                       line = ( "%d : %s = %s * exp( + i Phi_%03d ) -> phase = %.2f [deg] vs. %.2f [deg]" % (c,corr_matrix[r,c],corr_matrix_simul_solve[r,c],c,r_phase,r_phase_ev[c]) )
                       out_f.write( line + "\n" )
                       
                       cal_sol_per_ref_ant[ref_antenna][c] = r_phase
                 
       out_f.close()
       
       
   # save calibration solutions referenced to a specified reference_antenna - there was a lot of confusion due 
   # to having the same variable name ref_antenna in loops and globally ...    
   full_path = options.outdir + "/" + "calsolphase.txt"
   print "Saving file %s" % (full_path)   
       
   out_f_test = open( full_path , "w" )
   for ref_antenna in range(0,n_ants) :       
      calsol_phase_ant0 = cal_sol_per_ref_ant[ref_antenna][g_reference_antenna] # 0 -> g_reference_antenna in order to be able to save with respect to a specified reference antenna 
      print "DEBUG : ref_antenna = %d -> calsol_phase_ant0 = %.4f" % (ref_antenna,calsol_phase_ant0)
      
      line = "%03d : " % (ref_antenna)
      for a in range(0,len(cal_sol_per_ref_ant[ref_antenna])) :
          if b_save_with_respect_to_ant0 : # REALLY after addition of --reference_antenna option it can be more than just ant0 (in fact any other antenna)
             cal_sol_per_ref_ant[ref_antenna][a] = cal_sol_per_ref_ant[ref_antenna][a] - calsol_phase_ant0
             
             # normalise to -180,180 range 
             if cal_sol_per_ref_ant[ref_antenna][a] > 180.00 :
                 cal_sol_per_ref_ant[ref_antenna][a] = cal_sol_per_ref_ant[ref_antenna][a] - 360.00
             elif cal_sol_per_ref_ant[ref_antenna][a] < -180.00 :
                 cal_sol_per_ref_ant[ref_antenna][a] = cal_sol_per_ref_ant[ref_antenna][a] + 360.00
          
          # is flagged -> set NaN :       
          if a in flagged_antenna_list :
             cal_sol_per_ref_ant[ref_antenna][a] = numpy.NaN
                 
                 
          line += ('%+07.2f ' % (cal_sol_per_ref_ant[ref_antenna][a]))

      print "\tDEBUG : %s" % (line)      
      out_f_test.write( line + "\n" )
   out_f_test.close()


   # for each antenna calculate mean over equations solved using different refernece antennas :
   full_path = options.outdir + "/" + "calsolphase_mean.txt"
   print "Saving file %s" % (full_path)   
   out_f_test_mean = open( full_path , "w" )
   # IQR : https://docs.scipy.org/doc/scipy-0.19.0/reference/generated/scipy.stats.iqr.html
   out_f_test_mean.write("# ANT MEAN(over-ref-ant) STDDEV(over-ref-ant) MEDIAN(over-ref-ant) IQR/1.349 |Gain| Flagged UxTime\n") 
   for a in range(0,n_ants) :
       ant_cal = []
       for ref_antenna in range(0,n_ants) :
           if not numpy.isnan( cal_sol_per_ref_ant[ref_antenna][a] ) and ( math.fabs( cal_sol_per_ref_ant[ref_antenna][a] )>0.00001  or a==g_reference_antenna ): # ABS > 0.00001 or reference antenna 
               ant_cal.append( cal_sol_per_ref_ant[ref_antenna][a] )
       
       mean = numpy.nanmean(ant_cal)
       median = numpy.nanmedian(ant_cal)
       rms  = numpy.nanstd(ant_cal)
       iqr = stats.iqr(ant_cal,nan_policy='omit',scale='normal')
       
#       if numpy.isnan(mean) :
#          mean = 0
#       if numpy.isnan(median) :
#          median = 0
#       if numpy.isnan(rms) :
#          rms = 0
#       if numpy.isnan(iqr) :
#          iqr = 0
       
       flagged = 0 
       if a in flagged_antenna_list :
          flagged = 1   
       
       line = ( "%03d %+07.2f %+07.2f %+07.2f %+07.2f %+12.8f %d %.1f\n" % (a,mean,rms,median,iqr,g_abs[a],flagged,uxtime) )
       out_f_test_mean.write( line )
       
   out_f_test_mean.close()

   
   # save all matrices for inspection or further optimisation and solving of gain parameters :
   # def save_matrix( matrix , uvfitsfile, output_postfix ) :
   save_matrix( corr_matrix_simul       , base_uvfitsfile, "SimulMatrix", options.outdir )
   save_matrix( corr_matrix_simul_solve , base_uvfitsfile, "SimulTimesGainMatrix" , options.outdir )
   
#   outfile_corr_matrix_txt_re = uvfitsfile.replace(".uvfits","_CorrMatrix_RE.txt")
#   outfile_corr_matrix_txt_im = uvfitsfile.replace(".uvfits","_CorrMatrix_IM.txt")
#   outfile_corr_matrix_simul_re = uvfitsfile.replace(".uvfits","_SimulMatrix_RE.txt")
#   outfile_corr_matrix_simul_im = uvfitsfile.replace(".uvfits","_SimulMatrix_IM.txt")
#   outfile_corr_matrix_simulTimesGain_re = uvfitsfile.replace(".uvfits","_SimulTimesGainMatrix_RE.txt")
#   outfile_corr_matrix_simulTimesGain_im = uvfitsfile.replace(".uvfits","_SimulTimesGainMatrix_IM.txt")
#   
#   outfile_corr_matrix_fits_re = uvfitsfile.replace(".uvfits","_CorrMatrix_RE.fits")
#   outfile_corr_matrix_fits_im = uvfitsfile.replace(".uvfits","_CorrMatrix_IM.fits")
#   outfile_corr_matrix_simul_fits_re = uvfitsfile.replace(".uvfits","_SimulMatrix_RE.fits")
#   outfile_corr_matrix_simul_fits_im = uvfitsfile.replace(".uvfits","_SimulMatrix_IM.fits")
#   outfile_corr_matrix_simulTimesGain_fits_re = uvfitsfile.replace(".uvfits","_SimulTimesGainMatrix_RE.fits")
#   outfile_corr_matrix_simulTimesGain_fits_im = uvfitsfile.replace(".uvfits","_SimulTimesGainMatrix_IM.fits")
#
#   out_f_corr_matrix_txt_re = open( outfile_corr_matrix_txt_re , "w" )
#   out_f_corr_matrix_txt_im = open( outfile_corr_matrix_txt_im , "w" )
#   out_f_corr_matrix_simul_re = open( outfile_corr_matrix_simul_re , "w" )
#   out_f_corr_matrix_simul_im = open( outfile_corr_matrix_simul_im , "w" )
#   out_f_corr_matrix_simulTimesGain_re = open( outfile_corr_matrix_simulTimesGain_re , "w" )
#   out_f_corr_matrix_simulTimesGain_im = open( outfile_corr_matrix_simulTimesGain_im , "w" )
#   
#   # TODO : loop to save 
#   for ant1 in range(0,n_ants) :
#      line_corr_matrix_txt_re = ""   # corr_matrix RE
#      line_corr_matrix_txt_im = ""   # corr_matrix IM 
#      line_corr_matrix_simul_re = "" # corr_matrix_simul RE
#      line_corr_matrix_simul_im = "" # corr_matrix_simul IM 
#      line_corr_matrix_simulTimesGain_re = "" # corr_matrix_simul_solve RE
#      line_corr_matrix_simulTimesGain_im = "" # corr_matrix_simul_solve IM
#   
#      for ant2 in range(0,n_ants) :
#          vis_obs = corr_matrix[ant1,ant2]
#          vis_simul = corr_matrix_simul[ant1,ant2]
#          vis_simul_solve = corr_matrix_simul_solve[ant1,ant2]
#          
#          line_corr_matrix_txt_re += ("%+07.2f " % vis_obs.real)
#          line_corr_matrix_txt_im += ("%+07.2f " % vis_obs.imag)     
#          line_corr_matrix_simul_re += ("%+07.2f " % vis_simul.real)
#          line_corr_matrix_simul_im += ("%+07.2f " % vis_simul.imag)
#          line_corr_matrix_simulTimesGain_re += ("%+07.2f " % vis_simul_solve.real)
#          line_corr_matrix_simulTimesGain_im += ("%+07.2f " % vis_simul_solve.imag)
#      
#      out_f_corr_matrix_txt_re.write( line_corr_matrix_txt_re + "\n" )
#      out_f_corr_matrix_txt_im.write( line_corr_matrix_txt_im + "\n" )
#      out_f_corr_matrix_simul_re.write( line_corr_matrix_simul_re + "\n" )
#      out_f_corr_matrix_simul_im.write( line_corr_matrix_simul_im + "\n" )
#      out_f_corr_matrix_simulTimesGain_re.write(  line_corr_matrix_simulTimesGain_re + "\n" )
#      out_f_corr_matrix_simulTimesGain_im.write(  line_corr_matrix_simulTimesGain_im + "\n" )
#   
#   out_f_corr_matrix_txt_re.close()
#   out_f_corr_matrix_txt_im.close()
#   out_f_corr_matrix_simul_re.close()
#   out_f_corr_matrix_simul_im.close()
#   out_f_corr_matrix_simulTimesGain_re.close()
#   out_f_corr_matrix_simulTimesGain_im.close()
   
   
   

   
       

   # solve using Adrian's eigen values methods :
   # Step 1/ : Rv' = Rv * |g1|^2 ( Eq. 1 ) = corr_matrix_simul * |g1|^2 ( Eq. 1 )
   # Step 1a : Gr = just a matrix of gains divided by |g1| - normalised by magnitude of antenna 1 gain 
   # Step 2/ : find eigenvalues / unitary matrices Q_v for matrix Rv' = Qv * Lambda_v * Qv^H (Eq. 6 )
   # Step 3/ : find eigenvalues / unitary matrices Q_c for data matrix Rc = corr_matrix = Q_c * Lambda_c * Q_c^H (Eq. 5)
   # Step 4/ : Eq. 7 : G_r^-1 Rc (G_r^H)^-1 = Rv' =  G_r^-1 ( Q_c * Lambda_c * Q_c^H )  (G_r^H)^-1  = Qv * Lambda_v * Qv^H 
   #           Keep in mind that only magnitudes of G_r are known, but not phases !
   # Step 5/ : Assume that all antennas have the same magnitude of gain and get initial calibration solutions from the first eigenvectors (page 3 out of 9) :
   #           Gr^-1 Qc = Qv -> we can get gr2 , gr3 etc (with respect to ant1) by element by element division of the first eigenvector 
   # Step 6/ : Compare solutions to my other methods 
   # 
         
   
#   if True : 
##   plt.plot( numpy.abs(corr_matrix)**2 , label="test" )
#      plt.imshow( numpy.log10(numpy.abs(corr_matrix)**2), cmap='rainbow') # , vmin=mean-5*rms, vmax=mean+5*rms)
#      plt.colorbar()
#      plt.show()

   