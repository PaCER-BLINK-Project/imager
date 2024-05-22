"""
Script to dump visibilities into a correlation matrix FITS file 

Input: CASA ms 
Output: FITS files for REAL and IMAG part of visibilities 
"""

import sys,optparse, os
import re, math
import numpy
from string import Template
import pyfits

def write_corr_matrix_line( vis_file_re_f, line, line_count ) :
   line_str = ""
   for i in range(0,line_count) :
      vis = line[i]
      line_str += ("%.4f " % numpy.real(vis) )   
      
   vis_file_re_f.write( line_str + "\n" )

def parse_options():
   # create command line arguments
   parser=optparse.OptionParser()
   parser.set_usage("Usage: casapy -c get_corrmatrix_from_casa.py CASA.ms")
   parser.add_option('-o','--obsid',dest="obsid",default="-1", help="Observation ID",metavar="INTEGER",type="int")
   parser.add_option('-c','--channel',dest="channel",default="0", help="Channel index to dump (<=0 -> all) [default %default]",metavar="INTEGER",type="int")
   parser.add_option('-t','--time','--time_index','--time_step','--timestamp',dest="time_step",default="0", help="Time step (or integration) to dump (<=0 -> all) [default %default]",metavar="INTEGER",type="int")
   parser.add_option('-m','--max','--max_rows','--max_rows',dest="max_rows",default="100000", help="Maximum number of rows to read [default %default]",metavar="INTEGER",type="int")
   
   parser.add_option('--data_column','--column',dest="data_column",default="CORRECTED_DATA", help="Data column [default %default]")
   parser.add_option('--save_full_matrix','--save_conjugates','--conjugate',dest="save_conjugates",action="store_true", default=False, help="Save all values including conjugates [default %default]")
   parser.add_option('--polarisation_index','--pol','--polarisation',dest="polarisation",default=0,help="Polarisation index 0=xx,1=xy,2=yx,3=yy [default %default]",type="int")
   
   # parser.add_option('-s','--source',dest="source",default="HydA", help="Calibrator name",metavar="STRING")
   # parser.add_option('-a','--azim',dest="azim",default=288.100147, help="Externally passed azimuth of the source",metavar="FLOAT",type=float)
   # parser.add_option('-z','--za',dest="za",default=34.925223, help="Externally passed zenithal distance of the source",metavar="FLOAT",type=float)
 
   # parse command line arguments
   casa_index=sys.argv.index('-c')+2  # remove CASA options
   (options,args)=parser.parse_args(sys.argv[casa_index:])
   
   return (options, args)
   
if __name__ == '__main__':         
   
   (options, args) = parse_options()
   casams=args[-1]

   print("####################################################################################")
   print("PARAMETERS:")
   print("####################################################################################")
   print("CASA measurement set = %s" % (casams))
   print("Time step            = %d" % (options.time_step))
   print("Channel              = %d" % (options.channel))
   print("Max rows             = %d" % (options.max_rows))
   print("Data column in CASA  = %s" % (options.data_column))
   print("Polarisation         = %d" % (options.polarisation))
   print("####################################################################################")

   tb.open(casams + "/ANTENNA")
   Ant_names = tb.getcol('NAME')
   Ant_positions = tb.getcol('POSITION')
   tb.close()
   n_ant = len(Ant_names)
   n_timesteps = (n_ant*(n_ant+1))/2
   n_expected_baselines = (n_ant*(n_ant+1))/2

   print("Shape = %d x %d , number of antennas = %d -> number of recrods per timestep = %d" % (Ant_positions.shape[0],Ant_positions.shape[1],n_ant,n_timesteps))

   for ant_idx in range(0,len(Ant_names)) :
#      position = Ant_positions[ant_idx]
      position_x = Ant_positions[0,ant_idx]
      position_y = Ant_positions[1,ant_idx]
      position_z = Ant_positions[2,ant_idx]
      print("Antenna %03d : %s at (%.8f,%.8f,%.8f)" % (ant_idx,Ant_names[ant_idx],position_x,position_y,position_z))

# sys.exit(-1)

   tb.open(casams,nomodify=True)
   col=options.data_column # 'CORRECTED_DATA' # use calibrated data !!!
   # col='DATA'
# data=tb.getcol(col) #This is pol x channel x antenna

   ant_data1=tb.getcol('ANTENNA1',nrow=options.max_rows)
   ant_data2=tb.getcol('ANTENNA2',nrow=options.max_rows)
   data = tb.getcol(col,nrow=options.max_rows)
   uvw  = tb.getcol("UVW",nrow=options.max_rows)

   print("DEBUG : data.shape = %d x %d x %d , uvw.shape = %d x %d" % (data.shape[0],data.shape[1],data.shape[2],uvw.shape[0],uvw.shape[1]))

   # outputs fits files:
   # Real :
   hdu_re = pyfits.PrimaryHDU()
   hdu_re.data = -numpy.ones( (n_ant,n_ant))*numpy.nan
   hdulist_re = pyfits.HDUList([hdu_re])
   # Imaginary :
   hdu_im = pyfits.PrimaryHDU()
   hdu_im.data = -numpy.ones( (n_ant,n_ant))*numpy.nan
   hdulist_im = pyfits.HDUList([hdu_im])
   
   # UVW FITS FILES :
   hdu_u = pyfits.PrimaryHDU()
   hdu_u.data = -numpy.ones( (n_ant,n_ant))*numpy.nan
   hdulist_u = pyfits.HDUList([hdu_u])
   hdu_v = pyfits.PrimaryHDU()
   hdu_v.data = -numpy.ones( (n_ant,n_ant))*numpy.nan
   hdulist_v = pyfits.HDUList([hdu_v])
   hdu_w = pyfits.PrimaryHDU()
   hdu_w.data = -numpy.ones( (n_ant,n_ant))*numpy.nan
   hdulist_w = pyfits.HDUList([hdu_w])
   # END OF UVW 


   postfix= ( 'channel%03d_time%06d_pol%02d' % (options.channel,options.time_step,options.polarisation))
   out_vis_refits=casams.replace('.ms', ('_vis_real_%s.fits' % postfix))
   out_vis_imfits=casams.replace('.ms', ('_vis_imag_%s.fits' % postfix))
   out_u_fits=casams.replace('.ms', ('_u_%s.fits' % postfix))
   out_v_fits=casams.replace('.ms', ('_v_%s.fits' % postfix))
   out_w_fits=casams.replace('.ms', ('_w_%s.fits' % postfix))

   out_text_file = casams.replace('.ms', ('_vis_%s.txt' % postfix))
   vis_file_re_f = open( out_text_file ,"w" )

   line=numpy.ones(n_ant)*numpy.nan
   line_count=0

   # TODO : skip options.time_step
#   if options.time_step > 0 :
#      print("ERROR : option to select some specific timestamp is yet to be implemented")
#      sys.exit(-1)

   # skip first options.time_step steps to save timestep specified in options :
   time_step=0
   i = 0   
   if options.time_step > 0 :
      if ant_data1[0] == ant_data2[0] and ant_data1[0] == 0 :
         print("DEBUG : first row for antenna pair (0,0) -> will skip %d first timesteps" % (options.time_step))
      else :
         print("ERROR : first row is NOT ANTENNA PAIR (0,0) -> cannot continue")
         sys.exit(-1)
      
      i = 1
      while time_step < options.time_step :
         curr_ant1 = ant_data1[i]
         curr_ant2 = ant_data2[i]
      
         if curr_ant1 == curr_ant2 and curr_ant1 == 0 :
            print("DEBUG : skipped timestep = %d (row = %d)" % (time_step,i))
            time_step += 1
            
            if time_step == options.time_step :
               print("DEBUG : skipped %d timesteps starting at index = %d" % (time_step,i))               
               break # break is not to increase row index by 1 to start at antenna pair (0,0)
      
         i = i + 1 
      
      print("DEBUG : skipped %d  timesteps , current row index = %d" % (time_step,i))
   else :
      print("DEBUG : starting from time_step = 0 -> no skipping of rows is required")      
   
   
   prev_ant1=-1
   prev_ant2=-1
   baseline_count=0
   all_baselines=0
   end=False
   for i1 in range(i,len(ant_data1)) : # starting from i which is the starting row of specified integration (see --time option)
      curr_ant1 = ant_data1[i1]
      curr_ant2 = ant_data2[i1]
      
      if curr_ant1 > (prev_ant1+1):
         print("WARNING : missed antenna1 = %d" % (prev_ant1+1))
      
      if curr_ant1 >= prev_ant1 : 
         # save to FITS file :
         flag=""
         if curr_ant1 == curr_ant2 :
            flag = "auto"
            
         print("DEBUG saving visibility for antenna pair %d - %d = %.4f (%s)" % (ant_data1[i1],ant_data2[i1],numpy.real( data[options.polarisation,options.channel,i1] ), flag ))
         hdu_re.data[ant_data1[i1]][ant_data2[i1]] = numpy.real( data[options.polarisation,options.channel,i1] )
         hdu_im.data[ant_data1[i1]][ant_data2[i1]] = numpy.imag( data[options.polarisation,options.channel,i1] )         
         hdu_u.data[ant_data1[i1]][ant_data2[i1]] = uvw[0,i1]
         hdu_v.data[ant_data1[i1]][ant_data2[i1]] = uvw[1,i1]
         hdu_w.data[ant_data1[i1]][ant_data2[i1]] = uvw[2,i1]
         
         # optionally save conjugated and also -UVW :
         if options.save_conjugates :
            # save conjugates too:
            hdu_re.data[ant_data2[i1]][ant_data1[i1]] = numpy.real( data[options.polarisation,options.channel,i1] )
            hdu_im.data[ant_data2[i1]][ant_data1[i1]] = -numpy.imag( data[options.polarisation,options.channel,i1] )
            hdu_u.data[ant_data2[i1]][ant_data1[i1]] = -uvw[0,i1]
            hdu_v.data[ant_data2[i1]][ant_data1[i1]] = -uvw[1,i1]
            hdu_w.data[ant_data2[i1]][ant_data1[i1]] = -uvw[2,i1]
         
         # add value for text file:
         line[ant_data2[i1]] = numpy.real( data[options.polarisation,options.channel,i1] )
         line_count += 1

         if curr_ant2 > (prev_ant2+1) :
            print("WARNING : missed antenna2 = %d" % (prev_ant2+1))   
      
         if curr_ant2 > prev_ant2 :    
            # next antenna2 with antenna1 
            baseline_count += 1
            all_baselines  += 1      
         else :      
            # antenna2 starts again :            
#            write_corr_matrix_line( vis_file_re_f, line, line_count )
            line=numpy.ones(n_ant)*numpy.nan
            line_count=0
            
            # print number of baselines for previous ANT1 (as it has clearly just changed as curr_ant2 <= prev_ant2 :
            baselines_expected = n_ant - prev_ant1
            missing_baselines = (baselines_expected-baseline_count)
            print("ANTENNA1 %d has %d baselines (exepected %d -> %d missing)" % (prev_ant1,baseline_count,baselines_expected,missing_baselines))
            baseline_count = 1

#            if curr_ant1 == 0 and curr_ant2!=curr_ant1 :
#               print("WARNING : antenna1 = %d missed antenna2 = 0" % (curr_ant1))
         
      else :
         baselines_expected = n_ant - prev_ant1
         missing_baselines = (baselines_expected-baseline_count)
         print("ANTENNA1 %d has %d baselines (exepected %d -> %d missing)" % (prev_ant1,baseline_count,baselines_expected,missing_baselines))


         print("DEBUG : curr_ant1 = %d and prev_ant1 = %d -> END ?" % (curr_ant1,prev_ant1))
         break
         

      print("DEBUG : i1=%d -> ant1 = %d , ant2 = %d" % (i1,ant_data1[i1],ant_data2[i1]))
      prev_ant1 = curr_ant1
      prev_ant2 = curr_ant2

   hdulist_re.writeto(out_vis_refits,clobber=True)
   hdulist_im.writeto(out_vis_imfits,clobber=True)
   hdulist_u.writeto(out_u_fits,clobber=True)
   hdulist_v.writeto(out_v_fits,clobber=True)
   hdulist_w.writeto(out_w_fits,clobber=True)

   vis_file_re_f.close()
