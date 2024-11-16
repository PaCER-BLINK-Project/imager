#!/usr/bin/env python
"""
Read a calibration solution binary file as produced by Andre Offringa's tools
These files are not documented outside the source code, but you can find the relevant code by grepping for WriteSolution.
The fit_complex_gains may be useful more widely so it is kept independent from tha aocal stuff.
It should never be necessary to import the AOClass itself, rather it can be returned from the fromfile and zeros functions.

SCRIPT originally developed in MWA_Tools by John Morgan et al 
with couple of tweaks to use it for ASVO calibration database by Marcin Sokolowski
"""
# import pdb

import sys, os, struct, logging, glob, errno
from collections import namedtuple
import numpy as np
import cmath
import math
import pyfits

from optparse import OptionParser,OptionGroup

# should be optional :
import matplotlib.pyplot as pyplot

# MWA global variables :
N_FINE_CHANNELS_PER_COARSE=32 # assuming processing always 24*32 = 762 fine channels 
N_FINE_CHANNELS=768
COARSE_CHANNEL_WIDTH=1.28
FINE_CHANNEL=0.04 # in MHz 

HEADER_FORMAT = "8s6I2d"
HEADER_SIZE = struct.calcsize(HEADER_FORMAT)
HEADER_INTRO = "MWAOCAL\0"

Header = namedtuple("header",
                    "intro fileType structureType intervalCount antennaCount channelCount polarizationCount timeStart timeEnd")
Header.__new__.__defaults__ = (HEADER_INTRO, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

assert struct.calcsize("I") == 4, "Error, unexpected unsigned int size used by python struct"
assert struct.calcsize("c") == 1, "Error, unexpected char size used by python struct"
assert struct.calcsize("d") == 8, "Error, unexpected double size used by python struct"
if not sys.byteorder == 'little':
    logging.warn("byteorder=%s", sys.byteorder)
else:
    logging.debug("byteorder=%s", sys.byteorder)


def mkdir_p(path):
   try:
      os.makedirs(path)
   except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
         pass
      else: raise

def list_tile_name( table, out_file, flagged_only=False ) :
   f = open(out_file, 'w')
   
   n_inputs = table.shape[0]
   tiles_out={}
   for i in range(0,n_inputs):
     idx=table[i][0]
     tile_idx=table[i][1]
     tile_id=table[i][2]
     tile_name=table[i][3]
     tile_pol=table[i][4]
     delays=table[i][12]
     flag=table[i][7]
     
     tile_out_name="T%03d%s" % (tile_idx,tile_pol)

     # print "DEBUG : %s %s %d" % (tile_name,tile_pol,tile_idx)

     if tile_pol == "X" :
        # print "flagged_only = %s , flag = %d" % (flagged_only,flag)
        if not flagged_only or flag > 0  :
           outline = str(tile_name) + " " + str(tile_id) + " " + str(flag) + " " + str(tile_idx) + "\n"
           f.write(outline)
           
           tiles_out[tile_idx] = tile_name
        #else :
        #   print "Skipped tile %s" % tile_name


   f.close()

   return tiles_out

def fit_complex_gains(v, mode='model', amp_order=5, fft_pad_factor=8):
    """
    Fit amplitude & phases of a 1D array of complex values (v).

    Returns 1D array of complex values corresponding to the model

    NaN is treated as a flag. These values are excluded from the fit and
    persist in the returned model.

    Weight of each phase is taken to be abs(v)**-2 as is appropriate for
    interferometer calibration solutions

    A Coarse solution for the phases is determined via a padded FFT and average
    offset. This is then refined with a two-parameter least squares

    Amplitudes are fit using a polynomial, unless amp_order is set to <1, in
    which case amplitudes are preserved.
    """
    good = ~np.isnan(v)  # mask matching non-NaN values
    if sum(good) == 0:
        return v
    v_fft = np.fft.fft(np.nan_to_num(v * np.abs(v) ** -3), n=fft_pad_factor * len(v))
    v_index = np.arange(len(v), dtype=np.float)
    gradient = float(
        np.abs(v_fft).argmax()) / fft_pad_factor / 768.  # change in phase per increment of v due to phase wrap
    wrap = np.array(np.cos(2 * np.pi * v_index * gradient), dtype=np.complex128)
    wrap.imag = np.sin(2 * np.pi * v_index * gradient)

    # unwrap v, keeping only valid values
    u = v[good] / wrap[good]
    u_index = np.arange(len(v), dtype=np.float)[good]
    # centre on 0
    u_mean = np.average(u)
    u_mean /= np.abs(u_mean)
    u0 = u / u_mean
    if np.var(np.angle(u0)) > 1:
        logging.warn("high variance detected in phases, check output model!")
    # finally do least squares
    m, c = np.polyfit(u_index, np.angle(u0), 1, w=np.abs(u) ** -2)
    fit_complex = np.array(np.cos(v_index * m + c), dtype=np.complex128)
    fit_complex.imag = np.sin(v_index * m + c)
    # fit poly to amplitudes
    if amp_order > 0:
        amp_model_coeffs = np.polyfit(v_index[good], np.abs(v[good]), amp_order)
        amp_model = np.poly1d(amp_model_coeffs)(v_index)
    else:
        amp_model = np.abs(v)
    if mode == "model":
        return np.where(good, amp_model * fit_complex * wrap * u_mean, np.nan)
    elif mode == "clip":
        # check np.angle(u0) statistics
        # set "good" array to False where statistics are bad
        # return original array where not newly flagged
        # return np.where(good, v, np.nan)
        raise RuntimeError, "clip not implemented"
    else:
        raise RuntimeError, "mode %s not implemented" % mode


class AOCal(np.ndarray):
    """
    AOCAl stored as a numpy array (with start and stop time stored as floats)

    Array is of dtype complex128 with the following dimensions:

    - calibration interval
    - antenna
    - channel
    - polarisation (order XX, XY, YX, YY)

    The following attributes are made available for convenience, however they
    are not stored explicitly, just read from the array shape.

    aocal.n_int
    aocal.n_ant
    aocal.n_chan
    aocal.n_pol
    """

    def __new__(cls, input_array, time_start=0.0, time_end=0.0, header_string=None):
        """
        See http://docs.scipy.org/doc/numpy-1.10.1/user/basics.subclassing.html
        """
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.time_start = float(time_start)
        obj.time_end = float(time_end)
        obj.header_string = header_string
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.time_start = getattr(obj, 'time_start', None)
        self.time_end = getattr(obj, 'time_end', None)

    def __getattr__(self, name):
        if name == 'n_int':
            return self.shape[0]
        elif name == 'n_ant':
            return self.shape[1]
        elif name == 'n_chan':
            return self.shape[2]
        elif name == 'n_pol':
            return self.shape[3]
        elif name == 'time_start':
            # required to avoid infinite recursion
            return object.__getattribute__(self, 'time_start')
        elif name == 'time_end':
            # required to avoid infinite recursion
            return object.__getattribute__(self, 'time_end')
        else:
            raise AttributeError, "AOCal has no Attribute %s. Dimensions can be accessed via n_int, n_ant, n_chan, n_pol" % name

    def strip_edge(self, n_chan):
        """
        return a copy of the array with edge channels removed

        useful for printing without nans but don't write out as calibration solution!
        """
        return self[:, :, n_chan:-n_chan, :]

    def tofile(self, cal_filename):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError, "array must have 4 dimensions and be of type complex128"
        header = Header(intervalCount=self.shape[0], antennaCount=self.shape[1], channelCount=self.shape[2],
                        polarizationCount=self.shape[3], timeStart=self.time_start, timeEnd=self.time_end)
        with open(cal_filename, "wb") as cal_file:
            header_string = struct.pack(HEADER_FORMAT, *header)
            cal_file.write(header_string)
            logging.debug("header written")
            cal_file.seek(HEADER_SIZE, os.SEEK_SET)  # skip header. os.SEEK_SET means seek relative to start of file
            np.ndarray.tofile(self, cal_file)
            logging.debug("binary file written")

    def tostring(self):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError, "array must have 4 dimensions and be of type complex128"
        header = Header(intervalCount=self.shape[0], antennaCount=self.shape[1], channelCount=self.shape[2],
                        polarizationCount=self.shape[3], timeStart=self.time_start, timeEnd=self.time_end)

        header_string = struct.pack(HEADER_FORMAT, *header)
        logging.debug("header written")

        aocal_string = np.ndarray.tobytes(self)
        logging.debug("binary file written")

        return (header_string + aocal_string)

    def header_tostring( self, intervalCount, antennaCount, channelCount, polarizationCount ):
        header = Header(intervalCount=intervalCount, antennaCount=antennaCount, channelCount=channelCount,
                        polarizationCount=polarizationCount, timeStart=self.time_start, timeEnd=self.time_end)

        header_string = struct.pack(HEADER_FORMAT, *header)
        return header_string
    

    def fit(self, pols=(0, 3), mode='model', amp_order=5):
        if not (np.iscomplexobj(self) and self.itemsize == 16 and len(self.shape) == 4):
            raise TypeError, "array must have 4 dimensions and be of type complex128"
        fit_array = np.zeros(self.shape, dtype=np.complex128)
        for interval in xrange(self.shape[0]):
            for antenna in xrange(self.shape[1]):
                logging.debug("fitting antenna %d" % antenna)
                for pol in pols:
                    v = self[interval, antenna, :, pol]
                    if sum(~np.isnan(self[interval, antenna, :, pol])) > 0:
                        self[interval, antenna, :, pol] = fit_complex_gains(self[interval, antenna, :, pol],amp_order=amp_order)

    # def toJSON(self):
    #     return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)


def ones(n_interval=1, n_antennas=128, n_channel=3072, n_pol=4, time_start=0.0, time_end=0.0):
    """
    produce an aocal with all complex gains set to amp 1, phase 0.
    """
    return AOCal(np.ones((n_interval, n_antennas, n_channel, n_pol), dtype=np.complex128), time_start, time_end)


def zeros(n_interval=1, n_antennas=128, n_channel=3072, n_pol=4, time_start=0.0, time_end=0.0):
    """
    produce an aocal with all complex gains set to amp 0, phase 0.
    """
    return AOCal(np.zeros((n_interval, n_antennas, n_channel, n_pol), dtype=np.complex128), time_start, time_end)


def fromfile(cal_filename,debug=0):
    """
    Read AOCal from file.
    """
    header_string = None
    with open(cal_filename, "rb") as cal_file:
        header_string = cal_file.read(struct.calcsize(HEADER_FORMAT))
        header = Header._make(struct.unpack(HEADER_FORMAT, header_string))
        logging.debug(header)
        if debug > 0 :
            print "header = %s / header_string = %s" % (header,header_string)
        assert header.intro == HEADER_INTRO, "File is not a calibrator file"
        assert header.fileType == 0, "fileType not recognised. Only 0 (complex Jones solutions) is recognised in mwatools/solutionfile.h as of 2013-08-30"
        assert header.structureType == 0, "structureType not recognised. Only 0 (ordered real/imag, polarization, channel, antenna, time) is recognised in mwatools/solutionfile.h as of 2013-08-30"
        logging.debug("header OK")

        count = header.intervalCount * header.antennaCount * header.channelCount * header.polarizationCount
        assert os.path.getsize(cal_filename) == HEADER_SIZE + 2 * count * struct.calcsize(
            "d"), "File is the wrong size."
        logging.debug("file correct size")
        cal_file.seek(HEADER_SIZE, os.SEEK_SET)  # skip header. os.SEEK_SET means seek relative to start of file
        
        if debug > 0 :
            print "count = %d" % (count)

        data = np.fromfile(cal_file, dtype=np.complex128, count=count)
    # print "ORIGINAL data.shape   = %d" % (data.shape[0])        
    shape = [header.intervalCount, header.antennaCount, header.channelCount, header.polarizationCount]
    data = data.reshape(shape)
    # print "RE-SHAPED data.shape  = %d" % (data.shape[0])
    new_aocal = AOCal(data, header.timeStart, header.timeEnd, header_string)
    
    return new_aocal


def merge_bin_files( bin_file_list, outfile ) :
    ao_list=[]
    
    intervalCount = 1
    antennaCount  = 128
    channelCount  = 0
    polarizationCount = 4
    for bin_file in bin_file_list :
        ao = fromfile( bin_file )
        ao_list.append( ao.copy() )
        
        if ao.n_int != intervalCount :
            intervalCount = ao.n_int
        if ao.n_pol != polarizationCount :
            polarizationCount = ao.n_pol
        if ao.n_ant != antennaCount :
            antennaCount = ao.n_ant;  

        # sum channels from picket fence :
        channelCount += ao.n_chan

    print "Read %d ao files" % (len(ao_list))
    
    out_ao_data = np.dstack( (ao_list) )
    
    # def header_tostring( intervalCount, antennaCount, channelCount, polarizationCount ):
    header_string = ao_list[0].header_tostring( intervalCount, antennaCount, channelCount, polarizationCount )
    out_ao = AOCal( out_ao_data, ao_list[0].time_start, ao_list[0].time_end, header_string )
    
    out_ao.tofile( outfile )
    
    return out_ao

def final_n_channels( binfile, final_channel_number, outfile ) :
   ao = fromfile( binfile )
   
   if ao.n_chan != final_channel_number :
      if final_channel_number < ao.n_chan :
         n_avg = ao.n_chan / final_channel_number
         print "Reducing number of channels from %d to %d requires averaging of %d channels" % (ao.n_chan,final_channel_number,int(n_avg))
         average_channels( binfile, int(n_avg), outfile )
      else :
         n_split = final_channel_number / ao.n_chan
         print "Increasing number of channels from %d to %d requires split by %d channels" % (ao.n_chan,final_channel_number,int(n_split))
         split_bin_file( binfile, int(n_split), outfile )
   else :
      print "File %s already has %d channels -> nothing to be done" % (binfile,ao.n_chan)
 
   return False

def average_channels( binfile, avg_n_channels, outfile ) :
   ao = fromfile( binfile )
   rest = ( ao.n_chan % avg_n_channels )

   if avg_n_channels >= 0 and avg_n_channels < ao.n_chan and rest==0 :
      ao_out = ao.copy()
      ao_out_shape = np.copy(ao.shape)
      ao_out_shape[2] = ao_out_shape[2] / avg_n_channels
      ao_out = AOCal( np.zeros( ao_out_shape, dtype=ao.dtype ) , ao.time_start, ao.time_end, ao.header_tostring )
      ao_out.n_chan = ao.n_chan / avg_n_channels

      
      for out_ch in range(0,ao_out.n_chan) :
         for offset in range(0,avg_n_channels) :
            ao_out[:,:,out_ch,:] += ao[:,:,out_ch*avg_n_channels + offset,:]
            
         ao_out[:,:,out_ch,:] = ao_out[:,:,out_ch,:] / avg_n_channels

      if outfile.find(".bin") < 0 :
         outfile += ".bin"
                  
      ao_out.tofile( outfile )      
      print "Saved .bin file with %d channels averaged to %s" % (avg_n_channels,outfile)
      
   else :
      print "ERROR : avg_n_channels = %d ( <0 or >= %d) or rest = %d != 0" % (avg_n_channels,ao.n_chan,rest)   
 
   return False
   
def split_bin_file( binfile, split_by_n_channels, outfilebase ):
   ao = fromfile( binfile )
   rest = ( ao.n_chan % split_by_n_channels )

   if split_by_n_channels >= 0 and split_by_n_channels < ao.n_chan and rest==0:
      n_chunks = ao.n_chan / split_by_n_channels
         
      for i in range(0,n_chunks) :
         outfile = ( "%s_%03d.bin" % (outfilebase,i))
         ao_out = ao.copy()
         ao_out = ao[:,:,i*split_by_n_channels:(i+1)*split_by_n_channels,:]
         ao_out.n_chan = split_by_n_channels
         
         ao_out.tofile( outfile )

         print "%d : channels %d - %d saved to file %s" % (i,i*split_by_n_channels,(i+1)*split_by_n_channels,outfile)
      
      return True
   else :
#      print "ERROR : split_by_n_channels = %d (either <0 or >= %d or rest = %d != 0)" % (split_by_n_channels,ao.n_chan,rest)
      print "WARNING : experimental version spliting into two unequal parts"
      
      n_chunks = 2

      i = 0 
      outfile = ( "%s_%03d.bin" % (outfilebase,i))
      ao_out = ao.copy()
      ao_out = ao[:,:,0:split_by_n_channels,:]
      ao_out.n_chan = split_by_n_channels         
      ao_out.tofile( outfile )
      print "%d : channels %d - %d saved to file %s" % (i,i*split_by_n_channels,(i+1)*split_by_n_channels,outfile)

      i += 1
      outfile = ( "%s_%03d.bin" % (outfilebase,i))
      ao_out = ao.copy()
      ao_out = ao[:,:,split_by_n_channels:ao.shape[2],:]
      ao_out.n_chan = (ao.shape[2] - split_by_n_channels)
      ao_out.tofile( outfile )
      print "%d : channels %d - %d saved to file %s" % (i,i*split_by_n_channels,ao.shape[2],outfile)
      
      return True

      
   return False   
    

# rts_filename_pattern has to contain "node" otherwise it will not work as it will not find coarse channel correctly 
def rtsfile(metafitsfile, rts_filename_pattern="DI_JonesMatrices_node[0-9]*.dat", aocal_format=True):
    import astropy.io.fits as fits

    """
    Read DI Jones matrices from RTS output files and convert to "aocal" format.
    Assumes RTS solutions are one per coarse channel.
    Needs the associated metafits file to get the antenna ordering right.
    """
    # (Relative) comparison of RTS and OFFRINGA polarisation ordering:
    # OFFRINGA:  XX-R  XX-I  XY-R  XY-I  YX-R  YX-I  YY-R  YY-I
    # RTS:       YY-R  YY-I  YX-R  YX-I  XY-R  XY-I  XX-R  XX-R
    pol_map = [3, 2, 1, 0]

    # Antenna reording:
    hdu = fits.open(metafitsfile)
    ant_map = hdu['TILEDATA'].data['Antenna'][::2]  # only want each tile once
    hdu.close()

    # Assumptions:
    nintervals = 1
    npols = 4

    # Get file names
    rts_filenames = sorted(
        glob.glob(rts_filename_pattern))  # <-- Assumes file names are appropriately ordered by channel
    rts_filenames.reverse()
    nchannels = len(rts_filenames)

    for chan in range(len(rts_filenames)):
        rts_filename = rts_filenames[chan]

        # get coarse channel from the file name chan is not good as the order might be reversed !!! as it happen in my tests 
        idx=rts_filename.find("node")
        coarse_channel = int( rts_filename[idx+4:idx+4+3] ) - 1 

        
        with open(rts_filename, "r") as rts_file:
            Jref = float(
                rts_file.readline())  # Common factor of all gains is a single number in the first line of the file
            rts_file.readline()  # The second line contains the model primary beam Jones matrix (in the direction of the calibrator)
            lines = rts_file.readlines()

            # If first time through, get number of antennas and set up data array for solution
            if chan == 0:
                nantennas = len(lines)
                # Create numpy array structure
                data = np.empty((nintervals, nantennas, nchannels, npols,), dtype=np.complex128)
                data_out = np.empty((nintervals, nantennas, 768, npols,), dtype=np.complex128)
                data[:] = np.nan
                data_out[:] = np.nan
            else:
                assert len(lines) == nantennas, "Files contain different numbers of antennas"

            # Parse each line
            for ant in range(len(lines)):
                line = lines[ant]
                jones_str = line.split(",")
                assert len(jones_str) == 2 * npols, "Incorrect number of elements in Jones matrix"
                ant_idx = ant_map[ant]
                for pol in range(len(pol_map)):
                    p = pol_map[pol]
                    data[0, ant_idx, coarse_channel, pol] = float(jones_str[2 * p]) + float(jones_str[2 * p + 1]) * 1j

                    for finech in range(0, 32):
                        data_out[0, ant_idx, finech + coarse_channel * 32, pol] = data[0, ant_idx, coarse_channel, pol]

                        if aocal_format:
                            data_out[0, ant_idx, finech + coarse_channel * 32, pol] = 1.00 / data_out[ 0, ant_idx, finech + coarse_channel * 32, pol]

                    if ant == 2 and pol == 0:
                        print "%d %.4f" % (coarse_channel, abs(data_out[0, ant_idx, coarse_channel, pol]))

    new_aocal = AOCal(data_out, 0, 0)
    return new_aocal


# rts_filename_pattern has to contain "node" otherwise it will not work as it will not find coarse channel correctly 
def rtsfile_bandpass(metafitsfile, rts_filename_pattern="BandpassCalibration_node[0-9]*.dat", aocal_format=True,
                     use_fit=False,debug=-1):
    # def rtsfile_bandpass(metafitsfile, rts_filename_pattern="BandpassCalibration_node001.dat", aocal_format=True,use_fit=False):
    import astropy.io.fits as fits

    """
    Read DI Jones matrices from RTS output files and convert to "aocal" format.
    Assumes RTS solutions are one per coarse channel.
    Needs the associated metafits file to get the antenna ordering right.
    """
    # (Relative) comparison of RTS and OFFRINGA polarisation ordering:
    # OFFRINGA:  XX-R  XX-I  XY-R  XY-I  YX-R  YX-I  YY-R  YY-I
    # RTS:       YY-R  YY-I  YX-R  YX-I  XY-R  XY-I  XX-R  XX-R
    pol_map = [3, 2, 1, 0]

    # Antenna reording:
    hdu = fits.open(metafitsfile)
    ant_map = hdu['TILEDATA'].data['Antenna'][::2]  # only want each tile once
    hdu.close()

    # Assumptions:
    nintervals = 1
    npols = 4
    offset = 0
    if use_fit:
        offset = 1

        # Get file names
    rts_filenames = sorted(
        glob.glob(rts_filename_pattern))  # <-- Assumes file names are appropriately ordered by channel
    rts_filenames.reverse()
    nchannels = len(rts_filenames)

    for chan in range(len(rts_filenames)):
        rts_filename = rts_filenames[chan]
        
        # get coarse channel from the file name chan is not good as the order might be reversed !!! as it happen in my tests 
        idx=rts_filename.find("node")
        coarse_channel = int( rts_filename[idx+4:idx+4+3] ) - 1 
        
        with open(rts_filename, "r") as rts_file:
            freq_list = rts_file.readline()  # list of frequencies
            freq_list_str = freq_list.split(",")
            freq_list_float = np.fromstring(freq_list, sep=",")
            lines = rts_file.readlines()

            # If first time through, get number of antennas and set up data array for solution
            if chan == 0:
                # nantennas = len(lines)/8 
                nantennas = int(lines[len(lines) - 1].split(",")[0])
                # print "File %s : number of antennas = %d" % (rts_filename,nantennas)
                # Create numpy array structure
                data = np.empty((nintervals, nantennas, nchannels, npols,), dtype=np.complex128)
                data_out = np.empty((nintervals, nantennas, 768, npols,), dtype=np.complex128)
                data[:] = np.nan
                data_out[:] = np.nan
            else:
                nantennas_new = int(lines[len(lines) - 1].split(",")[0])
                assert nantennas_new == nantennas, "Files contain different numbers of antennas"

            print "File %s : nantennas = %d , lines = %d" % (rts_filename, nantennas, len(lines))
            # print lines[0]
            # print lines[1]
            # print lines[2]
            # print lines[3]
            # print "----------------------------------------"

            # Parse each line
            for ant in range(1, nantennas + 1):
                # ant_line_idx = ant*8 # 8 lines per antenna
                ant_line_idx = -1
                for k in range(0, len(lines)):
                    aaa = int(lines[k].split(",")[0])
                    if aaa == ant:
                        ant_line_idx = k
                        break

                if ant_line_idx < 0:
                    print "WARNING : antenna = %d not found -> skipped" % (ant)
                    continue

                print "antenna = %d" % (ant)
                line_j0 = lines[ant_line_idx + offset]
                line_j1 = lines[ant_line_idx + 2 + offset]
                line_j2 = lines[ant_line_idx + 4 + offset]
                line_j3 = lines[ant_line_idx + 6 + offset]
                
                if debug > 0 : 
                    print "DEBUG : lines %d , %d , %d , %d" % (ant_line_idx+offset,ant_line_idx+2+offset,ant_line_idx+4+offset,ant_line_idx+6+offset)
                    print lines[ant_line_idx+offset]
                    print lines[ant_line_idx+2+offset]
                    print lines[ant_line_idx+4+offset]
                    print lines[ant_line_idx+6+offset]
                    print "----------------------------------------------------------------------"

                j0 = np.fromstring(line_j0, sep=",")
                j1 = np.fromstring(line_j1, sep=",")
                j2 = np.fromstring(line_j2, sep=",")
                j3 = np.fromstring(line_j3, sep=",")

                # print "DEBUG : idx = %d , antenna=%d starts at line = %d [%d,%d,%d,%d]" % (ant_line_idx+offset,ant,ant_line_idx,j0[0],j1[0],j2[0],j3[0])
                # print "%s" % (line_j0)
                # print "%s" % (line_j1)
                # print "%s" % (line_j2)
                # print "%s" % (line_j3)

                ant_idx = ant_map[ant - 1]
                # for pol in range(len(pol_map)):
                #    p = pol_map[pol]
                # data[0,ant_idx,coarse_channel,pol] = float(jones_str[2*p]) + float(jones_str[2*p+1])*1j

                assert j0[0] == j1[0] and j1[0] == j2[0] and j2[0] == j3[
                    0], "Wrong antenna number %d vs. %d vs. %d vs. %d" % (j0[0], j1[0], j2[0], j3[0])

                # print "size(freq_list) = %d vs. len(j0) = %d" % (len(freq_list_float),len(j0))
                # 1, +1.000000,+0.000000, +4.353953,+0.600645, +3.444445,-1.117209, +0.863041,+0.023113, +0.881922,-0.081904, +0.964471,-0.189353, +0.912126,-0.240132, +0.878874,-0.009400,
                ant_idx_file = j0[0]
                for i in range(0, len(freq_list_float)):
                    finech = int(freq_list_float[i] / 0.04)
                    
                    if debug > 0 :
                       print "ANTIDX_%d : %.6f %.6f" % (ant_idx,j0[1 + 2 * i],j0[1 + 2 * i + 1])
                       
                    data_out[0, ant_idx, finech + coarse_channel * 32, 0] = j0[1 + 2 * i] * cmath.exp(j0[1 + 2 * i + 1] * 1j)
                    data_out[0, ant_idx, finech + coarse_channel * 32, 1] = j1[1 + 2 * i] * cmath.exp(j1[1 + 2 * i + 1] * 1j)
                    data_out[0, ant_idx, finech + coarse_channel * 32, 2] = j2[1 + 2 * i] * cmath.exp(j2[1 + 2 * i + 1] * 1j)
                    data_out[0, ant_idx, finech + coarse_channel * 32, 3] = j3[1 + 2 * i] * cmath.exp(j3[1 + 2 * i + 1] * 1j)

                    # for finech in range(0,32) :
                    #   data_out[0,ant_idx,finech+coarse_channel*32,pol] = data[0,ant_idx,coarse_channel,pol]

                    if aocal_format :
                        data_out[0,ant_idx,finech+coarse_channel*32,pol] = 1.00 / data_out[0,ant_idx,finech+coarse_channel*32,pol]

                    # if ant == 2 and pol==0 :
                    #   print "%d %.4f" % (coarse_channel,abs(data_out[0,ant_idx,coarse_channel,pol]))

    new_aocal = AOCal(data_out, 0, 0)
    return new_aocal


def dump_ao_calsolutions( calfile, options, ant=-1, do_phase=0, do_fit=0, do_reim=0, do_ampphase=0, out_basename_param="calsolutions", channels_str=None ) :
   out_channels=options.out_channels   
   band      = COARSE_CHANNEL_WIDTH
   half_band = COARSE_CHANNEL_WIDTH / 2.00
   ant_param = ant
   out_basename = out_basename_param
   
   start_cc=-1
   start_freq=0
   if options.channels2freq :
      if channels_str is None :
          if options.metafits is None :
              print "ERROR : metafits file name not provided, use option --metafits or -m -> cannot continue with option --ch2freq"
              sys.exit(-1)
          else :
              if os.path.exists(options.metafits) :
                  fits = pyfits.open(options.metafits)
                  header = fits[0].header   
                  try : 
                      channels_str=header['CHANNELS']
                  except :
                      print "ERROR : keyword CHANNELS not found in metafits file %s" % (options.metafits)                            
                  fits.close()
              else :
                  print "ERROR : metafits file %s does not exist -> cannot continue with option --ch2freq" % (options.metafits)    
                  sys.exit(-1)

      channels_list=channels_str.split(",")
      start_cc = int(channels_list[0])
      print "CHANNELS = %s -> start_cc = %d" % (channels_str,start_cc)

      # center of the first coarse channel is start_cc*band - half_band + FINE_CHANNEL/2
      start_freq = (start_cc*band) - half_band + FINE_CHANNEL/2
      print "start_cc = %d -> start_freq = %.4f MHz" % (start_cc,start_freq)
          
   
   caldata = fromfile( calfile )
   if do_fit > 0 :
      print "Fitting polynomial to amplitudes of calibration solutions"
      caldata.fit()
      calfile_fitted=calfile.replace('.bin', '_fit.bin')
      caldata.tofile( calfile_fitted )
      print "Saved fitted calibration solutions to file %s" % (calfile_fitted)

   # Header = namedtuple("header", "intro fileType structureType intervalCount antennaCount channelCount polarizationCount timeStart timeEnd")   
   print "Read calibration file calfile = %s" % (calfile)
   print "\tNumber of time steps    = %d" % (caldata.n_int)
   print "\tNumber of antenna       = %d" % (caldata.n_ant)
   print "\tNumber of channels      = %d" % (caldata.n_chan)
   print "\tNumber of polarisations = %d" % (caldata.n_pol)
   print "\tHeader string           = ||||%s||||" % (caldata.header_string)

   # if just one antenna specified ant>0 -> loop over just 1 ant :
   ant_low_range = ant
   ant_up_range  = ant+1 
   
   if ant < 0 :
      # dump all antennas :
      ant_low_range = 0
      ant_up_range  = caldata.n_ant
      print("DEBUG : dump all antennas in range %d - %d" % (ant_low_range,ant_up_range))


   filename_xx = out_basename + "_xx.txt"
   filename_xy = out_basename + "_xy.txt"
   filename_yx = out_basename + "_yx.txt"
   filename_yy = out_basename + "_yy.txt"

   file_xx = open( filename_xx , "w")
   file_xy = open( filename_xy , "w")
   file_yx = open( filename_yx , "w")
   file_yy = open( filename_yy , "w")


   for ant in range(ant_low_range,ant_up_range) :
       if ant >= 0 and ant<128 :
          # print "DEBUG : ant=%d , ant_param=%d, (ant_up_range-ant_low_range) = %d, find = %d" % (ant,ant_param,(ant_up_range-ant_up_range),out_basename_param.find("%"))
          if ant_param < 0 and (ant_up_range-ant_low_range)>=2 and out_basename_param.find("%")>=0 :
              out_basename = out_basename_param % (ant)
              print "Dumping antenna = %d from range %d - %d -> out_basename = %s" % (ant,ant_low_range,ant_up_range,out_basename)
              
       

          if caldata.n_chan != out_channels :
              caldata_new = zeros(n_interval=caldata.n_int, n_antennas=caldata.n_ant, n_channel=out_channels, n_pol=caldata.n_pol)
              avg_channels = 1
              if caldata.n_chan > out_channels :
                  avg_channels = ( caldata.n_chan / out_channels )

              print "WARNING : averaging of different number of channels (%d) than required %d -> averaging every %d channels" % (caldata.n_chan,out_channels,avg_channels)

              for integr_idx in range(0,caldata.n_int) :
                  for ant_idx in range(0,caldata.n_ant) :
                      for pol_idx in range(0,caldata.n_pol) :
                          for ch in range(0,out_channels) :

                              sum=0
                              cnt=0
                              for ch_fine in range(0,avg_channels):
                                  ch_idx = avg_channels*ch+ch_fine
                                  sum += caldata[integr_idx,ant_idx,ch_idx,pol_idx]
                                  cnt += 1

                              caldata_new[integr_idx,ant_idx,ch,pol_idx] = sum / cnt


              caldata = caldata_new.copy()
              caldata.n_chan = caldata_new.n_chan


          for ch in range(0,caldata.n_chan) :
             is_nan=False
             if options.swap_nans is not None :
                 if cmath.isnan( caldata[0,ant,ch,0] ) :
                     # caldata[0,ant,ch,0] = -10000
                     is_nan=True

             ch_freq = ch
             if options.channels2freq :
#                 ch_freq = start_freq + ch*FINE_CHANNEL
                  cc = ( ch  / N_FINE_CHANNELS_PER_COARSE ) 
                  ch_in_cc = ( ch % N_FINE_CHANNELS_PER_COARSE )
                  cc_absolute = int(channels_list[cc])
                  ch_freq = ( cc_absolute )*band -  half_band + FINE_CHANNEL/2 + ch_in_cc*FINE_CHANNEL
                  print "DEBUG : ch=%d -> cc = %d , ch_in_cc = %d -> (cc_absolute) = %d -> ch_freq = %.2f" % (ch,cc,ch_in_cc,cc_absolute,ch_freq)

             if is_nan and options.skip_nans :
                 print "INFO : nan value skipped in channel = %.4f (channel = %d)" % (ch_freq,ch)
                 continue

             # if gain_mag > 900 :
             # set nans :
             # caldata[0,ant,ch,0].real = math.nan
             # ...
             # ...
             # etc ...

             if do_ampphase > 0 :
                out_line_xx = "%.4f %.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,0]),cmath.phase(caldata[0,ant,ch,0])*(180.0/math.pi),ant)
                out_line_xy = "%.4f %.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,1]),cmath.phase(caldata[0,ant,ch,1])*(180.0/math.pi),ant)
                out_line_yx = "%.4f %.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,2]),cmath.phase(caldata[0,ant,ch,2])*(180.0/math.pi),ant)
                out_line_yy = "%.4f %.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,3]),cmath.phase(caldata[0,ant,ch,3])*(180.0/math.pi),ant)
             elif do_reim > 0 :
                out_line_xx = "%.4f %.4f %.4f %d\n" % (ch_freq,caldata[0,ant,ch,0].real,caldata[0,ant,ch,0].imag,ant)
                out_line_xy = "%.4f %.4f %.4f %d\n" % (ch_freq,caldata[0,ant,ch,1].real,caldata[0,ant,ch,1].imag,ant)
                out_line_yx = "%.4f %.4f %.4f %d\n" % (ch_freq,caldata[0,ant,ch,2].real,caldata[0,ant,ch,2].imag,ant)
                out_line_yy = "%.4f %.4f %.4f %d\n" % (ch_freq,caldata[0,ant,ch,3].real,caldata[0,ant,ch,3].imag,ant)
             else :
                # skip values = 1000 which are flagged channels due to RFI :
    #            if gain_mag < 900 :
                if do_phase > 0 :
                   if is_nan :
                       # if NaN -> change to -10000.00000000 so that it works the same way as my CASA based script :
                       out_line_xx = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_xy = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_yx = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_yy = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                   else :
                       out_line_xx = "%.4f %.4f %d\n" % (ch_freq,cmath.phase(caldata[0,ant,ch,0])*(180.0/math.pi),ch)
                       out_line_xy = "%.4f %.4f %d\n" % (ch_freq,cmath.phase(caldata[0,ant,ch,1])*(180.0/math.pi),ch)
                       out_line_yx = "%.4f %.4f %d\n" % (ch_freq,cmath.phase(caldata[0,ant,ch,2])*(180.0/math.pi),ch)
                       out_line_yy = "%.4f %.4f %d\n" % (ch_freq,cmath.phase(caldata[0,ant,ch,3])*(180.0/math.pi),ch)

                else :
                   if is_nan :
                       out_line_xx = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_xy = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_yx = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                       out_line_yy = "%.4f %.4f %d\n" % (ch_freq,-10000.00000000,ch)
                   else :
                       out_line_xx = "%.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,0]),ch)
                       out_line_xy = "%.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,1]),ch)
                       out_line_yx = "%.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,2]),ch)
                       out_line_yy = "%.4f %.4f %d\n" % (ch_freq,abs(caldata[0,ant,ch,3]),ch)

             file_xx.write( out_line_xx )
             file_xy.write( out_line_xy )
             file_yx.write( out_line_yx )
             file_yy.write( out_line_yy )

   file_xx.close()
   file_xy.close()
   file_yx.close()
   file_yy.close()

# 20181020 - changed to make it safe and use MEDIAN, otherwise a single outlier may completely spoil the mean :
#            alternatively remove all >2 ?
def mean_not_nan( array ) :
#   sum  = 0
#   sum2 = 0
#   cnt  = 0
#   
#   for i in range(0,array.size) :
#       if not np.isnan(array[i]) and abs(array[i])<1000000.00 :
#           sum  += array[i]
#           sum2 += (array[i]*array[i])
#           cnt  += 1
#      
#   mean = np.nan
#   rms  = np.nan        
#   if cnt > 0 :        
#      mean = sum / cnt
#      rms = math.sqrt( sum2/cnt - (mean*mean) )

   array_sorted = array.copy()
   array_sorted.sort()
   mean = array_sorted[array_sorted.size/2]
   
   # use interquarile range :
   q1 = int( array_sorted.size*0.25 ) 
   q3 = int( array_sorted.size*0.75 )
   rms  =  (array_sorted[q3] - array_sorted[q1]) / 1.349
   cnt = array_sorted.size  
   
   return (cnt,mean,rms)

                     
# max_rms - maximum allowed value of RMS                      
def calc_mean_rms( calfile, do_fit=0, verb=0, do_phase=False, max_rms=1.00, outfile=None ) :
   caldata_original = fromfile( calfile )     # to store original calibration solutions
   caldata          = caldata_original.copy() # to store fit
   residuals        = caldata_original.copy() # data to test quality 
   if outfile is None :
      outfile = calfile.replace(".bin","_stat.txt")

   if do_fit > 0 :
      print "Fitting polynomial to amplitudes of calibration solutions"
      caldata.fit( amp_order=do_fit )
      residuals = caldata_original - caldata

   rms_x = np.zeros( caldata_original.n_ant )
   rms_y = np.zeros( caldata_original.n_ant )
   mean_x = np.zeros( caldata_original.n_ant, dtype=complex )
   mean_y = np.zeros( caldata_original.n_ant, dtype=complex )
   out_cnt_x = np.zeros( caldata_original.n_ant )
   out_cnt_y = np.zeros( caldata_original.n_ant )
   out_nan_cnt_x = np.zeros( caldata_original.n_ant )
   out_nan_cnt_y = np.zeros( caldata_original.n_ant )
  
   ok_cnt_x = 0
   ok_cnt_y = 0
   
   for ant in range(0,residuals.n_ant) :
      # rms_x[ant] = residuals[0,ant,:,0].std()
      # rms_y[ant] = residuals[0,ant,:,3].std()
      # mean_x[ant] = residuals[0,ant,:,0].mean()
      # mean_y[ant] = residuals[0,ant,:,3].mean()
      
      print
      print "ANTENNA = %d" % (ant)
      
      # just a test :
      sum_x=0
      sum2_x=0
      sum_y=0
      sum2_y=0
      cnt_x=0
      cnt_y=0
      nan_cnt_x=0
      nan_cnt_y=0
      for ch in range(0,residuals.n_chan) :
         val_x = abs(residuals[0,ant,ch,0])
         val_y = abs(residuals[0,ant,ch,3])
         
         if do_phase :
             val_x = cmath.phase(residuals[0,ant,ch,0]) * (180.00/math.pi)
             val_y = cmath.phase(residuals[0,ant,ch,3]) * (180.00/math.pi)

        
         if verb > 0 : 
             print "\tDEBUG %d : %.8f %.8f" % (ch,val_x,val_y)

         if cmath.isnan(val_x) or abs(val_x)>1000000.00 :
             nan_cnt_x += 1
         else :
             sum_x  += val_x
             sum2_x += (val_x*val_x)
             cnt_x  += 1 

         if cmath.isnan(val_y) or abs(val_y)>1000000.00 :
             nan_cnt_y += 1
         else :
             sum_y  += val_y
             sum2_y += (val_y*val_y)
             cnt_y  += 1 

      
      test_mean_x = np.nan
      test_rms_x = np.nan      
      if cnt_x > 0 :
          test_mean_x = sum_x / cnt_x
          
          test_rms_x = np.nan
          sqrt_arg_x = sum2_x / cnt_x - test_mean_x*test_mean_x
          if sqrt_arg_x > 0 :
              test_rms_x = math.sqrt( sqrt_arg_x )
          else :
              print "\tERROR (X): negative value in SQRT( %.4f ) = SQRT( %.4f - %.4f)" % (sqrt_arg_x,(sum2_x / cnt_x),(test_mean_x*test_mean_x))
      else :
          print "\tERROR (X): no not-NaN values"

      
      test_mean_y = np.nan
      test_rms_y = np.nan      
      if cnt_y > 0 :
          test_mean_y = sum_y / cnt_y      
          
          sqrt_arg_y = sum2_y / cnt_y - test_mean_y*test_mean_y
          if sqrt_arg_y > 0 :      
              test_rms_y = math.sqrt( sum2_y / cnt_y - test_mean_y*test_mean_y )
          else :
              print "\tERROR (Y): negative value in SQRT( %.4f ) = SQRT( %.4f - %.4f)" % (sqrt_arg_y,(sum2_y / cnt_y),(test_mean_y*test_mean_y))
      else :
           print "\tERROR (Y): no not-NaN values"
      
      rms_x[ant] = test_rms_x
      rms_y[ant] = test_rms_y
      mean_x[ant] = test_mean_x
      mean_y[ant] = test_mean_y
      out_cnt_x[ant] = cnt_x
      out_cnt_y[ant] = cnt_y
      out_nan_cnt_x[ant] = nan_cnt_x
      out_nan_cnt_y[ant] = nan_cnt_y
      
      if not np.isnan(rms_x[ant]) and rms_x[ant] < max_rms :
         ok_cnt_x += 1

      if not np.isnan(rms_y[ant]) and rms_y[ant] < max_rms :
         ok_cnt_y += 1
       
      
      # print "ANT %03d : RMS X / Y = %.3f / %.3f ( test %.3f / %.3f ) , MEAN X / Y = %.3f / %.3f ( test %.3f / %.3f ) , non-NaN = %d / %d , NaN = %d / %d" % (ant,rms_x[ant],rms_y[ant],test_rms_x,test_rms_y,mean_x[ant],mean_y[ant],test_mean_x,test_mean_y,cnt_x,cnt_y,nan_cnt_x,nan_cnt_y)
      print "\tANT %03d : MEAN +/- RMS (X) = %.3f +/- %.3f , MEAN +/- RMS (Y) = %.3f +/- %.3f , non-NaN = %d / %d , NaN = %d / %d" % (ant,mean_x[ant],rms_x[ant],mean_y[ant],rms_y[ant],cnt_x,cnt_y,nan_cnt_x,nan_cnt_y)


   (cnt_ok_test_x,mean_mean_x,rms_mean_x) = mean_not_nan( mean_x )
   (cnt_ok_test2_x,mean_rms_x,rms_rms_x) = mean_not_nan( rms_x )           

   (cnt_ok_test_y,mean_mean_y,rms_mean_y) = mean_not_nan( mean_y )
   (cnt_ok_test2_y,mean_rms_y,rms_rms_y) = mean_not_nan( rms_y )           

   out_f = open( outfile , "w" )
   print
   line =  "# CALIBRATION SOLUTIONS STATISTICS :"
   print line
   out_f.write( line + "\n" )
   print
   line = "# X polarisation :"
   print line
   out_f.write( line + "\n" )
   
   line =  "# \t# ok tiles = %d ( %d / %d )" % (ok_cnt_x,cnt_ok_test_x,cnt_ok_test2_x)
   print line
   out_f.write( line + "\n" )

   line = "# \t<MEAN>     = %.4f +/- %.4f" %  (mean_mean_x,rms_mean_x)
   print line
   out_f.write( line + "\n" )

   line =  "# \t<RMS>      = %.4f +/- %.4f" %  (mean_rms_x,rms_rms_x)
   print line
   out_f.write( line + "\n\n" )

   print
   line = "# Y polarisation :"
   print line
   out_f.write( line + "\n" )

   line =  "# \t# ok tiles = %d ( %d / %d )" % (ok_cnt_y,cnt_ok_test_y,cnt_ok_test2_y)
   print line
   out_f.write( line + "\n" )

   line =  "# \t<MEAN>     = %.4f +/- %.4f" %  (mean_mean_y,rms_mean_y)
   print line
   out_f.write( line + "\n" )

   line =  "# \t<RMS>      = %.4f +/- %.4f" %  (mean_rms_y,rms_rms_y)
   print line
   out_f.write( line + "\n" )

   print
   
   out_f.write( "# ANT  CNT_OK_X   CNT_OK_Y   MEAN_X    MEAN_Y    RMS_X    RMS_X    NAN_CNT_X    NAN_CNT_Y\n")
   for ant in range(0,residuals.n_ant) : 
      line =    "%03d      %d        %d      %.4f    %.4f    %.4f   %.4f      %d           %d\n" % (
                ant,out_cnt_x[ant],out_cnt_y[ant],mean_x[ant],mean_y[ant],rms_x[ant],rms_y[ant],out_nan_cnt_x[ant],out_nan_cnt_y[ant])
      out_f.write( line )                                                                            
   
   out_f.close()

   return (ok_cnt_x,mean_x,rms_x,out_cnt_x,out_nan_cnt_x,mean_mean_x,rms_mean_x,mean_rms_x,rms_rms_x,
           ok_cnt_y,mean_y,rms_y,out_cnt_y,out_nan_cnt_y,mean_mean_y,rms_mean_y,mean_rms_y,rms_rms_y)
              
   
def plotcal( calfile, caldata=None, nx=16, ny=8, outdir="images/", do_show=True, min_y=0, max_y=1, do_fit=0, metafits=None, plotall=True, phase=0, obsid=-1, block_image=False, wrong_channels=0 ) :
   basename = calfile.replace(".bin","")
   if obsid <= 0 : 
       obsid = int(basename)
   outdir = "%s/%s/" % (outdir,basename)   
   mkdir_p( outdir )
   
   if phase > 0 :
      min_y = -190.0
      max_y = +190.0
   
   if metafits is None :
       metafits = "%s.metafits" % (basename)
      
   if not os.path.exists(metafits):
      # download if metafits does not exist :
#      wget_string = "wget http://mwa-metadata01.pawsey.org.au/metadata/fits/?obs_id=%d -O %d.metafits" % (obsid,obsid)
      wget_string = "wget http://ws.mwatelescope.org/metadata/fits?obs_id=%d -O %d.metafits" % (obsid,obsid)
      print "%s" % wget_string
      os.system(wget_string)
      
      
   tiles=None    
   if os.path.exists(metafits):
      out_tile_list="%s.tile_list" % (basename)
   
      fits = pyfits.open(metafits)
      table = fits[1].data
      tiles = list_tile_name( table, out_tile_list )
      fits.close()


   if caldata is None :
      caldata = fromfile( calfile )

   calfit = None 
   if do_fit > 0 :       
      calfit = caldata.copy()   
      print "Fitting polynomial to amplitudes of calibration solutions"
      calfit.fit( amp_order=do_fit )
      
   fig  = pyplot.figure(0,figsize=(20,10))
#   pyplot.title( calfile )
   
   channels = np.arange(0,caldata.n_chan)

   if plotall :   
       drawn=0
       page=1
       images_per_page = nx*ny
       image_on_page   = 0
   
       print "Drawing %d x %d = %d antennas per page (%d ants in %s)" % (nx,ny,images_per_page,caldata.n_ant,calfile)
   
       for ant in range(0,caldata.n_ant) :
           ant_str = "%03d" % (ant)
           if tiles is not None :
               ant_str = "%03d / %s" % (ant,tiles[ant])
           # if image_on_page == 0 :
           #     title = calfile + (" page %02d" % page )
           #     pyplot.text( caldata.n_chan, max_y*1.2, title )
           #    pyplot.title( calfile + (" page %02d" % page ) )
       
           ax = fig.add_subplot( nx, ny , (image_on_page+1) ) # create axes within the figure : we have 1 plot in X direction, 1 plot in Y direction and we select plot 1
#           ax = fig.add_subplot( 1, 1 , ant+1)
       #    ax.set_xlim(0,caldata.n_chan)
           ax.set_ylim(bottom=min_y,top=max_y)


           if ( image_on_page % nx ) == 0 :
              # ax.set_xlabel('Channel')
              if phase > 0 :
                  ax.set_ylabel('Phase [deg]')
              else :
                  ax.set_ylabel('Amplitude')
           if ( image_on_page / nx ) == (ny-1) :
              ax.set_xlabel('Channel')
       
           val_x = abs(caldata[0,ant,:,0])
           val_y = abs(caldata[0,ant,:,3])
           
           if phase > 0 :
               val_x =  np.angle(caldata[0,ant,:,0],deg=True) 
               val_y =  np.angle(caldata[0,ant,:,3],deg=True) 

       
#           print "DEBUG : ant = %d -> size = %d / %d -> min/max = %.4f / %.4f" % (ant,val_x.size,val_y.size,np.nanmin(val_x),np.nanmax(val_x))
           if val_x.size <= 0 or ( np.isnan( np.min(val_x) ) and np.isnan( np.max(val_x) ) ):
               ax.set_xlim(0,caldata.n_chan)
       

           ax.plot( channels, val_x , color="blue" , marker="x", markersize=1 ) # green=Y
           ax.plot( channels, val_y , color="green", marker="x", markersize=1 ) # green=Y
#           ax.scatter( channels, val_x , color="blue", marker="+", markersize=1 )  # blue =X
#           ax.scatter( channels, val_y , color="green", marker="+", markersize=1 ) # green=Y
           pyplot.text( caldata.n_chan*0.4 , max_y*0.8, ant_str )
       
           if image_on_page == 0 :
               title = calfile + (" page %02d" % page )
               pyplot.text( caldata.n_chan, max_y*1.2, title )
          #    pyplot.title( calfile + (" page %02d" % page ) )

       
           if calfit is not None :
               fit_x = abs(calfit[0,ant,:,0])
               fit_y = abs(calfit[0,ant,:,3])
               
               if phase > 0 :
                   fit_x = np.angle(calfit[0,ant,:,0],deg=True)
                   fit_y = np.angle(calfit[0,ant,:,3],deg=True)
          
               ax.plot( channels, fit_x , color="black", linewidth=1 )
               ax.plot( channels, fit_y , color="yellow", linewidth=1 )
       
           drawn += 1
           image_on_page += 1
       
           if image_on_page == images_per_page and image_on_page < caldata.n_ant:
               # save page 
               pngfile = "%s/%s_amp_%.02d.png" % (outdir,basename,page)
               if phase > 0 :
                   pngfile = "%s/%s_phase_%.02d.png" % (outdir,basename,page)

               print "Saved plot of %d antennas to file %s" % (image_on_page,pngfile)

               pyplot.savefig(pngfile)               
               page += 1 
               image_on_page = 0
               if do_show :
                   pyplot.show(block=False)

               pyplot.pause(1)
               pyplot.clf()

           # pyplot.title( calfile + (" page %02d" % page ) )
           
#       if ant > 10 :       
#          break

       # plot with all tiles in 1 image :
       if (nx*ny) >= 128 :
           pngfile = "%s/%s_amp.png" % (outdir,basename)
           if phase > 0 :
              pngfile = "%s/%s_phase.png" % (outdir,basename)

           pyplot.savefig(pngfile)
           print "Saved file %s" % (pngfile)
   
           if do_show :
               pyplot.show(block=block_image)

#  calc_statistics and show :
   (ok_cnt_x,mean_x,rms_x,out_cnt_x,out_nan_cnt_x,mean_mean_x,rms_mean_x,mean_rms_x,rms_rms_x,
    ok_cnt_y,mean_y,rms_y,out_cnt_y,out_nan_cnt_y,mean_mean_y,rms_mean_y,mean_rms_y,rms_rms_y) = calc_mean_rms( calfile, do_fit=options.do_fit, do_phase=0 )

   pyplot.clf()  
   
   mean_mean_x_str = "%.3f" % (mean_mean_x)
   if len(mean_mean_x_str) > 20 :
       mean_mean_x_str = "NaN (too long)"
   mean_rms_x_str  = "%.3f" % (mean_rms_x)
   if len(mean_rms_x_str) > 20 :
       mean_rms_x_str = "NaN (too long)"

   mean_mean_y_str = "%.3f" % (mean_mean_y)
   if len(mean_mean_y_str) > 20 :
       mean_mean_y_str = "NaN (too long)"
   mean_rms_y_str  = "%.3f" % (mean_rms_y)
   if len(mean_rms_y_str) > 20 :
       mean_rms_y_str = "NaN (too long)"
   
   title_color = 'black'
   quality_str = "OK"   
   if ok_cnt_x < 0.75*caldata.n_ant or ok_cnt_y < 0.75*caldata.n_ant :
      title_color = 'red'
      quality_str = "BAD"
      
   title = "%s , %s ANT = X : %d / %s / %s , Y : %d / %s / %s" % (basename, quality_str, ok_cnt_x, mean_mean_x_str, mean_rms_x_str, ok_cnt_y, mean_mean_y_str, mean_rms_y_str )      
                                                                                     
   pyplot.title( title, color=title_color, fontsize=1, y=1.05 )
   pyplot.axis('off')
   for ant in range(0,caldata.n_ant) :
       ant_str = "%03d" % (ant)
       if tiles is not None :
           ant_str = "%03d / %s" % (ant,tiles[ant])
       ax = fig.add_subplot( 16, 8 , (ant+1) )
       
       ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[], aspect=0.2) # 0.3

       if ant == 0 :
           pyplot.text( 0.5, 1.5, title, fontsize=20, color=title_color )           

       mean_x_str = "%.3f" % (mean_x[ant])
       rms_x_str  = "%.3f" % (rms_x[ant])
       if len(mean_x_str) > 10 or len(rms_x_str) > 10 :
           mean_x_str = "NaN"
           rms_x_str  = "NaN (too long)"
       info = "X : %s +/- %s" % (mean_x_str,rms_x_str)

       # to assess quality of calibration solutions we only consider channels known to be good.
       # if some channels (edges of coarse channels) are excluded we don't want to consider them :
       # so we check number of bad channels as (NaN-channels - wrong_channels) to exclude the channels we know are bad :
       ant_good_channels = caldata.n_chan - wrong_channels
       color='black'       
       if np.isnan(mean_x[ant]) or np.isnan(rms_x[ant]) or mean_mean_x>2 or mean_rms_x>1 or out_cnt_x[ant]<(ant_good_channels*0.75) or (out_nan_cnt_x[ant]-wrong_channels)>(0.25*ant_good_channels) or \
          np.isnan(mean_y[ant]) or np.isnan(rms_y[ant]) or mean_mean_y>2 or mean_rms_y>1 or out_cnt_y[ant]<(ant_good_channels*0.75) or (out_nan_cnt_y[ant]-wrong_channels)>(0.25*ant_good_channels) :          
          print "ANT%d is red becuase : mean>2 or rms>1 or is NaN or %d < %d or %d < %d OR %d > %d or %d > %d" % (ant,out_cnt_x[ant],(ant_good_channels*0.75),out_cnt_y[ant],(ant_good_channels*0.75),(out_nan_cnt_x[ant]-wrong_channels),(0.25*ant_good_channels),(out_nan_cnt_y[ant]-wrong_channels),(0.25*ant_good_channels))
          print "Values for X : %.2f , %.2f , %.2f,  %.2f , %d , %d" % (np.isnan(mean_x[ant]),np.isnan(rms_x[ant]),mean_mean_x,mean_rms_x,out_cnt_x[ant], (out_nan_cnt_x[ant]-wrong_channels))
          print "Values for Y : %.2f , %.2f , %.2f,  %.2f , %d , %d" % (np.isnan(mean_y[ant]),np.isnan(rms_y[ant]),mean_mean_y,mean_rms_y,out_cnt_y[ant], (out_nan_cnt_y[ant]-wrong_channels))

          color = 'red'
          
       
       pyplot.text( 0., 0.1, info, fontsize=10, color=color )       

       mean_y_str = "%.3f" % (mean_y[ant])
       rms_y_str  = "%.3f" % (rms_y[ant])
       if len(mean_y_str) > 10 or len(rms_y_str) > 10 :
           mean_y_str = "NaN"
           rms_y_str  = "NaN (too long)"
       info = "Y : %s +/- %s" % (mean_y_str,rms_y_str)

       pyplot.text( 0., 0.45, info, fontsize=10, color=color )
       desc = "%s / %d / %d" % (ant_str,min(out_cnt_x[ant],out_cnt_y[ant]),max(out_nan_cnt_x[ant],out_nan_cnt_y[ant]))
       
       if out_cnt_x[ant] < ant_good_channels*0.75 or ( (out_nan_cnt_x[ant] - wrong_channels) > 0.25*ant_good_channels ) or out_cnt_y[ant]<ant_good_channels*0.75 or ( (out_nan_cnt_y[ant]-wrong_channels) > 0.25*ant_good_channels ) :
           print "ANT%d is red becuase : %d < %d or %d < %d OR %d > %d or %d > %d" % (ant,out_cnt_x[ant],(ant_good_channels*0.75),out_cnt_y[ant],(ant_good_channels*0.75),(out_nan_cnt_x[ant] - wrong_channels),(0.25*ant_good_channels),(out_nan_cnt_y[ant] - wrong_channels),(0.25*ant_good_channels))           
           color = 'red'
           
       pyplot.text( 0.1, 0.8, desc, fontsize=8, fontweight='bold', color=color )
 
   fig.tight_layout(pad=0)


#   fig, axs = pyplot.subplots( 16, 8, figsize=(20, 10))
#   ant=0
#   for ax in axs.flat:
#      ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[], aspect=1)
#      ant_str = "%03d" % (ant)
#      pyplot.text( 0.7, 0.6, ant_str )
#      
#      ant += 1
#      x1, y1 = 0.3, 0.2
#      x2, y2 = 0.8, 0.6
#      ax.plot([x1, x2], [y1, y2], ".")
#      
#      print ant_str
##      ax.annotate("",
##                xy=(x1, y1), xycoords='data',
##                xytext=(x2, y2), textcoords='data',
##                arrowprops=dict(arrowstyle="->",
##                                color="0.5",
##                                shrinkA=5, shrinkB=5,
##                                patchA=None,
##                                patchB=None,
##                                connectionstyle="angle3,angleA=90,angleB=0",
##                                ),
##               )
#   
#   fig.tight_layout(pad=0)

   # pyplot.axis('off')
#   for ant in range(0,caldata.n_ant) :
#       ant_str = "%03d" % (ant)
#       ax = fig.add_subplot( 16, 8 , (ant+1) )
       
       # https://matplotlib.org/gallery/userdemo/connectionstyle_demo.html#sphx-glr-gallery-userdemo-connectionstyle-demo-py
#       ax.set(xlim=(0, 1), ylim=(0, 1), xticks=[], yticks=[], aspect=1)
##       ax.set_xlim(0,1)
##       ax.set_ylim(bottom=0,top=1)
##       ax.axis('off')
       
#       pyplot.text( 0.7, 0.6, ant_str )
#
#       info = "%.4f +/- %.4f" % (mean_x[ant],rms_x[ant])
#       pyplot.text( 0.1, 0.1, info, fontsize=5 )
   
   # pyplot.axis('off')    
   # pyplot.gca().axes.get_xaxis().set_visible(False)
   # pyplot.gca().axes.get_yaxis().set_visible(False)
   
   pngfile = "%s/%s_stat_amp.png" % (outdir,basename)
   if phase > 0 :
       pngfile = "%s/%s_stat_phase.png" % (outdir,basename)
   pyplot.savefig(pngfile)
   print "Saved file %s" % (pngfile)
   if do_show :
       pyplot.show(block=block_image)
   else :
       pyplot.pause(1)
      

                     
if __name__ == "__main__":                     
    calfile="v.bin"
    if len(sys.argv) > 1:
       calfile = sys.argv[1]
    
    calfile_list = None   
    if calfile.find(",") > 0 :
          calfile_list_tmp = calfile.split(",")
          if len(calfile_list_tmp) : 
              calfile_list = []
              for file in calfile_list_tmp :
                 if len(file) > 0 and file.find(".bin") > 0 :
                     calfile_list.append( file )
             
          


    usage="Usage: %prog [options]\n"
    usage+='\tReads calibration solutions from AOCAL binary file (.bin)\n'
    parser = OptionParser(usage=usage,version=1.00)
    parser.add_option('-c','--out_channels',dest="out_channels",default=768, help="Number of output channels [default %default]",type="int")
    parser.add_option('-n','--swap_nans',dest="swap_nans",default=None,help="Swap NaNs to a different value [default %default]")
    parser.add_option('--skip_nans',dest="skip_nans",action="store_true",default=False,help="Skip NaNs  [default %default]")
    parser.add_option('--ch2freq',dest="channels2freq",action="store_true",default=False,help="Convert channels to frequency, requires metafits file to be provided (option --metafits) [default %default]")
    parser.add_option('-m','--metafits',dest="metafits",default=None, help="Metafits file [default %default]")
    # 
    parser.add_option('-a','--ant',dest="ant",default=-1, help="Dump only 1 antenna, if =-1 -> all  [default %default]",type="int")
    parser.add_option('-o','--obsid',dest="obsid",default=-1, help="Obsid  [default %default]",type="int")
    parser.add_option('-p','--phase','--do_phase',dest="do_phase",default=0, help="Dump phase  [default %default , 0 - means dump amplitudes]",type="int")
    parser.add_option('-f','--do_fit','--fit',dest="do_fit",default=0, help="Do fitting  [default %default]",type="int")
    parser.add_option('-i','--do_reim','--reim',dest="do_reim",default=0, help="Save real/imaginary [default %default]",type="int")
    parser.add_option('--do_ampphase','--ampphase',dest="do_ampphase",default=0, help="Save amp/phase [default %default]",type="int")
    parser.add_option('-e','--action','--execute',dest="action",default="dump", help="Execute action [default %default], dump - dumps to txtfiles (if no antenna specified or --ant=-1 -> dumps all), calc_rms, plot, merge")
    parser.add_option('--plot',dest="do_plot",action="store_true",default=False,help="Do plot in addition to other actions [default %default]")
    parser.add_option('--nx',dest="nx",default=16, help="Plot nx  [default %default]",type="int")
    parser.add_option('--ny',dest="ny",default=8,  help="Plot ny  [default %default]",type="int")
    parser.add_option('--min_y',dest="min_y",default=0, help="Min Y value  [default %default]",type="float")
    parser.add_option('--max_y',dest="max_y",default=1.5, help="Max Y value  [default %default]",type="float")
    parser.add_option('--edge_channels','--wrong_channels',dest="wrong_channels",default=192, help="Number of wrong channels (already excised at cotter stage, assuming edgewidth=160 -> 4*24*4 = 192  [default %default]",type="int")
    parser.add_option('--out_basename',dest="out_basename",default="calsolutions", help="Outfile basename  [default %default]")

    parser.add_option('--channels_list','--channels',dest="channels_str",default=None, help="String with list of coarse channels  [default %default]")
    parser.add_option('--n_coarse_channels',dest="n_coarse_channels",default=24, help="Number of coarse channels [default %default]",type="int")
    parser.add_option('--n_fine_channels_per_coarse',dest="n_fine_channels_per_coarse",default=32, help="Number of fine channels per coarse [default %default]",type="int")
    parser.add_option('--coarse_channel_width',dest="coarse_channel_width",default=1.28, help="Width of a single coarse channel [default %default]",type="float")
  
    parser.add_option('--merged_file','--merge',dest="merged_bin_file",default=None,help="Merge a coma separated list of .bin files (from picket fence observations) to a single file , use option value to provide a name of output file [default %default]")
    parser.add_option('--split_file','--split','--split_by_n_channels',dest="split_by_n_channels",default=-1,help="Allows to split AO .bin file into serveral consisting of specified number of channels [default %default]",type="int")
    parser.add_option('--outname','--outbasename','--outfile',dest="outbasename",default=None,help="Output file name base [default %default]");
# def average_channels( binfile, avg_n_channels, outfile ) :
    parser.add_option('--avg_n_channels','--average_channes','--avg_n',dest="average_n_channels",default=-1,help="Average N channels [default %default and <0 -> no averaging]",type="int");    
    parser.add_option('--final_n_channels','--final_channes','--final_ch',dest="final_n_channels",default=-1,help="Make sure final number of channels is as specified [default %default and <0 -> no averaging]",type="int");
    
    (options,args)=parser.parse_args(sys.argv[1:])
    if options.final_n_channels > 0 :
       options.action = "final_n_channels"
    if options.average_n_channels > 0 :
       options.action = "average_n_channels"
    if options.split_by_n_channels > 0 :
       options.action = "split_by_n_channels"
    
    COARSE_CHANNEL_WIDTH = options.coarse_channel_width
    N_FINE_CHANNELS_PER_COARSE = options.n_fine_channels_per_coarse
    FINE_CHANNEL = COARSE_CHANNEL_WIDTH / N_FINE_CHANNELS_PER_COARSE    
    
    if options.merged_bin_file is not None :
        options.action = "merge"    

    band      = COARSE_CHANNEL_WIDTH
    half_band = COARSE_CHANNEL_WIDTH / 2.00

    print "#######################################################"
    print "PARAMETERS :"
    print "#######################################################"
    print "calfile      = %s" % (calfile)
    print "obsid        = %d" % (options.obsid)
    print "out_basename = %s" % (options.out_basename)
    print "action       = %s" % (options.action)
    print "plot         = %s (%d x %d)" % (options.do_plot,options.nx,options.ny)
    print "swap_nans    = %s" % (options.swap_nans)
    print "skip_nans    = %s" % (options.skip_nans)
    print "ant          = %d" % (options.ant)
    print "do_phase     = %d" % (options.do_phase)
    print "do_fit       = %d" % (options.do_fit)
    print "do_reim      = %d" % (options.do_reim)
    print "do_do_ampphase = %d" % (options.do_ampphase)
    print "out_channels = %d" % (options.out_channels)
    print "channels2freq = %s" % (options.channels2freq)
    print "metafits     = %s" % (options.metafits)
    print "band / half_band = %.4f / %.4f MHz" % (band,half_band)
    print "Y values range = %.4f - %.4f" % (options.min_y,options.max_y)
    print "wrong channels = %d" % (options.wrong_channels)
    print "channels_str   = %s" % (options.channels_str)
    print "COARSE_CHANNEL_WIDTH = %.2f MHz" % (COARSE_CHANNEL_WIDTH)
    print "N_FINE_CHANNELS_PER_COARSE = %d" % (N_FINE_CHANNELS_PER_COARSE)
    print "FINE_CHANNEL               = %.4f MHz" % (FINE_CHANNEL)
    print "merged_bin_file = %s" % (options.merged_bin_file)
    print "split_by_n_channels = %d" % (options.split_by_n_channels)
    print "average N channels  = %d" % (options.average_n_channels)
    print "Final N channels    = %d" % (options.final_n_channels)
    print "out file basename   = %s" % (options.outbasename)
    print "#######################################################"

    if options.action == "dump" :     
        dump_ao_calsolutions( calfile=calfile, options=options, ant=options.ant, do_phase=options.do_phase, do_fit=options.do_fit, do_reim=options.do_reim, do_ampphase=options.do_ampphase,out_basename_param=options.out_basename,  channels_str=options.channels_str )    
    elif options.action == "calc_rms" or options.action == "rms" :
        (ok_cnt_x,mean_x,rms_x,out_cnt_x,out_nan_cnt_x,mean_mean_x,rms_mean_x,mean_rms_x,rms_rms_x,
         ok_cnt_y,mean_y,rms_y,out_cnt_y,out_nan_cnt_y,mean_mean_y,rms_mean_y,mean_rms_y,rms_rms_y) = calc_mean_rms( calfile, do_fit=options.do_fit, do_phase=options.do_phase )         
    elif options.action == "plot" :
        plotcal( calfile, nx=options.nx, ny=options.ny, do_fit=options.do_fit, phase=options.do_phase, min_y=options.min_y, max_y=options.max_y, wrong_channels=options.wrong_channels, obsid=options.obsid ) 
    elif options.action == "merge" and calfile_list is not None :
        merge_bin_files( calfile_list, options.merged_bin_file )       
    elif options.split_by_n_channels > 0 :
        # def split_bin_file( binfile, split_by_n_channels, outfilebase ):
        split_bin_file( calfile, options.split_by_n_channels, options.outbasename )
    elif options.average_n_channels > 0 :
        # def average_channels( binfile, avg_n_channels, outfile ) :
        average_channels( calfile, options.average_n_channels, options.outbasename )
    elif options.final_n_channels > 0 :
        final_n_channels( calfile, options.final_n_channels, options.outbasename )
    else :   
        print "ERROR : unknown action = %s" % (options.action)
        
        
# test in /home/msok/Desktop/MWA/doc/ASVO/data/2013/201306/1056386176/Chris_Jordan_RTS/20180801/1050871872/rts/RTS_NATIVE/ :
# 20180921 - latest test to save RTS cal. solutions to aocal file :
#   rtscal = rtsfile_bandpass('1050871872.metafits',aocal_format=False,debug=1)
#   rtscal.tofile('1050871872_rts_bandpass.bin')
   
#   rtscal = rtsfile('1050871872.metafits',aocal_format=False)
#   rtscal.tofile('1050871872_rts.bin')
#
#   rtscal = rtsfile('1056386176.metafits',aocal_format=False)
#   rtscal.tofile('1056386176_rts.bin')

#    rtscal = rtsfile_bandpass('1056386176.metafits',aocal_format=False)
#    rtscal.tofile('1056386176_bandpass.bin')
