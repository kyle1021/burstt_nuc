#!/usr/bin/env python

from gnuradio import gr,blocks
import pmt
import sys, struct, os.path
from gnuradio.blocks import parse_file_metadata
import numpy as np

print_out = True


def load_blade_data(filename, nseg=-1, nskip=0, nchan=64, seg_byte=8000000, size_byte=4,
        meta=True, inline=True, offset=1, verbose=1, vector=1):
    '''
    load data saved by GnuRadio
    filename:: filanem of the data file
    nseg:: how many segments to read. -1 means everything in the file
    nskip:: number of segments to skip 
    nchan:: vector length == number of channenls of the FFT
    seg_byte:: number of bytes in a segment (after each header)
    size_byte:: number of bytes per number (I or Q)
        4 --> for float32
        1 --> for char
    meta:: True is metadata is saved; False if no metadata
    inline:: True if metadata is contained in the same file
             False if metadata is in a separate file 'filename.hdr'
    offset:: number of bytes between header and data
    verbose:: print out each header info
                0 --> no details (aka False)
                1 --> basic info
                2 --> more details (aka True)
    vector:: True --> data is FFT vectors, rx_rate in metadata is slowed
             False --> data is raw stream, rx_rate is original
    '''

    header_len = parse_file_metadata.HEADER_LENGTH
    fmt0 = ''

    if (isinstance(verbose, bool) and verbose == True):
        verbose = 2
    elif (verbose == False):
        verbose = 0
    if (verbose > 1):
        print_out = True
    else:
        print_out = False

    if (meta == False or inline==False):
        if (nseg == -1):
            nitem = -1
        else:
            nitem = nseg * seg_byte     # number of bytes to read, no header
        skip = nskip * seg_byte        # number of bytes to skip
        offset = 0
    else:
        nitem = seg_byte    # number of bytes to read after each header
        skip = nskip * (header_len + offset + seg_byte)    # number of bytes to skip

    if (size_byte == 1):
        fmt = 'b'   # signed character
        rescale = 1./127.
    elif (size_byte == 4 or size_byte == 8):
        fmt = 'f'   # float
        rescale = 1.
    elif (size_byte == 2):  # the native DAC format 'SC16Q11'
        fmt = 'h'   # 2-byte short (signed integer)
        fmt0 = '<'  # little-endian
        rescale = 1.


    data = []
    timestamp = []

    fs = os.path.getsize(filename)
    nsamp = fs // (size_byte * 2)
    if (verbose > 0):
        print('== parameters for read_meta2 ==')
        print(' total number of samples:', nsamp)
        print(' total number of segments:', (fs//seg_byte))
        print(' byte per I/Q:', size_byte)
        print(' bytes in a segment:', seg_byte)
        print(' num of channel:', nchan)
        print(' num of seg to read:', nseg)
        print(' num of seg to skip:', nskip)
        print('=====')

    fh = open(filename, 'rb')
    fh.seek(skip)

    seg_count = 0
    while(seg_count < nseg or nseg==-1):
        if (seg_count%100 == 0):
            if (verbose > 0):
                print('  segment reached:', seg_count+nskip)

        if (meta):
            header_str = fh.read(header_len)
            if (len(header_str) < header_len):
                break
            header = pmt.deserialize_str(header_str)
            header_info = parse_file_metadata.parse_header(header, print_out)
            tsec = header_info['rx_time']
            rate = header_info['rx_rate']
            dt = 1./rate
            if (vector == False):
                dt *= nchan
            nspec = seg_byte // (size_byte * 2 * nchan) # number of spectra in this segment
            seg_time = tsec + np.arange(nspec)*dt
            #print('tsec=', tsec, 'dt=', dt, 'seg_time:', seg_time[0], seg_time[-1], 'nspec=',nspec)

        try:
            #seg_data = np.fromfile(fh, dtype=np.complex64, count=nitem, sep='', offset=offset)
            if (offset > 0):
                dummy = fh.read(offset)
            buf = fh.read(nitem)
            nbyte = len(buf)
            ncplx = nbyte // (size_byte * 2)
            #seg_data = np.asarray(struct.unpack('%d%s'%(ncplx*2, fmt), buf)).astype(np.float32).view(np.complex64)
            seg_data = np.asarray(struct.unpack('%s%d%s'%(fmt0, ncplx*2, fmt), buf)).astype(np.float32).view(np.complex64)

            if (meta==False):
                seg_data = seg_data[:seg_data.shape[0]//nchan*nchan].reshape((-1,nchan))
                data.append(seg_data)
                timestamp.append(np.arange(seg_data.shape[0], dtype=float))
                break
            else:
                if (len(seg_data)//nchan < nspec):
                    break
                #seg_data = seg_data.reshape((-1, nchan))
                seg_data = seg_data[:seg_data.shape[0]//nchan*nchan].reshape((-1,nchan))
                data.append(seg_data)
                timestamp.append(seg_time)
        except:
            print('error. skip this segment.')

        seg_count += 1

    fh.close()
    if (verbose > 0):
        print('total num of segments:', seg_count)

    data = np.asarray(data).reshape((-1, nchan))
    data *= rescale
    timestamp = np.asarray(timestamp).reshape(-1)
    #timestamp /= 2. # problem with timing, need to be fixed later
    return data, timestamp

