#!/usr/bin/env python

import sys, struct, os.path
import numpy as np


def load_cli_data(filename, rate, nseg=-1, nskip=0, nchan=1, seg_byte=4096, size_byte=2, meta=False, dual=False, verbose=1):
    '''
    load data saved by bladeRF-cli

    <mandatory, no defaults>
    filename:: filanem of the data file
    rate:: sampling rate in Hz

    <optional>
    nseg:: how many segments to read. -1 means everything in the file
    nskip:: number of segments to skip 
    nchan:: vector length == number of channenls of the FFT
    seg_byte:: number of bytes in a segment (after each header)
    size_byte:: number of bytes per number (I or Q)
        4 --> for float32
        1 --> for char
        2 --> native DAC format
    meta:: True is metadata is saved; False if no metadata
    dual:: True if data is savesd in mimo mode (channel 1 and 2 interleaved)
           False if only one channel is saved
    verbose:: 0 --> no details
              1 --> basic info
              2 --> additional info
            (older style verbose=True is converted to verbose=2; verbose=False is converted to verbose=0)
    '''
    #print('in load_cli_data, rate:', rate)

    header_len = 16 # sc16_q11_meta spec mentions 16 bytes as metadata
    fmt0 = ''
    dt = 1./rate

    if (meta == False):
        if (nseg == -1):
            nitem = -1
        else:
            nitem = seg_byte     # number of bytes to read, no header
        skip = nskip * seg_byte        # number of bytes to skip
    else:
        nitem = seg_byte    # number of bytes to read after each header
        skip = nskip * (header_len + seg_byte)    # number of bytes to skip

    if (size_byte == 1):
        fmt = 'b'   # signed character
        rescale = 1./127.
    elif (size_byte == 4):
        fmt = 'f'   # float
        rescale = 1.
    elif (size_byte == 2):  # the native DAC format 'SC16Q11'
        fmt = 'h'   # 2-byte short (signed integer)
        fmt0 = '<'  # little-endian
        rescale = 1.

    if (verbose == False):
        verbose = 0
    elif (verbose == True):
        verbose = 2


    fs = os.path.getsize(filename)
    nsamp = fs // (size_byte * 2)

    if (verbose > 0):
        print('== parameters for load_cli_data ==')
        print(' total number of samples:', nsamp)
        print(' total number of segments:', fs//(seg_byte))
        print(' byte per I/Q:', size_byte)
        print(' bytes in a segment:', seg_byte)
        print(' num of channel:', nchan)
        print(' num of seg to read:', nseg)
        print(' num of seg to skip:', nskip)
        print(' dual channel:', dual)
        print(' sampling rate (MSps):', rate/1e6)
        print('=====')


    fh = open(filename, 'rb')
    fh.seek(skip)


    if (dual):
        data0 = []
        data1 = []
        nInp = 2
    else:
        data = []
        nInp = 1

    # avoid a problem of splitting inputs
    while (seg_byte < nchan*size_byte*2*nInp):
        seg_byte *= 2
        nseg = nseg//2
        print('warning: seg_byte changed to', seg_byte)

    seg_dt = nitem / (size_byte*2) * dt
    seg_dt /= nInp
    #print('seg_dt', seg_dt)
    if (nitem > -1):
        nitem = seg_byte

    chunk = seg_byte
    while (chunk < 1000000):
        chunk *= 10
    seg_chunk = chunk // seg_byte

    timestamp = []
    seg_count = 0
    while(seg_count < nseg or nseg==-1):
        if (seg_count%seg_chunk == 0 and verbose>1):
            print('  segment reached:', seg_count+nskip)

        # tsec may be extracted from meta data when it is implemented
        # for now, set tsec according to the segment number
        tsec = (nskip + seg_count) * seg_dt
        if (meta):
            header_str = fh.read(header_len)
            if (len(header_str) < header_len):
                print('header incomplete segment:', seg_count)
                break

        try:
            #seg_data = np.fromfile(fh, dtype=np.complex64, count=nitem, sep='', offset=offset)
            buf = fh.read(nitem)
            nbyte = len(buf)
            if (len(buf) < nitem or len(buf)==0):
                print('not enough data left:', seg_count)
                break
            ncplx = nbyte // (size_byte * 2)
            #seg_data = np.asarray(struct.unpack('%d%s'%(ncplx*2, fmt), buf)).astype(np.float32).view(np.complex64)
            seg_data = np.asarray(struct.unpack('%s%d%s'%(fmt0, ncplx*2, fmt), buf)).astype(np.float32).view(np.complex64)
        except:
            print('error. skip this segment.')
            break

        if (meta):
            if (len(seg_data)//nchan < nspec):
                print('incomplete segment, skip:', seg_count)
                break

        if (dual):
            ntime = ncplx // 2
            seg_data0 = seg_data[0:ntime*2:2]
            seg_data0 = seg_data0[:seg_data0.shape[0]//nchan*nchan].reshape((-1,nchan))
            data0.append(seg_data0)
            seg_data1 = seg_data[1:ntime*2:2]
            seg_data1 = seg_data1[:seg_data1.shape[0]//nchan*nchan].reshape((-1,nchan))
            data1.append(seg_data1)
        else:
            ntime = ncplx
            seg_data = seg_data[:seg_data.shape[0]//nchan*nchan].reshape((-1,nchan))
            data.append(seg_data)

        nspec = ntime // nchan # number of spectra in this segment (or whole data when nitem==-1)
        seg_time = tsec + np.arange(nspec)*dt*nchan # each spectra time is increment by dt*nchan
        #if (seg_count<20):
        #    print('ncplx, ntime, nspec, nchan, tsec', ncplx, ntime, nspec, nchan, tsec)
        timestamp.append(seg_time)


        seg_count += 1

    fh.close()
    #print('total num of segments:', seg_count)

    timestamp = np.asarray(timestamp).reshape(-1)
    #timestamp /= 2. # problem with timing, need to be fixed later

    if (dual):
        data0 = np.asarray(data0).reshape((-1, nchan))
        data0 *= rescale
        data1 = np.asarray(data1).reshape((-1, nchan))
        data1 *= rescale
        return data0, data1, timestamp
    else:
        data = np.asarray(data).reshape((-1, nchan))
        data *= rescale
        return data, timestamp

