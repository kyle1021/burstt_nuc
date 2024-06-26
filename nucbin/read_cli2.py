#!/usr/bin/env python

import sys, struct, os.path
import numpy as np
import loadh5 as lh
from scipy.linalg import svd, eigh
import multiprocessing as mp
import time
from calibrate_func import *


def load_cli_hdr(hdrname, C0=0, nSamp=16384000, verbose=1):
    '''
    with the help of header file (.h5), load the bladeRF-meta data with better precision
    the return data will be a masked array, where dropped data is filled with zeros and masked

    input:
    hdrname::       header file name
    C0::            retrieving data from clock count C0
                    note: the file begins with C0=0, regardless of actual fpga_tick saved
    nSamp::         num of samples to retrieve

    return:
    data::          masked array, shape(nSamp,nInp)
                    nInp is the number of input channels (1 or 2), determined from the header file
    clock::         uint64 array, shape(nSamp)
                    clock count starting from 0 at the beginning of the file
    '''

    seg_samp = 32768    # num of samples per seg

    #name = hdrname.strip('.h5')
    name = hdrname.replace('.h5', '')
    fbin = '%s.meta' % name

    if (not os.path.isfile(hdrname)):
        print('error finding header file:', hdrname)
        return None
    if (not os.path.isfile(fbin)):
        print('error finding bindary file:', fbin)
        return None

    attrs = lh.getAttrs(hdrname)
    nInp = int(attrs.get('nInp'))
    nclk = seg_samp // nInp
    if (nInp == 2):
        dual = True
    elif (nInp == 1):
        dual = False

    tick = lh.getData(hdrname, 'fpga_tick')    # the free running fpga tick of each segment
    tick -= tick[0]                         # removing the starting tick
    toff = tick.astype(np.int64) - C0
    i0 = np.argmin(np.abs(toff))    # i0 is the segment id closes to C0
    if (toff[i0]>0 and i0>0):       # tick[i0] should be smaller than C0
        i0 -= 1
    # toff[i0] should be negative now.
    noff = int(np.ceil(np.abs(toff[i0])/nclk))  # extra num of segments needed to fill the nSamp requirement

    nseg = int(np.ceil(nSamp/nclk)) + noff  # get slightly more data to ensure we can save nSamp
    data, clock = load_cli_data2(fbin, nseg=nseg, nskip=i0, dual=dual, verbose=verbose)
    i1 = C0 - int(tick[i0]) # begin of return data
    i2 = i1 + nSamp         # end of return data
    if (i2 > len(clock)):
        i2 = len(clock)

    return data[i1:i2], clock[i1:i2]



def load_cli_data2(filename, nseg=1000, nskip=0, seg_byte=131072, size_byte=2, meta=True, dual=True, verbose=1):
    '''
    load data saved by bladeRF-meta
    when using meta=True, non-consecutive data is filled with zeros
    note that returning clock does not require the knowledge of sample_rate

    return data:
        data[nSamp, nInp], masked
        clock[nSamp]
        

    <mandatory, no defaults>
    filename:: filanem of the data file

    <optional>
    nseg:: how many segments to read. -1 means everything in the file
    nskip:: number of segments to skip 
    seg_byte:: number of bytes in a segment (after each header)
    size_byte:: number of bytes per number (I or Q)
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

    header_len = 11 # 64bit timestamp + 16bit actual_count + 8bit overrun
    fmt0 = ''

    if (meta == False):
        if (nseg == -1):
            nitem = -1
        else:
            nitem = seg_byte     # number of bytes to read, no header
        skip = nskip * seg_byte        # number of bytes to skip
    else:
        nitem = seg_byte    # number of bytes to read after each header
        skip = nskip * (header_len + seg_byte)    # number of bytes to skip
        skip += 11          # there is a dummy header at the beginning, which has no sample following it

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
    maxseg = fs // seg_byte

    if (verbose > 0):
        print('== parameters for load_cli_data ==')
        print(' total number of samples:', nsamp)
        print(' total number of segments:', maxseg)
        print(' byte per I/Q:', size_byte)
        print(' bytes in a segment:', seg_byte)
        print(' num of seg to read:', nseg)
        print(' num of seg to skip:', nskip)
        print(' dual channel:', dual)
        print('=====')


    fh = open(filename, 'rb')
    if (meta):
        fh.seek(11)
        tick0 = struct.unpack('<Q', fh.read(8))[0]  # get the first tick in the file
    else:
        tick0 = 0
    if (verbose>1):
        print('file starting tick:', tick0)

    fh.seek(skip)   # go to the beginning of the segment (skipping nskip segments)


    if (nitem > -1):
        nitem = seg_byte
    else: # (nitem = nseg) == -1
        nitem = seg_byte
        nseg = (maxseg - nskip)
        

    if (dual):
        nInp = 2
    else:
        nInp = 1

    seg_samp = seg_byte // (size_byte * 2)  # num of samples in a segment
    nclk = seg_samp // nInp                 # number of clocks in a segment
    ntotal = nseg * nclk                    # total length of data to return
    clock = np.arange(ntotal, dtype=np.uint64)
    data = np.ma.zeros((ntotal, nInp), dtype=np.complex64)
    data.mask = np.ones_like(data, dtype=bool)  # all masked by default. unmask when valid data is entered.


    seg_count = 0
    while(seg_count < nseg):
        progress = 100. * seg_count/nseg
        if (seg_count%1000 == 0 and verbose>1):
            print('  %2d%%, segment reached:' % progress, seg_count+nskip)

        tick_id = (nskip + seg_count) * nclk    # the intended starting tick of this segment
        if (meta):
            try:
                tick = struct.unpack('<Q', fh.read(8))[0]
                acct = struct.unpack('<H', fh.read(2))[0]
                ovrn = struct.unpack('<B', fh.read(1))[0]
                if (ovrn>0):
                    print(tick, acct, ovrn)
                tick -= tick0           # removing tick offset between files
            except:
                print('error, skip this segment')
                break

        else:
            tick = tick_id
            acct = seg_samp
            ovrn = 0

        if (seg_count == 0):
            tick1 = tick        # staring tick of first read segment
            clock += tick1
        itk = tick - tick1      # starting index of each segment
        nvalid = acct // nInp   # valid number of samples for each Inp
        if (itk >= ntotal):     # segment outside return range
            break

        try:
            buf = fh.read(nitem)
            nbyte = len(buf)
            if (len(buf) < nitem or len(buf)==0):
                print('not enough data left:', seg_count)
                break
            ncplx = nbyte // (size_byte * 2)
            seg_data = np.asarray(struct.unpack('%s%d%s'%(fmt0, ncplx*2, fmt), buf)).astype(np.float32).view(np.complex64)
            seg_data = seg_data.reshape((-1,nInp))
            if (itk+nvalid < ntotal):
                iend = itk + nvalid
            else:
                iend = ntotal
                nvalid = ntotal - itk
                
            data[itk:iend,:] = seg_data[:nvalid,:]
            data.mask[itk:iend,:] = False
            
        except:
            print('error. skip this segment.')
            break

        seg_count += 1

    fh.close()
    print('total num of segments:', seg_count)

    return data, clock


def load_cli_data(filename, rate, nseg=-1, nskip=0, nchan=1, seg_byte=131072, size_byte=2, meta=False, dual=False, verbose=1):
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

    header_len = 11 # 64bit timestamp + 16bit actual_count + 8bit overrun
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
        skip += 11          # there is a dummy header at the beginning, which has no sample following it

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


    timestamp = []
    if (dual):
        data0 = []
        data1 = []
        nInp = 2
    else:
        data = []
        nInp = 1

    if (meta):
        acount = []     # both have one entry per segment
        overrun = []

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
    while (chunk < 100000000):
        chunk *= 10
    seg_chunk = chunk // seg_byte

    seg_count = 0
    while(seg_count < nseg or nseg==-1):
        if (seg_count%seg_chunk == 0 and verbose>1):
            print('  segment reached:', seg_count+nskip)

        # tsec may be extracted from meta data when it is implemented
        # for now, set tsec according to the segment number
        tsec = (nskip + seg_count) * seg_dt
        if (meta):
            tick = struct.unpack('<Q', fh.read(8))[0]
            acct = struct.unpack('<H', fh.read(2))[0]
            ovrn = struct.unpack('<B', fh.read(1))[0]
            if (ovrn>0):
                print(tick, acct, ovrn)

            dt = 1./rate
            tsec = tick * dt
            nspec = seg_byte // (size_byte * 2 * nchan) # number of spectra in this segment
            seg_time = tsec + np.arange(nspec)*dt
            #print('tsec=', tsec, 'dt=', dt, 'seg_time:', seg_time[0], seg_time[-1], 'nspec=',nspec)
            acount.append(acct)
            overrun.append(ovrn)

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
    if (meta):
        acount = np.asarray(acount)
        overrun = np.asarray(overrun)

    if (dual):
        data0 = np.asarray(data0).reshape((-1, nchan))
        data0 *= rescale
        data1 = np.asarray(data1).reshape((-1, nchan))
        data1 *= rescale
        if (meta):
            return data0, data1, timestamp, acount, overrun
        else:
            return data0, data1, timestamp
    else:
        data = np.asarray(data).reshape((-1, nchan))
        data *= rescale
        if (meta):
            return data, timestamp, acount, overrun
        else:
            return data, timestamp


def extractMeta(fname, nProc=1, verbose=1):
    '''
    extract the full metadata from saved data
    verbose=1       print out some info
    verbose=0       turn off info

    nProc   number of parallel processes, max mp.cpu_count()
    '''

    mProc = mp.cpu_count()
    if (nProc > mProc):
        nProc = mProc
        print('extractMeta: reset nProc from %d to max:%d'%(nProc,mProc))

    header_len = 11
    seg_byte = 131072

    fs = os.path.getsize(fname)
    #fh = open(fname, 'rb')
    #dummy = fh.read(header_len) # first meta is dummy
    pos0 = header_len
    max_seg = (fs-pos0)//seg_byte
    args = []
    if (nProc == 1):
        args.append([0, fname, pos0, max_seg, seg_byte, verbose])
    else:
        exSeg = max_seg % nProc
        nSeg = int((max_seg-exSeg)/nProc)  # num of seg per process
        if (exSeg > 0):
            nSeg += 1    # make the distribution more even, last one may be shorter
        for i in range(nProc):
            seg1 = int(i*nSeg)
            seg2 = min(max_seg, seg1+nSeg)
            pos1 = pos0 + int(seg1 * (seg_byte+header_len))
            args.append([i, fname, pos1, seg2-seg1, seg_byte, verbose])

    #print('args:', args)
    #sys.exit()


    with mp.Pool(nProc) as p:
        chunk_results = p.starmap(readChunk, args)

    clock = []
    acount = []
    overrun = []
    for cr in chunk_results:
        clock.extend(cr[0])
        acount.extend(cr[1])
        overrun.extend(cr[2])

    return np.asarray(clock, dtype=np.uint64), np.asarray(acount, dtype=np.uint16), np.asarray(overrun, dtype=np.uint8)


def readChunk(idx, fname, pos, nSeg, seg_byte, verbose):
    t1 = time.time()
    clock   = []
    acount  = []
    overrun = []
    progress = 0.
    seg2prompt = 720000 # roughly 300sec under 40MHz dual input
    if (verbose>0):
        print('readChunk: idx, pos, nSeg, seg_byte, verbose:', idx, pos, nSeg, seg_byte, verbose)

    fh = open(fname, 'rb')
    fh.seek(pos)
    for i in range(nSeg):
        progress = 100. * i / nSeg
        if (i % seg2prompt == 0):
            if (verbose>0):
                print('progress%d: %3d%%, seg:%d' % (idx, progress, i), time.asctime())
            
        try:
            tick = struct.unpack('<Q', fh.read(8))[0]
            acct = struct.unpack('<H', fh.read(2))[0]
            ovrn = struct.unpack('<B', fh.read(1))[0]
            clock.append(tick)
            acount.append(acct)
            overrun.append(ovrn)
            if (ovrn>0 and verbose>1):
                print(tick, acct, ovrn)

            fh.seek(fh.tell()+seg_byte)

        except:
            print('finished.')
            break

    fh.close()
    t2 = time.time()
    if (verbose>0):
        print('idx:', idx, 'time:', t2-t1, 'sec')

    return [clock, acount, overrun]


def logParse(flog):
    '''
    parse the bladeRF-cli command log (by e.g. run_blade_meta.pl)
    and return a dictionary that can be saved as header
    '''

    fh = open(flog, 'r')
    lines = fh.readlines()

    freq = []
    rate = []
    gain = []
    chan = []
    dur  = 0.
    for l in lines:
        if ('RX1 Frequency' in l):
            ff = float(l.split()[2])
            freq.append(ff)
        elif ('RX2 Frequency' in l):
            ff = float(l.split()[2])
            freq.append(ff)
        elif ('RX1 sample rate' in l):
            rr = float(l.split()[9])
            rate.append(rr)
        elif ('RX2 sample rate' in l):
            rr = float(l.split()[9])
            rate.append(rr)
        elif ('RX1 overall gain' in l):
            gg = float(l.split()[5])
            gain.append(gg)
        elif ('RX2 overall gain' in l):
            gg = float(l.split()[5])
            gain.append(gg)
        elif ('Channels:' in l):
            chan = l.strip(',').split()[1:]
        elif ('Duration:' in l):
            dur = float(l.split()[1])

    attrs = {
        'freq':freq,
        'rate':rate,
        'gain':gain,
        'channels':chan,
        'nInp':len(chan),
        'duration':dur
            }

    return attrs

