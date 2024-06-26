#!/usr/bin/env python3

import os.path
from read_meta3 import *
import matplotlib.pyplot as plt

inp = sys.argv[0:]
pg  = inp.pop(0)


#-- defaults --
size_byte = 4   # default to 4 for now
meta   = True   # whether binary file contains metadata
inline = True   # whether metadata is embedded in binary file
nchan  = 16     # number channels in the FFT
nseg   = -1     # number of segments to read (-1 for whole file)
nskip  = 0      # number of segments to skip
seg_byte = 8000000 # number of bytes per segment (after each header)
is_raw = False  # whether the input file is the raw data
verbose = False
vector = True


#-- usage --
usage = '''
make a waterfall plot of the binary data
saved from GnuRadio 'file sink' or 'file_meta_sink' block
the default is a 'file_meta_sink' with 'inline' metadata

syntax:
    %s <binary_file> [options]

where options are:
    --no-meta       the binary file contains no metadata
    -n N            read N segments instead of the whole file
    -s N            number of segments to skip
    -c nchan        set number of channels (default: %d)
    --size SIZE     number of bytes per I or Q variable (default: %d)
                    (i.e. 4 bytes or 1 byte)
    --raw           file stores raw data before FFT
                    an additional FFT plot is produced
    --seg_byte SEG_BYTE
                    sets how many bytes to read in between metadata blocks
                    (default: %d)
    -v              print out metadata for each block


''' % (pg, nchan, size_byte, seg_byte)


#-- parsing arguements --
if (len(inp) < 1):
    sys.exit(usage)

while (inp):
    k = inp.pop(0)
    if (k == '--no-meta'):
        meta = False
    elif (k == '-n'):
        nseg = int(inp.pop(0))
    elif (k == '-s'):
        nskip = int(inp.pop(0))
    elif (k == '-c'):
        nchan = int(inp.pop(0))
    elif (k == '-v'):
        verbose = True
    elif (k == '--seg_byte'):
        seg_byte = int(inp.pop(0))
    elif (k == '--raw'):
        is_raw = True
        vector = False
    elif (k == '--size'):
        size_byte = int(inp.pop(0))
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        filename = k
        if not(os.path.isfile(filename)):
            sys.exit('error finding data: %s' % filename)


augname = 'skip%d' % nskip
if (nseg == -1):
    augname += '_all'
else:
    augname += '_seg%d' % nseg

data, timestamp = load_blade_data(filename, nseg=nseg, nskip=nskip, nchan=nchan, seg_byte=seg_byte, size_byte=size_byte, meta=meta, verbose=verbose, vector=vector)
#print(data.shape, timestamp.shape)
#print(timestamp[:50])
if (False): # debug
    seg_len = int(1000000 / nchan) 
    for i in range(10):
        j0 = i * seg_len
        j1 = (i+1) * seg_len - 1
        print(i, timestamp[j0], timestamp[j1])


fig, ax = plt.subplots(1,1,figsize=(12,6))
chan = np.arange(nchan)
#s = ax.pcolormesh(timestamp, chan, 20.*np.log10(np.abs(data.T)))
#s = ax.imshow(20.*np.log10(np.abs(data.T)), origin='lower', extent=[timestamp[0], timestamp[-1], chan[0], chan[-1]])
s = ax.imshow(np.abs(data.T), origin='lower', extent=[timestamp[0], timestamp[-1], chan[0], chan[-1]])
ax.set_aspect('auto')
ax.set_xlabel('time (sec)')
if (is_raw):
    ax.set_ylabel('time idx')
else:
    ax.set_ylabel('chan')
cb = plt.colorbar(s, ax=ax)
#cb.set_label('Power (dB)')
cb.set_label('Amplitude')
ax.set_title('%s, %s' % (filename, augname))

png = '%s.%s.png' % (filename, augname)
fig.savefig(png)
plt.close(fig)

if (is_raw):
    data = np.fft.fft(data, axis=1)
    data = np.fft.fftshift(data, axes=1)

    fig, ax = plt.subplots(1,1,figsize=(12,6))
    s = ax.imshow(np.abs(data.T), origin='lower', extent=[timestamp[0], timestamp[-1], chan[0], chan[-1]])
    ax.set_aspect('auto')
    ax.set_xlabel('time (sec)')
    ax.set_ylabel('chan')
    cb = plt.colorbar(s, ax=ax)
    cb.set_label('Amplitude')
    ax.set_title('%s, %s, FFT' % (filename, augname))

    png = '%s.%s.fft.png' % (filename, augname)
    fig.savefig(png)
    plt.close(fig)

    if (False): # skip this part
        fig, ax = plt.subplots(1,1,figsize=(12,6))
        rfi = np.zeros(nchan, dtype=bool)
        ## add RFI channel = True here
        # rfi[ch1:ch2] = True
        chavg = data[:,~rfi].mean(axis=1)
        ax.plot(timestamp, chavg.real, label='real')
        ax.plot(timestamp, chavg.imag, label='imag')
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('Real/Imag (arbitrary)')
        ax.set_title('%s, %s, ch-avg' % (filename, augname))
        ax.legend()

        png = '%s.%s.ch-avg.png' % (filename, augname)
        fig.savefig(png)
        plt.close(fig)

