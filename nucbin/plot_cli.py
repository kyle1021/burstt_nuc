#!/usr/bin/env python3

import os.path
from read_cli import *
import matplotlib.pyplot as plt

inp = sys.argv[0:]
pg  = inp.pop(0)


#-- defaults --
size_byte = 2   # default to 4 for now
seg_byte = 4096 # number of bytes per segment (after each header)
nchan  = 1024   # number channels in the FFT
nseg   = -1     # number of segments to read (-1 for whole file)
nskip  = 0      # number of segments to skip
meta   = False  # whether binary file contains metadata
verbose = False
vector = False
dual = True     # data file is dual-channel
fcen = 0.       # central freq in MHz


#-- usage --
usage = '''
make a waterfall plot of the binary data
saved by bladeRF-cli

syntax:
    %s <binary_files> -r <samp_rate> [options]

    binary_files:: filename of the binary data
                    multiple files can be given
    samp_rate:: sampling rate in Msps (no default)
                this controls the timestamps dt

options are:
    --single        single-channel data file
                    (default: dual-channel data file)

    --fcen MHz      central frequency in MHz
                    (default: 0)

    -n N            read N segments instead of the whole file
                    (default: whole file)

    -s N            number of segments to skip (default: %d)

    -c nchan        set number of channels (default: %d)

    --size SIZE     number of bytes per I or Q variable (default: %d)
                    (i.e. 4 bytes or 1 byte)

    --seg_byte SEG_BYTE
                    sets how many bytes to read in between metadata blocks
                    (default: %d)

    -v              print out metadata for each block

''' % (pg, nskip, nchan, size_byte, seg_byte)


#-- parsing arguements --
if (len(inp) < 3):
    sys.exit(usage)

files = []
while (inp):
    k = inp.pop(0)
    if (k == '-n'):
        nseg = int(inp.pop(0))
    elif (k == '-s'):
        nskip = int(inp.pop(0))
    elif (k == '-r'):
        samprate = float(inp.pop(0))
        samprate *= 1e6
    elif (k == '-c'):
        nchan = int(inp.pop(0))
    elif (k == '-v'):
        verbose = True
    elif (k == '--fcen'):
        fcen = float(inp.pop(0))
    elif (k == '--seg_byte'):
        seg_byte = int(inp.pop(0))
        vector = False
    elif (k == '--size'):
        size_byte = int(inp.pop(0))
    elif (k == '--single'):
        dual = False
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        files.append(k)


freq = np.fft.fftfreq(nchan, d=1/samprate)
freq = np.fft.fftshift(freq)    # Hz
freq /= 1e6
freq += fcen    # MHz

augname = 'skip%d' % nskip
if (nseg == -1):
    augname += '_all'
else:
    augname += '_seg%d' % nseg


for filename in files:
    if not(os.path.isfile(filename)):
        print('error finding data: %s' % filename)
        next
    print('== plotting: %s\n' % filename)
    png = '%s.%s.view.png' % (filename, augname)

    ## loading data ##
    if (dual):
        nInp = 2
        tmp1, tmp2, timestamp = load_cli_data(filename, samprate, nseg=nseg, nskip=nskip, seg_byte=seg_byte, size_byte=size_byte, meta=meta, dual=dual)
        data = [tmp1, tmp2]
    else:
        nInp = 1
        tmp, timestamp = load_cli_data(filename, samprate, nseg=nseg, nskip=nskip, seg_byte=seg_byte, size_byte=size_byte, meta=meta, dual=dual)
        data = [tmp]


    if (False):  # debug
        print('timestamp.shape:', timestamp.shape)
        for i in range(len(data)):
            print(i, 'data.shape:', data[i].shape)
        sys.exit()


    ## plotting ##
    fig, sub = plt.subplots(4,nInp,figsize=(12,9),squeeze=False)
    for i in range(nInp):
        d1 = data[i].flatten()      # 1D voltage-vs-time
        nspec = len(d1) // nchan    # num of spec snapshots
        d2 = d1[:nspec*nchan].reshape((nspec, nchan))
        d2 = np.fft.fft(d2, axis=1)
        d2 = np.fft.fftshift(d2, axes=1)    # 2D spec

        #-- waterfall --
        ax = sub[0,i]
        s = ax.imshow(np.abs(d2.T), origin='lower', extent=[timestamp[0], timestamp[-1], freq[0], freq[-1]])
        ax.set_aspect('auto')
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('freq (MHz)')
        cb = plt.colorbar(s, ax=ax)
        cb.set_label('Amplitude')
        #ax.set_title('%s, %s, FFT' % (filename, augname))
        ax.set_title('input %d' % i)

        #-- spectrum --
        ax = sub[1,i]
        spec = np.abs(d2).mean(axis=0)
        ax.plot(freq, spec)
        ax.set_xlabel('freq (MHz)')
        ax.set_ylabel('amplitude')

        #-- t-series, zoom --
        ax = sub[2,i]
        nshow = 128
        ax.plot(timestamp[:nshow], d1[:nshow].real, label='real')
        ax.plot(timestamp[:nshow], d1[:nshow].imag, label='imag')
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('ADC')
        ax.legend()

        #-- t-series, full --
        ax = sub[3,i]
        nshow = -1
        ax.plot(timestamp[:nshow], d1[:nshow].real, label='real')
        ax.plot(timestamp[:nshow], d1[:nshow].imag, label='imag')
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('ADC')

    fig.tight_layout(rect=[0,0.03,1,0.95])
    fig.suptitle('%s' % (filename,))
    fig.savefig(png)
    plt.close(fig)

