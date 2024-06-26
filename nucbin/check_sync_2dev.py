#!/usr/bin/env python

from read_cli2 import *
import matplotlib.pyplot as plt
import time
import scipy.signal as sg
import gc
from subprocess import call


# data related
nInp    = 2         # dual channel
sampHz  = 10e6      # sample rate in Hz
dt      = 1/sampHz  # dt in sec
seg_byte  = 131072    # typical segment length in bytes
size_byte = 4       # sample size (I+Q) in bytes
tpSeg   = seg_byte/size_byte/nInp*dt   # segment lengthe in time

# analysis related
#tWin    = 0.512     # window size in time
tWin    = 1.024     # window size in time
sepWin  = tWin/2    # window separation
nseg    = np.ceil(tWin/tpSeg).astype(int)     # num of seg to read from each file, at least as long as tWin
nshort  = nseg      # cross-corr data length, also in num of seg
tRef    = 5.12      # reference window time, in sec
winlim  = [9,12]   # window range  # +/- 2.5s
#winlim  = [10,30]   # window range  # +/- 2.5s
#winlim  = [0, 40]   # about +/- 5s
thres   = 0.04       # a threshold for valid correlation peaks
odir    = '.'

## define pairs of data to correlate ##
pairs = [
    [0, 1],     # same board, bladeRF1
    [0, 2],     # cross-board
#    [1, 3]
    [1, 3],     # cross-board
    [2, 3]      # same board, bladeRF2
]
#pairs = [[0,1]] # for testing


inp = sys.argv[0:]
pg  = inp.pop(0)

usage   = '''
time-domain cross-correlation of 4 channels from 2 devices
syntax:
    %s <file1> <file2> [options]

options are:
    -o <DIR>        a directory to save outputs (default .)
    --ref Tsec      set reference time in sec
                    note: the search window is currently set to [2.56s:7.68s] 
    --win w1 w2     specify the search window range
                    each window is separated by 0.256s (half the window length)
''' % pg

if (len(inp) < 2):
    sys.exit(usage)

print('debug:', inp)

#files = ['cli1-crab-test-220119.bin', 'cli2-crab-test-220119.bin']
files = []
while (inp):
    k = inp.pop(0)
    if (k == '-o'):
        odir = inp.pop(0)
        call('mkdir -p %s'%odir, shell=True)
    elif (k == '--ref'):
        tRef = float(inp.pop(0))
    elif (k == '--win'):
        winlim[0] = int(inp.pop(0))
        winlim[1] = int(inp.pop(0))
    else:
        files.append(k)
        files.append(inp.pop(0))

for f in files:
    if (not os.path.isfile(f)):
        sys.exit('file not founc: %s' % f)
print('processing files:', files)


## loading reference data ##
refwin  = int(tRef/sepWin)  # reference window
nskip   = int(tRef/tpSeg)
data0, data1, timestamp1, acount, overrun = load_cli_data(files[0], sampHz, nskip=nskip, nseg=nseg, dual=True, verbose=0, meta=True)
tref = timestamp1[0]
refData = [data0.flatten(), data1.flatten()]

data0, data1, timestamp1, acount, overrun = load_cli_data(files[1], sampHz, nskip=nskip, nseg=nseg, dual=True, verbose=0, meta=True)
refData.append(data0.flatten())
refData.append(data1.flatten())
del data0, data1
gc.collect()


## begin looping over combinations ##
for loop in range(len(pairs)):
    refID, varID = pairs[loop]
    print('## ref input:', refID, '  test input:', varID, " ##\n")
    png = '%s/fft_sync2d_%d%d.png' % (odir, refID, varID)
    log = '%s/fft_sync2d_%d%d.log' % (odir, refID, varID)
    LOG = open(log, 'w')
    print('#window_i,  peak_amp, peak_off_cycle', file=LOG)

    samp0 = refData[refID]
    norm0 = np.sum(samp0*samp0.conjugate())
    nlen = samp0.size
    #print('samp0.shape', nlen)

    if (varID < 2):
        fname = files[0]
        inpID = varID
    else:
        fname = files[1]
        inpID = varID - 2


    ## setup figure ##
    fig, ax = plt.subplots(figsize=(10,6))
    toff = (timestamp1[-1]-timestamp1[0])/2.

    #cross2d = []
    time_y = []
    #for chunk_i in range(nchunk):
    for chunk_i in range(winlim[0],winlim[1]):
        idx0 = int(nshort/2 * chunk_i)
        #data0, data1, timestamp2 = load_cli_data(fname, sampHz, nskip=idx0, nseg=nseg, dual=True, verbose=0)
        data0, data1, timestamp2, acount2, overrun2 = load_cli_data(fname, sampHz, nskip=idx0, nseg=nseg, dual=True, meta=True, verbose=0)
        varData = [data0.flatten(), data1.flatten()]
        del data0, data1
        gc.collect()
        time_y.append(timestamp2[0])
        chunk = varData[inpID]
        norm1 = np.sum(chunk*chunk.conjugate())
        #print('chunk.shape', chunk.shape)

        t0 = time.time()
        cross = sg.correlate(chunk, samp0, mode='same')
        dsec = time.time() - t0
        #print('time used for corr: %.2f sec' % dsec)

        cross /= np.sqrt(norm0*norm1)
        amp = np.abs(cross)
        iPeak = np.argmax(amp)
        aPeak = np.max(amp)
        print('chunk_i, idx0, peak_arg:', chunk_i, idx0, iPeak)
        print('max(abs(cross)):', aPeak)
        if (aPeak > thres):
            print('%02d'%chunk_i, '%.3f'%aPeak, '%d'%(iPeak-nlen//2), file=LOG)
        #cross2d.append(cross)

        #t0 = time.time()
        #cross2 = sg.correlate(data0.flatten(), data1.flatten(), mode='same')
        #print('time used for 2nd corr:', time.time()-t0)
        #print('peak loc', np.argmax(np.abs(cross2)))

        if (chunk_i%2==0):
            off_y = 0.
        else:
            off_y = 0.2
        ax.plot(timestamp2-toff, amp*1.0+off_y)

        del cross
        gc.collect()
    #cross2d = np.array(cross2d)

    #s = ax.imshow(np.abs(cross2d), origin='lower', aspect='auto', extent=[timestamp1[0], timestamp1[-1], time_y[0], time_y[-1]])
    #cb = plt.colorbar(s, ax=ax)
    #ax.set_ylim([time_y[0], time_y[-1]+1.])
    #ax.set_yscale('log')
    ax.set_xlabel('window time (sec)')
    ax.set_ylabel('norm.corr (w. offset)')
    ax.set_ylim([-0.25, 1.45])
    ax.set_title('refCH:%d, refTime:%f; varCH:%d' % (refID, tref, varID))

    #ax.plot(time_y, time_y, 'k-', marker='o')

    fig.savefig(png)
    plt.close(fig)

    LOG.close()

