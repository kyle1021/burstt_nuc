#!/usr/bin/env python
import sys, os.path, gc
from read_cli2 import *
from scipy.signal import correlate
import matplotlib.pyplot as plt

#files = ['sun_220212_071100_dev1.meta', 'sun_220212_071100_dev2.meta']
samprate = 10e6

nchunk = 1000   # num of seg in a chunk, 1seg=16384samples
mlim = [0,3]    # check the first 10 chunks
nlen = 16384*nchunk
odir = '.'

pairs = [
            [0,1],
            [0,2],
            [1,3],
            [2,3]
        ]


inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
check data synch by correlating chunks of inputs multiple times
chunk size is 1000 segs = 16384*1000 samples

** note: both device files should be passed in the command line

syntax:
%s <files> [options]

options are:
-o ODIR         # output dir to save the plot

''' % (pg, )

if (len(inp) < 2):
    sys.exit(usage)


files = []
while (inp):
    k = inp.pop(0)
    if (k == '-o'):
        odir = inp.pop(0)
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        files.append(k)


all_corr = []       # length = num of chunks
for m in range(mlim[0], mlim[1]):
    print('chunk', m)
    nskip = m * nchunk

    data = []       # length = num of inputs
    norm = []
    for f in files:
        d0, d1, ts, ac, ov = load_cli_data(f, samprate, nseg=nchunk, nskip=nskip, meta=True, dual=True, verbose=0)
        d0 = d0.flatten()
        d1 = d1.flatten()
        data.extend([d0, d1])
        norm.extend([np.sum(d0*d0.conjugate()), np.sum(d1*d1.conjugate())])
        del d0, d1
        gc.collect()

    chunk_corr = [] # length = num of pairs
    for i, j in pairs:
        corr = correlate(data[i], data[j], mode='same') # length = nchunk
        corr /= np.sqrt(norm[i]*norm[j])
        corr = np.abs(corr)
        pk  = corr.max()
        pki = corr.argmax() - nlen/2
        print('peak%d-%d, max, off: %f, %d' % (i,j,pk,pki))
        chunk_corr.append(corr)

    all_corr.append(chunk_corr)

all_corr = np.array(all_corr)       # shape = (nmlim, npair, nlen)
avg_corr = all_corr.mean(axis=0)    # shape = (npair, nlen)

f1, sub1 = plt.subplots(2,2,sharex=True,sharey=True,figsize=(8,6))
s1 = sub1.flatten()
x = np.arange(nlen) - nlen/2
p = -1
for i,j in pairs:
    p += 1
    ax = s1[p]
    pname = '%d-%d' % (i,j)
    for m in range(all_corr.shape[0]):
        cc = 'C%d'%m
        pki = all_corr[m,p].argmax() - nlen/2
        ax.plot(x, all_corr[m,p], color=cc, label='%d'%pki)
    pki = avg_corr[p].argmax() - nlen/2
    ax.plot(x, avg_corr[p], color='k', label='mean=%d'%pki)
    ax.legend()
    ax.set_xlabel('offset (samples)')
    ax.set_ylabel('norm.corr')
    ax.set_title('inputs: %s' % pname)

p1 = '%s/mean_sync.png' % odir
f1.savefig(p1)
plt.close(f1)


