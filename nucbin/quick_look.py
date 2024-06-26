#!/usr/bin/env python

from read_cli2 import *
from loadh5 import *
import matplotlib.pyplot as plt
from subprocess import call
from astropy.stats import sigma_clip

rate0  = 20e6
nSamp = 16384000
#tskips = [0, 5, 10]
tskips = [0]
nchan = 1024

inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
extract metadata and save as .h5 if it does not exist
generate quick look plots for times specified in tskips (sec)
quick look includes:
    - time-avg spectrum
    - sample histogram

usage: %s <files>

note: files can be a list of either .log, .meta, or .h5 (or mix)

''' % (pg,)


if (len(inp)<1):
    sys.exit(usage)

files = []
while (inp):
    k = inp.pop(0)
    files.append(k)

files.sort()

# generate header files (skip existing ones)
call('cli2header.py %s' % (' '.join(files)), shell=True)

chan = np.arange(nchan)

for f in files:
    name = f.replace('.log','').replace('.meta','').replace('.h5','')
    fh5 = '%s.h5' % name
    odir = '%s.check' % name
    call('mkdir -p %s' % odir, shell=True)
    print('plotting', fh5)
    attrs = getAttrs(fh5)

    nInp = attrs.get('nInp')
    rate = attrs.get('rate') # list, Hz
    fcen = attrs.get('freq') # list, Hz
    gain = attrs.get('gain') # list, dB
    if (len(gain)==0):
        gain = list(gain)
        for ch in range(nInp):
            gain.append(0)

    #freq = np.fft.fftfreq(nchan, d=1/rate[0])  # in Hz
    #freq = np.fft.fftshift(freq)
    #fMHz = freq / 1e6   # in MHz

    for t in tskips:
        C0 = int(rate[0]*t)
        data, clock = load_cli_hdr(fh5, C0=C0, nSamp=nSamp)

        png = '%s/spec_hist_%.3fs.png' % (odir, t)
        f1, s1 = plt.subplots(nInp, 2, figsize=(10,4*nInp), squeeze=False)

        for ch in range(nInp):
            # histgram plot
            d = data[:,ch]
            ax = s1[ch, 1]
            tup = ax.hist(d.real[d.mask==False], bins=100, range=(-2048,2048), histtype='step', density=True, label='real')
            tup = ax.hist(d.imag[d.mask==False], bins=100, range=(-2048,2048), histtype='step', density=True, label='imag')
            f1024 = np.count_nonzero(np.abs(d.real)>1024) + np.count_nonzero(np.abs(d.imag)>1024)
            f1024 /= (len(d)*2)
            ax.text(0.05,0.90,'f(>1024)=%.1f%%'%(f1024*100), transform=ax.transAxes)
            ax.legend()
            ax.set_xlabel('I,Q (ADC)')
            ax.set_ylabel('Prob')
            ax.set_title('CH%d, gain=%.0fdB' % (ch,gain[ch]))

            # time-avg spectrum
            nlen = len(d)
            d = d[:nlen//nchan*nchan].reshape((-1,nchan))
            d = np.fft.fft(d, axis=1)
            d = np.fft.fftshift(d, axes=1)
            spec = d.mean(axis=0)
            ampld = np.ma.abs(d).mean(axis=0)
            clip_d = sigma_clip(np.abs(d), axis=0)
            clip_amp = clip_d.mean(axis=0)
            ax = s1[ch, 0]
            ax.plot(chan, ampld, label='mean(abs(raw))')
            ax.plot(chan, clip_amp, label='mean(abs(clip))')
            ax.plot(chan, np.abs(spec), label='abs(mean(raw))')
            #ax.set_yscale('log')
            ax.legend()
            ax.set_xlabel('channel')
            ax.set_ylabel('abs(spec)')
            ax.set_title('CH%d, fcen=%.0fMHz, bw=%.0fMHz' % (ch, fcen[ch]/1e6, rate[ch]/1e6))

        if (nInp > 1):
            s1[0,0].get_shared_y_axes().join(s1[0,0], s1[1,0])

        f1.tight_layout(rect=[0,0.03,1,0.95])
        f1.suptitle('%s, t=%.3fs' % (fh5, t))
        f1.savefig(png)
        plt.close(f1)


