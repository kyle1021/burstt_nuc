#!/usr/bin/env python

from read_cli2 import *
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from glob import glob
from loadh5 import *
from astropy.time import Time
import re

home = os.getenv('HOME')

nchan0 = 1024   # for the reference
nchan = nchan0  # for 40M BW
fbp = '%s/local/bin/blade_bandpass.h5' % home


inp = sys.argv[0:]
pg  = inp.pop(0)

name0 = 'data'
dev   = '*'
yset  = False
fref  = ''
bpname= 'median'
bpoff = 0.

usage = '''
plots spectrum of all files in a folder (assumed with different central freq)
onto the same SED plot
apply a median bandpass correction

note: this version does not require the .h5 headers

sytax: %s <DIR(s)> [options]

options are:
--name NAME         beginning part of the filename (%s)
--dev DEV           specify one device id (e.g. dev1, dev3, ...)
                    (default: %s) 
--ylim ymin ymax    in dB (default: auto)
--ref <SED....spec.h5>  include a reference spectrum in the plot
--bw20              specify the data has only 20MHz BW
                    to plot with 512 channels

''' % (pg, name0, dev)


if (len(inp) < 1):
    sys.exit(usage)

dirs = []
while (inp):
    k = inp.pop(0)
    if (k == '--name'):
        name0 = inp.pop(0)
    elif (k == '--dev'):
        dev = inp.pop(0)
    elif (k == '--ylim'):
        ymin = float(inp.pop(0))
        ymax = float(inp.pop(0))
        yset = True
    elif (k == '--ref'):
        fref = inp.pop(0)
    elif (k == '--bw20'):
        bpname = 'median_20M'
        nchan = 512
        bpoff = -20.
    elif (k.startswith('-')):
        sys.exit('unknown option: %s'%k)
    else:
        dirs.append(k)
dirs.sort()
print('to process:', dirs)


med_bp = getData(fbp, bpname) # median bandpass in dB (1024 ch), with 20dB gain
med_bp += bpoff

if (fref != ''):
    if ('blade_terminate.h5' in fref):
        spec = 'LNAon_40dB'
        refspec = getData(fref, spec)   # shape (ntune, nchan)
    else:
        spec = 'specdB'
        refspec = getData(fref, spec)   # shape (ntune, nInp, nchan)
        refspec = refspec.mean(axis=1)  # avg of dual-input and reduce index (ntune, nchan)
    reffreq = getData(fref, 'freqMHz')  # shape (ntune, nchan)


if (dev == '*'):
    trail = 'all'
else:
    trail = dev


for idir in dirs:
    head = os.path.basename(idir)
    png1 = '%s/SED_%s_%s.png' % (idir, head, trail)
    fspec = png1.replace('.png', '.spec.h5')

    f1, s1 = plt.subplots(2,1,figsize=(15,9), sharey=True, sharex=True)

    for run in range(1):
        yoff = 0    # figure trace offset in y
        g = 0       # extra gain adjust in dB
        edge = 4    # num of edge channels to mask

        files = glob('%s/%s*_%s.meta'%(idir, name0, dev))
        print('%s/%s*_%s.meta'%(idir, name0, dev))
        print("files = ", files)
        files.sort()
        nfiles = len(files)

        allfreqs = []
        allspecs = []
        alltimes = []
        devs = []
        #cc = 'C%d'%run
        for fi, f in enumerate(files):
            print('loading: %s, (%d/%d)' % (f, fi+1, nfiles))
            cc = 'C%d' % (fi%10)

            print('meta size:', os.path.getsize(f))
            if (os.path.getsize(f) < 16384):
                continue

            flog = f.replace('.meta', '.log')
            ii = f.find('_dev')
            #idev = f[ii+1:ii+5]
            idev = re.search("_(dev\d*)\.", flog).groups()[0]
            devs.append(idev)

            idtstr = f[ii-13:ii]
            alltimes.append(datetime.strptime(idtstr, '%y%m%d_%H%M%S'))

            #attrs = lh.getAttrs(f)
            attrs = logParse(flog)
            fcen = attrs['freq'][0]/1e6 # to MHz
            rate = attrs['rate'][0]/1e6 # to MHz
            gain_log = attrs['gain']    # gain recorded in log files, a list of nInp elements, in dB
            freq = np.fft.fftfreq(nchan, d=1/rate)
            freq = np.fft.fftshift(freq)
            freq += fcen
            allfreqs.append(freq)

            #data, clock = load_cli_hdr(f, C0=0, nSamp=16384000, verbose=0)
            data, clock = load_cli_data2(f, nseg=1000, nskip=0, verbose=0)
            filespec = []
            for ch in range(2):
                g = gain_log[ch] - 20   # g = extra gain above 20dB, adopted in the med_bp
                ftdata = maskedFFT(data[:,ch], nchan)
                spec = np.ma.abs(ftdata).mean(axis=0)
                spec = 20.*np.ma.log10(spec)
                spec -= g
                spec -= med_bp
                ax = s1[ch]
                ax.grid(1)
                if (fi == 0):
                    ax.plot(freq[edge:nchan-edge], spec[edge:nchan-edge]+yoff, color=cc, label='0.0sec')
                else:
                    ax.plot(freq[edge:nchan-edge], spec[edge:nchan-edge]+yoff, color=cc)
                filespec.append(spec)
            allspecs.append(filespec)

    if (fref != ''):
        ntune = reffreq.shape[0]
        for fi in range(ntune):
            for ch in range(2):
                ax = s1[ch]
                if (fi == 0):
                    #ax.plot(reffreq[fi,edge:nchan-edge], refspec[fi,ch,edge:nchan-edge]+yoff, color='gray', label='LNAterm')
                    ax.plot(reffreq[fi,edge:nchan0-edge], refspec[fi,edge:nchan0-edge]+yoff, color='gray', label='LNAterm')
                else:
                    ax.plot(reffreq[fi,edge:nchan0-edge], refspec[fi,edge:nchan0-edge]+yoff, color='gray')

    for ch in range(2):
        ax = s1[ch]
        ax.set_xlabel('freq (MHz)')
        ax.set_ylabel('ampld (dB)')
        ax.set_title('ch%d'%ch)
        ax.legend()
        if (yset):
            ax.set_ylim([ymin, ymax])

    udev = np.unique(devs)
    sptitle = '%s,' % idir
    for ud in udev:
        sptitle += ' %s'%ud

    f1.tight_layout(rect=[0,0.03,1,0.95])
    f1.suptitle(sptitle)
    f1.savefig(png1)
    plt.close(f1)

    adoneh5(fspec, np.array(allfreqs), 'freqMHz')
    adoneh5(fspec, np.array(allspecs), 'specdB')

    allepoch = Time(alltimes, format='datetime').to_value('unix')
    adoneh5(fspec, np.array(allepoch), 'unixtime')

