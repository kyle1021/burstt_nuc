#!/usr/bin/env python

from read_cli2 import *
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from subprocess import call


samprate = 20e6
fcen_MHz = 560
seg_byte = 131072
size_byte = 2   # 2Bytes per I or Q
nAnt = 2        # 2 for dual channel
dur = 2.0       # waterfall duration in sec
tplots = [0.]   # plot times in sec
nchan = 1024
nchan2 = 256
#chlim = [1600, 2000]    # zoom in window corresponding to 8192 chan
odir = 'water_view'

inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
plot the normalized waterfall spectrum
i.e. the 2D spectrum is normalized by the time-avg spectrum and subtract by 1
this method enhances the RFI in the 'signal'-free region (e.g. regions next to a strong TV band)

note: this version assumes metadata is included in the file

syntax:
    %s <files> [options]

    multiple files can be specified, provided they have the same rate, etc.

options are
-r RATE         set the sampling rate (in MHz) of the file (default: %f MHz)
-f FCEN         set the central frequency in MHz (only for plot labels)
-o DIR          save the output to another directory (default: %s)
-d DUR          set the duration (in sec) of the waterfall (default: %.3fs)
-t 't_sec1 t_sec2 ...'
                make multiple plots at the specified times (in sec)
                note: quote multiple numbers and use space to separate
                (default: [0])

''' % (pg, samprate/1e6, odir, dur)

if (len(inp) < 1):
    sys.exit(usage)

files = []
while (inp):
    k = inp.pop(0)
    if (k == '-r'):
        rate_MHz = float(inp.pop(0))
        samprate = rate_MHz * 1e6
    elif (k == '-f'):
        fcen_MHz = float(inp.pop(0))
    elif (k == '-o'):
        odir = inp.pop(0)
    elif (k == '-d'):
        dur = float(inp.pop(0))
    elif (k == '-t'):
        strplots = inp.pop(0)
        tplots = [float(x) for x in strplots.split()]
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        files.append(k)

call('mkdir -p %s' % odir, shell=True)

tperseg = seg_byte / (size_byte*2*nAnt) / samprate  # segment duration in sec
nseg = int(dur//tperseg)     # rounded number of seg to reach dur

for rpt in range(len(tplots)):
    nskip = int(tplots[rpt] / tperseg)  # rounded number of seg to skip
    print('plotting at: %fs' % tplots[rpt])
    print(' nseg, nskip:', nseg, nskip)

    rawdata = []
    ftdata = []     # for data with nchan
    rawtime = []
    alltime = []
    for fname in files:
        print('loading data:', fname)
        #data0, data1, time1 = load_cli_data(fname, samprate, nskip=nskip, nseg=nseg, seg_byte=seg_byte, nchan=nchan, dual=True, verbose=0)
        data0, data1, time1, acount1, overrun1 = load_cli_data(fname, samprate, nskip=nskip, nseg=nseg, seg_byte=seg_byte, nchan=1, dual=True, verbose=0, meta=True)
        tbegin = time1[0]
        time1 -= tbegin
        print('tbegin:', tbegin)
        ft0 = np.fft.fft(data0.reshape((-1,nchan)), axis=1)
        ft0 = np.fft.fftshift(ft0, axes=1)
        ft1 = np.fft.fft(data1.reshape((-1,nchan)), axis=1)
        ft1 = np.fft.fftshift(ft1, axes=1)
        tspec1 = time1[0::nchan]
        rawdata.extend([data0, data1])
        ftdata.extend([ft0, ft1])
        rawtime.extend([time1, time1])
        alltime.extend([tspec1, tspec1])

    nInp = len(ftdata)
    freq = np.fft.fftfreq(nchan, d=1/samprate)
    freq = np.fft.fftshift(freq)    # Hz
    freq /= 1e6                     # MHz
    ch   = np.arange(nchan)

    ftdata2 = []    # for data with nchan2
    alltime2 = []
    for chi in range(nInp):
        data0 = rawdata[chi]
        data0 = data0.reshape((-1,nchan2))
        ft0 = np.fft.fft(data0, axis=1)
        ft0 = np.fft.fftshift(ft0, axes=1)
        ftdata2.append(ft0)
        #time2 = alltime[chi][0] + np.arange(data0.shape[0])*nchan2/samprate
        time2 = rawtime[chi][0::nchan2]
        alltime2.append(time2)
    ch2 = np.arange(nchan2)
    freq2 = np.fft.fftfreq(nchan2, d=1/samprate)
    freq2 = np.fft.fftshift(freq2)    # Hz
    freq2 /= 1e6                     # MHz


    png1 = '%s/waterfall_%.4fs.png' % (odir, tbegin)
    print('saving ', png1)
    f1, sub1 = plt.subplots(nInp, 4, figsize=(16,2.5*nInp))
    for chi in range(nInp):
        spec = np.abs(ftdata[chi]).mean(axis=0)
        ax = sub1[chi,0]
        ax.plot(ch, spec)
        ax.set_xlabel('chan')
        ax.set_ylabel('abs(spec).mean')
        ax.set_title('ch%d spec %dMHz' % (chi, samprate/1e6))

        nspec = np.abs(ftdata[chi]) / spec.reshape((1,-1))
        nspec -= 1.
        ax = sub1[chi,1]
        rnspec = nspec[:(nspec.shape[0]//200)*200,:].reshape((-1,200,nchan)).mean(axis=1)
        #s = ax.imshow(rnspec.T, origin='lower', norm=colors.SymLogNorm(linthresh=1e-3, linscale=1e-3, vmin=-1, vmax=1), extent=[alltime[chi][0], alltime[chi][-1], ch[0], ch[-1]], aspect='auto')
        rtime = alltime[chi][0::200]
        rtime = rtime[:rnspec.shape[0]]
        #print('debug:', rtime.shape, ch.shape, rnspec.T.shape)
        s = ax.pcolormesh(rtime, ch, rnspec.T, norm=colors.SymLogNorm(linthresh=1e-3, linscale=1e-3, vmin=-1, vmax=1))
        cb = plt.colorbar(s, ax=ax)
        ax.set_xlabel('time (sec)')
        ax.set_ylabel('chan')
        tres1 = nchan*200/samprate*1e3  # time resolution 1 in ms
        tspn1 = dur     # time span 1 in s
        ax.set_title('ch%d frac.spec, 1024ch, %dsec@%.2fms' % (chi, tspn1, tres1))

        spec2 = np.abs(ftdata2[chi]).mean(axis=0)
        nspec2 = np.abs(ftdata2[chi]) / spec2.reshape((1,-1))
        nspec2 -= 1.
        ax = sub1[chi,2]
        #s = ax.imshow(nspec2[:80,:].T, origin='lower', norm=colors.SymLogNorm(linthresh=1e-3, linscale=1e-3, vmin=-10, vmax=10), extent=[alltime2[chi][0]*1000, alltime2[chi][80]*1000, ch2[0], ch2[-1]], aspect='auto')
        s = ax.pcolormesh(alltime2[chi][:80], ch2, nspec2[:80,:].T, norm=colors.SymLogNorm(linthresh=1e-3, linscale=1e-3, vmin=-1, vmax=1))
        tres2 = nchan2/samprate*1e6     # time resolution 2 in mus
        tspn2 = nchan2*80/samprate*1e3  # time span 2 in ms
        cb = plt.colorbar(s, ax=ax)
        ax.set_xlabel('time (ms)')
        ax.set_ylabel('chan')
        ax.set_title('ch%d frac.spec, 256ch, %dms@%.1fmus' % (chi, tspn2, tres2))

        ax = sub1[chi,3]
        ax.plot(freq, np.std(rnspec, axis=0), label='std@%.2fms' % tres1)
        ax.plot(freq2, np.std(nspec2[:80,:], axis=0), label='std@%.1fmus' % tres2)
        #ax.plot(ch, np.mean(nspec, axis=0)+0.5, label='mean+0.5')
        ax.legend()
        ax.set_xlabel('freq (MHz)')
        ax.set_ylabel('std(frac.spec)')
        ax.set_title('ch%d std' % chi)

    sub1[0,0].get_shared_y_axes().join(*sub1[:,0])
    sub1[0,3].get_shared_y_axes().join(*sub1[:,3])
    f1.tight_layout(rect=[0,0.03,1,0.95])
    #f1.suptitle('%s, %s, 550--570MHz@20Msps, tbegin=%fs' % (files[0], files[1], tbegin))
    f1.suptitle('%s, fcen=%dMHz, %dMsps, tbegin=%fs' % (files[0], fcen_MHz, samprate/1e6, tbegin))
    f1.savefig(png1)
    plt.close(f1)


