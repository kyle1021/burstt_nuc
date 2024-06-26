#!/usr/bin/env python

from loadh5 import *
from read_cli2 import *
import matplotlib.pyplot as plt
import gc
from astropy.time import Time

home = os.getenv('HOME')


def toPower(x, mode='vol'):
    '''
    convert linear voltage to power in dB
    if mode=='vis', convert visibility to power in dB
    '''
    if (mode=='vol'):
        p = 20.*np.log10(x)
    elif (mode=='vis'):
        p = 10.*np.log10(x)
    return p


nWin   = 1      # number of windows to analysis
winDur = 1.     # window duration in sec
winSep = 10     # window separation in sec
winOff = 0.     # overall window offset in sec.
fMHz   = 0.
srcFlux = 0.
do_append = False

nChan  = 1024   # default num of spec channel
nInp   = 2      # default num of input

inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
process the raw data at selected time windows and produce .vish5
this version takes 1 window from each file
and files can be filtered by central freq
add variance in addition to time-avg results

syntax:
    %s <dirs> [options]

    <dirs> are folders that each contains one loop of multiple frequency files

options are:
    -n nWin     # specify the number of windows to process per file (%d)
    -d winDur   # specify the window duration in seconds (%.1f)
    -s winSep   # specify the separation between windows in seconds (%.1f)
    -o winOff   # specify the overall window offset in seconds (%.1f)
    -c nChan    # specify the number of spectral channels (%d)
    -f 'freq1 [freq2 ...]'
                # filter only files that match the freq in MHz (%.0f)
                # specifying multiple freq (quoted) will repeat the analysis at 
                # different freq

    --append    # append new window to existing arrays

''' % (pg, nWin, winDur, winSep, winOff, nChan, fMHz)

if (len(inp) < 1):
    sys.exit(usage)


dirs = []
fMHzs = []
while (inp):
    k = inp.pop(0)
    if (k == '-n'):
        nWin = int(inp.pop(0))
    elif (k == '-d'):
        winDur = float(inp.pop(0))
    elif (k == '-s'):
        winSep = float(inp.pop(0))
    elif (k == '-c'):
        nChan = int(inp.pop(0))
    elif (k == '-o'):
        winOff = float(inp.pop(0))
    elif (k == '-f'):
        tmp = inp.pop(0).split()
        for tt in tmp:
            fMHzs.append(float(tt))
    elif (k == '--append'):
        do_append = True
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        dirs.append(k)


if (len(fMHzs) == 0):
    sys.exit("please specify a frequency to process with '-f'")

print('dirs to process:', dirs)
nDir = len(dirs)

for fMHz in fMHzs:
    bmark1 = time.time()
    print('central freq: %.0fMHz' % fMHz)

    files = []
    gains = []
    rates = []
    fcens = []
    for d in dirs:
        tmp = glob('%s/*.h5' % d)
        if (len(tmp) == 0):
            print('skip folder without .h5 files:', d)
            continue
        for f in tmp:
            if (f.endswith('.spec.h5')):
                continue
            attrs = getAttrs(f)
            fcen = attrs['freq'][0]
            gain = attrs['gain'][0]
            rate = attrs['rate'][0]
            if (np.isclose(fcen/1e6, fMHz)):
                files.append(f)
                gains.append(gain)
                rates.append(rate)
                fcens.append(fcen)

    nFile = len(files)
    if (nFile == 0):
        sys.exit('no valid files for fcen=%.0f' % fcen0)
    #print('fcen=%.0fMHz, files:', files)
    print('nFile =', nFile)
    fcen = np.median(fcens)
    gain = np.median(gains)
    rate = np.median(rate)

    fname = 'analyze_%.0f' % fMHz
    fout = fname + '.vish5'
    odir = fname + '.check'
    if (not os.path.isdir(odir)):
        call('mkdir %s'%odir, shell=True)

    # bandpass
    fbp = '%s/local/bin/blade_bandpass.h5' % home   # the freq-independent bandpass file, for nChan=1024 / 40MHz BW
    medbp = getData(fbp, 'median')          # bandpass in dB
    #linbp = 10**(medbp/10.)                 # in linear visbility scale
    gref = 20                               # the reference gain when bandpass was measured

    bp = medbp + gain - gref    # bandpass under the current gain, in dB
    linbp = 10**(bp/20)         # linear bandpass for voltage

    BWch = rate / nChan     # channel bandwidth in Hz
    Tint = 1./rate * nChan  # instantaneous integration time in sec

    fh5 = files[0]
    idx = fh5.find('_dev')
    str_open = fh5[idx-13:idx]
    dt0 = datetime.strptime(str_open, '%y%m%d_%H%M%S')
    at0 = Time(dt0, format='datetime')  # local time
    ut0 = at0 - 8/24    # convert to UTC
    unix0 = ut0.to_value('unix')


    # update last attrs
    attrs['rate'] = [rate, rate]
    attrs['gain'] = [gain, gain]
    attrs['freq'] = [fcen, fcen]
    # additional info to write
    attrs['winDur'] = winDur
    attrs['nChan'] = nChan
    attrs['bpfile'] = fbp
    attrs['unix_utc_open'] = unix0
    putAttrs(fout, attrs)


    IFreq = np.fft.fftfreq(nChan, d=1/rate)
    IFreq = np.fft.fftshift(IFreq)  # IF in Hz
    RFreq = IFreq + fcen            # RF in Hz
    freq  = RFreq / 1e6             # RF in MHz

    savdata = {}
    savmask = {}
    keys = ['winSec', 'winNFT', 'winSpec', 'winVar', 'winCoeff', 'winMedCoeff']
    print('initializing:', keys)
    redo = True
    if (os.path.isfile(fout) and do_append):
        redo = False
        for k in keys:
            tmp = getData(fout, k)
            if (tmp is None):
                print('... %s not found in %s'%(k, fout))
                redo = True     # any key not found will result in redo
                break
            else:
                savdata[k] = tmp.tolist()
                if (isinstance(tmp, np.ma.MaskedArray)):
                    savmask[k] = tmp.mask.tolist()
                else:
                    savmask[k] = None

    if (redo):
        for k in keys:
            savdata[k] = []
            savmask[k] = []
            if (k == 'winSec'):
                savmask[k] = None



    # processing files
    for i in range(nFile):
        fh5 = files[i]
        print('processing:', fh5, '(%d/%d)'%(i+1,nFile))

        idx = fh5.find('_dev')
        str_open = fh5[idx-13:idx]
        dt_open = datetime.strptime(str_open, '%y%m%d_%H%M%S')
        off_open = (dt_open - dt0).total_seconds()  # offset in file opening time, in seconds

        ## processing windows
        for j in range(nWin):
            print('window:(%d/%d)'%(j+1, nWin))
            tWin = winOff + j*winSep

            sav_tWin = tWin + off_open
            if (sav_tWin in savdata['winSec']):
                print('skip existing window:', sav_tWin, fh5)
                continue

            savdata['winSec'].append(sav_tWin)
            C0 = int(tWin * rate)
            nSamp = int(winDur * rate)
            nFrame = nSamp // nChan
            nSamp = nFrame * nChan    # regularize the nSamp to multiples of nChan

            tag = 't%08.1f' % sav_tWin

            data, clock = load_cli_hdr(fh5, C0=C0, nSamp=nSamp, verbose=0)

            FTdata = np.ma.array(np.zeros((nInp, nFrame, nChan), dtype=complex), mask=False)
            for k in range(nInp):
                FTdata[k] = maskedFFT(data[:,k])
                #print('main: (nmask)', np.count_nonzero(FTdata.mask[k]))

            nFTdata = FTdata / linbp.reshape((1,1,nChan))   # convert to dBm/ch
            #nFTspec = np.ma.abs(nFTdata).mean(axis=1)          # shape (nInp, nChan)
            nFTspec = np.ma.median(np.ma.abs(nFTdata), axis=1)          # shape (nInp, nChan)
            savdata['winNFT'].append(nFTspec)
            savmask['winNFT'].append(nFTspec.mask)



            ## vis vs. freq plot

            png1 = '%s/single_vis.%s.png' % (odir, tag)
            f1, s1 = plt.subplots(2,2,figsize=(12,9), sharex=True)

            for k in range(nInp):
                ax = s1[0,k]
                #ax.plot(freq, toPower(np.ma.abs(nFTdata[k]).mean(axis=0)))
                ax.plot(freq, toPower(nFTspec[k]))
                ax.set_ylabel('power (dBm)')
                ax.set_title('ch%d'%k)

            xcorr = nFTdata[0] * nFTdata[1].conjugate()
            xsr   = np.ma.median(xcorr.real, axis=0)
            xsi   = np.ma.median(xcorr.imag, axis=0)
            xspec = xsr + 1.j*xsi
            savdata['winSpec'].append(xspec)
            savmask['winSpec'].append(xspec.mask)


            ax = s1[1,0]
            ax.plot(freq, toPower(np.ma.abs(xspec),mode='vis'))
            ax.set_ylabel('vis.power (dBm)')
            ax.set_xlabel('Freq (MHz)')
            ax.set_title('corr01')

            ax = s1[1,1]
            ax.plot(freq, np.ma.angle(xspec))
            ax.set_ylabel('vis.phase (rad)')
            ax.set_xlabel('Freq (MHz)')
            ax.set_title('corr01')

            s1[0,0].get_shared_y_axes().join(s1[0,0], s1[0,1], s1[1,0])

            f1.tight_layout(rect=[0,0.03,1,0.95])
            f1.suptitle('%s, %s'%(fh5, tag))
            f1.savefig(png1)
            plt.close(f1)


            ## variance vs freq plot

            png2 = '%s/single_var.%s.png' % (odir, tag)
            f2, s2 = plt.subplots(1,1)

            if (srcFlux > 0):
                calFac = np.ma.abs(xspec) / srcFlux                  # TargetSpec / calFac --> targetFlux
            else:
                calFac = 1.
            var = xcorr.var(axis=0)
            vmask = xcorr.mask.all(axis=0)
            var.mask = vmask
            print(var, vmask)
            print(var.shape, vmask.shape, xcorr.shape)
            print(np.count_nonzero(xcorr.mask, axis=0))
            var /= (calFac*calFac)                # Jy^2, or mW^2 if calFac==1
            savdata['winVar'].append(var)
            savmask['winVar'].append(vmask)

            ax = s2
            ax.plot(freq, var)
            ax.set_ylabel('variance (mW^2)')
            ax.set_yscale('log')
            ax.set_xlabel('Freq (MHz)')
            ax.set_title('%s, %s'%(fh5, tag))
            f2.tight_layout()
            f2.savefig(png2)
            plt.close(f2)



            ## corr-coefficient plot

            png3 = '%s/single_corr-coeff.%s.png' % (odir, tag)
            f3, s3 = plt.subplots(1,1)

            coeff = xcorr / (np.ma.abs(nFTdata[0]) * np.ma.abs(nFTdata[1]))
            #coeff2 = xcorr / (nFTspec[0] * nFTspec[1])
            tmpr = np.ma.median(coeff.real, axis=0)
            tmpi = np.ma.median(coeff.imag, axis=0)
            coeffspec = tmpr + 1j*tmpi
            #tmpr = np.ma.median(coeff2.real, axis=0)
            #tmpi = np.ma.median(coeff2.imag, axis=0)
            #coeffspec2 = tmpr + 1j*tmpi

            #coeff.fill_value = 0.
            #par = np.percentile(np.ma.abs(coeff.filled()), [0,15,50,85,100], axis=0)
            #med_coeff = np.median(par[2])
            med_coeff  = np.ma.median(np.ma.abs(coeffspec))
            #med_coeff2 = np.ma.median(np.ma.abs(coeffspec2))
            savdata['winCoeff'].append(coeffspec)
            savmask['winCoeff'].append(coeffspec.mask)
            savdata['winMedCoeff'].append(med_coeff)
            #savmask['winMedCoeff'].append(med_coeff.mask)

            ax = s3
            #ax.fill_between(freq, par[0], par[4], color='b', alpha=0.2, label='min/max')
            #ax.fill_between(freq, par[1], par[3], color='b', alpha=0.6, label='15%/85% pct')
            #ax.plot(freq, par[2], color='b', label='med')
            ax.plot(freq, np.ma.abs(coeffspec), color='b', label='method1')
            #ax.plot(freq, np.ma.abs(coeffspec2), color='g', label='method2')
            ax.text(0.05, 0.90, 'spec.med =%.3g'%med_coeff, transform=ax.transAxes)
            #ax.text(0.05, 0.80, 'spec.med2=%.3g'%med_coeff2, transform=ax.transAxes)
            #ax.legend()
            ax.set_ylabel('corr-coeff')
            ax.set_yscale('log')
            ax.set_xlabel('Freq (MHz)')
            ax.set_title('%s, %s'%(fh5, tag))
            f3.tight_layout()
            f3.savefig(png3)
            plt.close(f3)


            del data, clock, FTdata, nFTdata
            del xcorr, xsr, xsi, xspec, coeff
            gc.collect()


    winSec = np.array(savdata['winSec'])
    adoneh5(fout, winSec, 'winSec')
    putAttrs(fout, {'note':'window time in seconds'}, dest='winSec')

    winNFT = np.ma.array(savdata['winNFT'], mask=savmask['winNFT'])
    adoneh5(fout, winNFT, 'winNFT')
    #putAttrs(fout, {'note':'window voltage data after bandpass calibration'}, dest='winNFT')
    putAttrs(fout, {'note':'window time-avg antenna abs after bandpass calibration'}, dest='winNFT')

    winSpec = np.ma.array(savdata['winSpec'], mask=savmask['winSpec'])
    adoneh5(fout, winSpec, 'winSpec')
    putAttrs(fout, {'note':'window time-avg of the xcorr spectrum'}, dest='winSpec')

    adoneh5(fout, freq, 'freq')
    putAttrs(fout, {'note':'RF in MHz'}, dest='freq')

    winVar = np.ma.array(savdata['winVar'], mask=savmask['winVar'])
    adoneh5(fout, winVar, 'winVar')
    putAttrs(fout, {'note':'window time-avriance of the xcorr spectrum'}, dest='winVar')

    winCoeff = np.ma.array(savdata['winCoeff'], mask=savmask['winCoeff'])
    adoneh5(fout, winCoeff, 'winCoeff')
    putAttrs(fout, {'note':'window correlation coefficients taking the r-/i- median in time',
                    'median':savdata['winMedCoeff']}, dest='winCoeff')

    bmark2 = time.time()
    print('process %d windows in: %f sec' % (nFile, bmark2-bmark1))

