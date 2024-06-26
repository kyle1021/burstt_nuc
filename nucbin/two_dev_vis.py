#!/usr/bin/env python

from loadh5 import *
from read_cli2 import *
import matplotlib.pyplot as plt
import gc
from astropy.time import Time
from scipy.signal import correlate
from astropy.stats import sigma_clip

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
devs   = []     # no defaults
do_append = False
auto_align = False
med_off0 = 0
tChunk = 1.0    # split the window into chunks

nChan  = 1024   # default num of spec channel
nInp   = 4      # default num of input
nBl    = int(nInp * (nInp-1) / 2)

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
                # NOTE: if winDur > tChunk(%.1fsec), then the integration is split into chunks
                # the chunks are then averaged before plotting and saving to .vish5
    -s winSep   # specify the separation between windows in seconds (%.1f)
    -o winOff   # specify the overall window offset in seconds (%.1f)
    -c nChan    # specify the number of spectral channels (%d)
    -f 'freq1 [freq2 ...]'
                # filter only files that match the freq in MHz (%.0f)
                # specifying multiple freq (quoted) will repeat the analysis at 
                # different freq

    --dev 'devN devM'
                # the device id of two units
                # e.g. dev1 dev2
                # (data are assumed to be stored in subfolders as the device id)

    --append    # append new window to existing arrays

    --align     # auto-align cross-dev offset

    --fixed OFFSET
                # specify a fixed offset between the two devices, in samples

''' % (pg, nWin, winDur, tChunk, winSep, winOff, nChan, fMHz)

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
    elif (k == '--dev'):
        tmp = inp.pop(0)
        devs = tmp.split()
    elif (k == '--append'):
        do_append = True
    elif (k == '--align'):
        auto_align = True
    elif (k == '--fixed'):
        med_off0 = int(inp.pop(0))
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        dirs.append(k)


if (len(fMHzs) == 0):
    sys.exit("please specify a frequency to process with '-f'")
if (len(devs) == 0):
    sys.exit("please specify two device ids with '-d'")


#print('dirs to process:', dirs)
#nDir = len(dirs)
dev1, dev2 = devs
dirs_byDev = {dev1:[], dev2:[]}
for d in dirs:
    if (d.startswith(dev1)):
        dirs_byDev[dev1].append(d)
    elif (d.startswith(dev2)):
        dirs_byDev[dev2].append(d)
print('num of dirs: %s/%d, %s/%d' % (dev1, len(dirs_byDev[dev1]), dev2, len(dirs_byDev[dev2])))
print(dirs_byDev[dev1], dirs_byDev[dev2])


for fMHz in fMHzs:
    bmark1 = time.time()
    print('central freq: %.0fMHz' % fMHz)

    files1 = []
    files2 = []
    gains = []
    rates = []
    fcens = []
    for d in dirs_byDev[dev1]:
        tmp = glob('%s/*.h5' % d)
        #print(tmp)
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
                files1.append(f)
                gains.append(gain)
                rates.append(rate)
                fcens.append(fcen)
                f2 = f.replace(dev1, dev2)
                files2.append(f2)

    nFile = len(files1)
    if (nFile == 0):
        sys.exit('no valid files for fcen=%.0f' % fMHz)
    #print('fcen=%.0fMHz, files:', files)
    print('nFile =', nFile)
    fcen = np.median(fcens)
    gain = np.median(gains)
    rate = np.median(rate)

    if (np.isclose(rate, 20e6)):
        bpname = 'median_20M'
        bpoff = -20.
    else:
        bpname = 'median'
        bpoff = 0.

    fname = 'analyze2_%.0f' % fMHz
    fout = fname + '.vish5'
    odir = fname + '.check'
    if (not os.path.isdir(odir)):
        call('mkdir %s'%odir, shell=True)

    # bandpass
    fbp = '%s/local/bin/blade_bandpass.h5' % home   # the freq-independent bandpass file, for nChan=1024 / 40MHz BW
    medbp = getData(fbp, bpname)          # bandpass in dB
    medbp += bpoff
    #linbp = 10**(medbp/10.)                 # in linear visbility scale
    gref = 20                               # the reference gain when bandpass was measured

    #bp = medbp + gain - gref    # bandpass under the current gain, in dB
    #linbp = 10**(bp/20)         # linear bandpass for voltage

    BWch = rate / nChan     # channel bandwidth in Hz
    Tint = 1./rate * nChan  # instantaneous integration time in sec

    fh51 = files1[0]
    idx = fh51.find('_dev')
    str_open = fh51[idx-13:idx]
    dt0 = datetime.strptime(str_open, '%y%m%d_%H%M%S')
    at0 = Time(dt0, format='datetime')  # local time
    ut0 = at0 - 8/24    # convert to UTC
    unix0 = ut0.to_value('unix')
    attrs1 = getAttrs(fh51)

    fh52 = files2[0]
    attrs2 = getAttrs(fh52)


    # update last attrs
    attrs['rate'] = list(tuple(attrs1['rate'])+tuple(attrs2['rate']))
    attrs['gain'] = list(tuple(attrs1['gain'])+tuple(attrs2['gain']))
    attrs['freq'] = list(tuple(attrs1['freq'])+tuple(attrs2['freq']))
    attrs['nInp'] = nInp
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
    keys = ['winSec', 'winNFT', 'winSpec', 'winVar', 'winCoeff', 'winMedCoeff', 'winAlign', 'winFile']
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
            if (k == 'winAlign'):
                savmask[k] = None
            if (k == 'winFile'):
                savmask[k] = None



    # processing files
    for i in range(nFile):
        fh51 = files1[i]
        fh52 = files2[i]
        print('processing:', fh51, fh52, '(%d/%d)'%(i+1,nFile), '%.1fMHz'%fMHz)
        savdata['winFile'].append(fh51)

        idx = fh51.find('_dev')
        str_open = fh51[idx-13:idx]
        dt_open = datetime.strptime(str_open, '%y%m%d_%H%M%S')
        off_open = (dt_open - dt0).total_seconds()  # offset in file opening time, in seconds

        fgains = np.array(attrs['gain'])    # assume gains do not vary with time

        bp = medbp.reshape((1,-1)) + fgains.reshape((4,1)) - gref    # bandpass under the current gain, in dB
        linbp = 10**(bp/20)         # linear bandpass for voltage, shape(nInp, nChan)

        ## processing windows
        for j in range(nWin):
            print('window:(%d/%d)'%(j+1, nWin))
            tWin = winOff + j*winSep

            sav_tWin = tWin + off_open
            if (sav_tWin in savdata['winSec']):
                print('skip existing window:', tWin, fh5)
                continue

            savdata['winSec'].append(sav_tWin)
            tag = 't%08.1f' % sav_tWin


            ## initialize chunks of a window
            chkdata = {}
            chkmask = {}
            for k in ['winAlign']:
                chkdata[k] = []
                chkmask[k] = None
            for k in ['winNFT', 'winSpec', 'winCoeff', 'winMedCoeff', 'winVar']:
                chkdata[k] = []
                chkmask[k] = []

            C0_chunk = int(tWin * rate)
            if (winDur <= tChunk):
                nSamp = int(winDur * rate)
                nFrame = nSamp // nChan
                nSamp = nFrame * nChan          # regularize the nSamp to multiples of nChan
                nChunk = 1
            else:
                nSamp = int(tChunk * rate)
                nFrame = nSamp // nChan
                nSamp = nFrame * nChan          # regularize the nSamp to multiples of nChan
                nChunk = int(winDur // tChunk)  # keep only integer number of full chunks

            for chki in range(nChunk):
                if (chki == 0):
                    C0 = C0_chunk
                else:
                    C0 = C0_prev + nSamp        # next chunk right after the previous chunk
                C0_prev = C0
                print('...chunk', chki, 'C0', C0)

                # if fixed is set, a non-zero med_off is provided
                C0_1 = 0
                C0_2 = 0
                if (med_off0 > 0):
                    C0_1 = med_off0
                elif (med_off0 < 0):
                    C0_2 = -med_off0
                
                print(dev1)
                data1, clock = load_cli_hdr(fh51, C0=C0+C0_1, nSamp=nSamp, verbose=0)
                print(dev2)
                data2, clock = load_cli_hdr(fh52, C0=C0+C0_2, nSamp=nSamp, verbose=0)
                data = np.ma.array([data1[:,0], data1[:,1], data2[:,0], data2[:,1]])
                data.mask = [data1.mask[:,0], data1.mask[:,1], data2.mask[:,0], data2.mask[:,1]]
                del data1, data2, clock
                gc.collect()

                rfi_mask = np.zeros(nChan, dtype=bool)
                if (auto_align):
                    print('auto-align:')
                    # check only cross-dev alignment
                    data.fill_value = 0.
                    winLen = data.shape[1]
                    midLen = winLen // 2
                    pairs = [[0,2],[1,3]]

                    #rfi_mask[768:] = True
                    #rfi_mask[:128] = True
                    #rfi_mask[250:275] = True
                    #rfi_mask[500:550] = True
                    #rfi_mask[750:825] = True
                    FTdata = np.ma.array(np.zeros((nInp, nFrame, nChan), dtype=complex), mask=False)
                    for k in range(nInp):
                        FTdata[k] = maskedFFT(data[k], nChan)
                    #xspec = np.ma.array(np.zeros((nBl, nChan), dtype=complex), mask=False)

                    ## auto-determine rfi_mask for the purpose of timestrema filtering
                    b = -1
                    for ai in range(nInp-1):
                        for aj in range(ai+1,nInp):
                            b += 1
                            xtmp = (FTdata[ai] * FTdata[aj].conjugate()).mean(axis=0)
                            #xspec[b] = xtmp.mean(axis=0)
                            clip_spec = sigma_clip(np.ma.abs(xtmp), sigma=4.0)
                            mask = clip_spec.mask
                            rfi_mask = np.logical_or(rfi_mask, mask)
                    
                    sfdata = np.ma.zeros(data.shape)
                    sfdata.mask = data.mask
                    sfdata.fill_value = 0.
                    for ii in range(nInp):
                        sfdata[ii] = streamFilter(data[ii], rfi_mask)

                    #norm2 = np.array([x.var() for x in data])
                    norm2 = np.array([x.var() for x in sfdata])


                    #peaks = []
                    offs  = []
                    for ipa, pa in enumerate(pairs):
                        ai, aj = pa
                        vi = sfdata[ai].filled()
                        vj = sfdata[aj].filled()
                        corr_arr = correlate(vi, vj, mode='same')
                        ncorr_arr = np.abs(corr_arr) / np.sqrt(norm2[ai]*norm2[aj]) / winLen
                        corr_peak = ncorr_arr.max()
                        corr_off  = ncorr_arr.argmax() - midLen
                        print('... ai,aj: peak, off', ai,aj,corr_peak,corr_off)
                        if (corr_peak > 0.001):
                            #peaks.append(corr_peak)
                            offs.append(corr_off)
                    if len(offs) > 0:
                        med_off = int(np.median(offs))
                    else:
                        med_off = 0
                    print('... med off:', med_off)
                    chkdata['winAlign'].append(med_off)
                    del sfdata
                    gc.collect()


                    C0_1 = 0
                    C0_2 = 0
                    if (med_off > 0):
                        C0_1 = med_off
                    elif (med_off < 0):
                        C0_2 = -med_off

                    # reload data
                    del data
                    gc.collect()
                    data1, clock = load_cli_hdr(fh51, C0=C0+C0_1, nSamp=nSamp, verbose=0)
                    data2, clock = load_cli_hdr(fh52, C0=C0+C0_2, nSamp=nSamp, verbose=0)
                    data = [data1[:,0], data1[:,1], data2[:,0], data2[:,1]]
                    del data1, data2, clock
                    gc.collect()

                else:   # no auto-align
                    chkdata['winAlign'].append(med_off0)


                FTdata = np.ma.array(np.zeros((nInp, nFrame, nChan), dtype=complex), mask=False)
                for k in range(nInp):
                    FTdata[k] = maskedFFT(data[k], nChan)
                    #print('main: (nmask)', np.count_nonzero(FTdata.mask[k]))
                del data
                gc.collect()

                nFTdata = FTdata / linbp.reshape((nInp,1,nChan))   # convert to dBm/ch
                #nFTspec = np.ma.abs(nFTdata).mean(axis=1)          # shape (nInp, nChan)
                nFTspec = np.ma.median(np.ma.abs(nFTdata), axis=1)          # shape (nInp, nChan)
                chkdata['winNFT'].append(nFTspec)                   # auto-correlation
                chkmask['winNFT'].append(nFTspec.mask)
                del FTdata
                gc.collect()



                ## cross-correlation
                xcorr = np.ma.array(np.zeros((nBl, nFrame, nChan), dtype=complex), mask=False)
                xspec = np.ma.array(np.zeros((nBl, nChan), dtype=complex), mask=False)
                coeff = np.ma.array(np.zeros((nBl, nFrame, nChan), dtype=complex), mask=False)
                coeffspec = np.ma.array(np.zeros((nBl, nChan), dtype=complex), mask=False)
                b = -1
                for ai in range(nInp-1):
                    for aj in range(ai+1, nInp):
                        b += 1

                        xtmp = nFTdata[ai] * nFTdata[aj].conjugate()    # shape(nFrame, nChan)
                        xcorr[b] = xtmp
                        xcorr.mask[b] = xtmp.mask

                        xsr   = np.ma.median(xtmp.real, axis=0)
                        xsi   = np.ma.median(xtmp.imag, axis=0)
                        taxtmp = xsr + 1.j*xsi
                        xspec[b] = taxtmp
                        xspec.mask[b] = taxtmp.mask

                        xc     = xtmp / (np.ma.abs(nFTdata[ai]) * np.ma.abs(nFTdata[aj]))
                        coeff[b] = xc
                        coeff.mask[b] = xc.mask

                        xcr   = np.ma.median(coeff[b].real, axis=0)
                        xci   = np.ma.median(coeff[b].imag, axis=0)
                        taxc = xcr + 1.j*xci
                        coeffspec[b] = taxc
                        coeffspec.mask[b] = taxc.mask
                        #print(taxtmp, taxc)

                chkdata['winSpec'].append(xspec)
                chkmask['winSpec'].append(xspec.mask)
                chkdata['winCoeff'].append(coeffspec)
                chkmask['winCoeff'].append(coeffspec.mask)

                med_coeff  = np.ma.median(np.ma.abs(coeffspec), axis=1)
                chkdata['winMedCoeff'].append(med_coeff)
                chkmask['winMedCoeff'].append(med_coeff.mask)

                var = xcorr.var(axis=1) # along nFrame
                var.mask = xcorr.mask.all(axis=1)
                chkdata['winVar'].append(var)   # shape(nBl, nChan)
                chkmask['winVar'].append(var.mask)

                del xcorr, coeff, nFTdata
                gc.collect()


            ## averaging the chunks in a window
            for k in chkdata.keys():
                if (savmask[k] is None):
                    tmp = np.array(chkdata[k]).mean(axis=0)
                    savdata[k].append(tmp)
                else:
                    tmp = np.ma.array(chkdata[k], mask=chkmask[k]).mean(axis=0)
                    savdata[k].append(tmp)
                    savmask[k].append(tmp.mask)
                    # recover the plot variables
                    if (k == 'winNFT'):
                        nFTspec = tmp
                    elif (k == 'winVar'):
                        var = tmp
                    elif (k == 'winSpec'):
                        xspec = tmp
                    elif (k == 'winCoeff'):
                        coeffspec = tmp
                    elif (k == 'winMedCoeff'):
                        med_coeff = tmp
            del chkdata, chkmask, tmp
            gc.collect()


            ## plotting for each window
            blname = []
            for ai in range(nInp-1):
                for aj in range(ai+1, nInp):
                    b += 1
                    corr = '%d-%d' % (ai, aj)
                    blname.append(corr)

            ## vs freq plots

            for pt in range(2):     # amp and pha plots
                if (pt == 0):
                    png1 = '%s/correlate_amp.%s.png' % (odir, tag)
                    ylab = 'power (dBm)'
                elif (pt == 1):
                    png1 = '%s/correlate_pha.%s.png' % (odir, tag)
                    ylab = 'phase (rad)'

                f1, s1 = plt.subplots(4,4,figsize=(16,12), sharex=True, sharey=True)
                b = -1
                for ii in range(nInp):
                    for jj in range(nInp):
                        ax = s1[ii,jj]
                        if (ii>jj):
                            ax.remove()
                        elif (ii == jj):
                            ax.set_xlabel('freq (MHz)')
                            ax.set_ylabel('power (dBm)')
                            ax.set_title('ch%d'%ii)
                            if (pt==0):
                                ax.plot(freq, toPower(nFTspec[ii]))
                        else:
                            b += 1
                            ax.set_title(blname[b])
                            if (pt==0):
                                ax.plot(freq, toPower(np.ma.abs(xspec[b]),mode='vis'))
                                ax.plot(freq[rfi_mask], toPower(np.ma.abs(xspec[b]), mode='vis')[rfi_mask])
                            elif (pt==1):
                                ax.plot(freq, np.ma.angle(xspec[b]))
                                ax.plot(freq[rfi_mask], np.ma.angle(xspec[b])[rfi_mask])
                                ax.set_ylim([-3.3,3.3])

                f1.tight_layout(rect=[0,0.03,1,0.95])
                f1.suptitle('%s, %s\n%s'%(fh51, fh52, tag))
                f1.savefig(png1)
                plt.close(f1)


            ## variance vs freq plot

            png2 = '%s/correlate_var.%s.png' % (odir, tag)
            f2, s2 = plt.subplots(nInp-1,nInp-1,figsize=(12,9), sharex=True, sharey=True)

            b = -1
            for ii in range(nInp-1):
                for jj in range(nInp-1):
                    ax = s2[ii,jj]
                    if (jj<ii):
                        ax.remove()
                    else:
                        b += 1
                        ax.plot(freq, var[b])
                        ax.set_yscale('log')
                        ax.set_title(blname[b])
                        if (ii==jj):
                            ax.set_ylabel('variance (mW^2)')
                            ax.set_xlabel('Freq (MHz)')

            f2.tight_layout(rect=[0,0.03,1,0.95])
            f2.suptitle('%s, %s\n%s'%(fh51, fh52, tag))
            f2.savefig(png2)
            plt.close(f2)


            ## corr-coefficient plot

            png3 = '%s/correlate_coeff.%s.png' % (odir, tag)
            f3, s3 = plt.subplots(nInp-1, nInp-1, figsize=(12,9), sharex=True, sharey=True)

            b = -1
            for ii in range(nInp-1):
                for jj in range(nInp-1):
                    ax = s3[ii,jj]
                    if (jj<ii):
                        ax.remove()
                    else:
                        b += 1
                        ax.plot(freq, np.ma.abs(coeffspec[b]), color='b', label='method1')
                        ax.text(0.05, 0.90, 'spec.med =%.3g'%med_coeff[b], transform=ax.transAxes)
                        ax.set_yscale('log')
                        ax.set_title(blname[b])
                        if (ii==jj):
                            ax.set_ylabel('corr-coeff')
                            ax.set_xlabel('Freq (MHz)')

            f3.tight_layout(rect=[0,0.03,1,0.95])
            f3.suptitle('%s, %s\n%s'%(fh51, fh52, tag))
            f3.savefig(png3)
            plt.close(f3)

            del xspec, coeffspec, var
            gc.collect()


    winSec = np.array(savdata['winSec'])
    adoneh5(fout, winSec, 'winSec')
    putAttrs(fout, {'note':'window time in seconds'}, dest='winSec')

    winNFT = np.ma.array(savdata['winNFT'], mask=savmask['winNFT'])
    adoneh5(fout, winNFT, 'winNFT')
    #putAttrs(fout, {'note':'window voltage data after bandpass calibration'}, dest='winNFT')
    putAttrs(fout, {'note':'window time-avg antenna auto-corr after bandpass calibration'}, dest='winNFT')

    winSpec = np.ma.array(savdata['winSpec'], mask=savmask['winSpec'])
    adoneh5(fout, winSpec, 'winSpec')
    putAttrs(fout, {'note':'window time-avg of the x-corr spectrum'}, dest='winSpec')

    adoneh5(fout, freq, 'freq')
    putAttrs(fout, {'note':'RF in MHz'}, dest='freq')

    winVar = np.ma.array(savdata['winVar'], mask=savmask['winVar'])
    adoneh5(fout, winVar, 'winVar')
    putAttrs(fout, {'note':'window time-avriance of the x-corr spectrum'}, dest='winVar')

    winCoeff = np.ma.array(savdata['winCoeff'], mask=savmask['winCoeff'])
    adoneh5(fout, winCoeff, 'winCoeff')
    putAttrs(fout, {'note':'window correlation coefficients taking the r-/i- median in time',
                    'median':savdata['winMedCoeff']}, dest='winCoeff')

    adoneh5(fout, blname, 'blname')

    sav_names = [x.encode() for x in savdata['winFile']]
    adoneh5(fout, sav_names, 'winFile')
    adoneh5(fout, savdata['winAlign'], 'winAlign')

    bmark2 = time.time()
    print('process %d windows in: %f sec' % (nFile, bmark2-bmark1))

