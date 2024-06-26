#!/usr/bin/env python

from read_cli2 import *
from loadh5 import *
from scipy.signal import correlate
from use_matplot import *
from astropy.time import Time
from astropy.stats import sigma_clip
import gc
import time

t1 = time.time()

## compare a number segments between two devices to determine the clock offset
## requires the header file for data reading
## also produce the visibility (as a function of spectral channels)
## save in a visibility file


inp = sys.argv[0:]
pg  = inp.pop(0)

# specify frequency
fMHz = 620
samprate = 40e6

# correlation windows
segSamp = 16384     # sample length per channel in dual-channel mode
winDur  = 0.5       # window length in seconds
nWin    = 24        # num of windows to check
winSep  = 300.      # 5min; window separation in seconds
winOff  = 0.    # 6h; overall window offset in seconds
                    # data started at 4am, skip 6h to start at 10am and analyze 2 hours
verbose = 0
# data directoies
idir    = {'dev7':'sun_%.0fMHz_dev7'%fMHz, 'dev8':'sun_%.0fMHz_dev8'%fMHz}    # key = device name (e.g. dev1, dev2, ...); value = directory (here 'dev1/', 'dev2/')
devs    = list(idir.keys())
nAnt    = int(2*len(devs))
nBl     = int(nAnt*(nAnt-1)/2)
# define cross-delay
# e.g.: 4 inputs are read into 0,1,2,3
#pairs   = [[0,1],[0,2],[1,3],[2,3]]
pairs   = [[0,2],[1,3]]
# spectrum
nchan   = 1024
chan = np.arange(nchan)

freq = np.fft.fftfreq(nchan, d=1/samprate)
freq = np.fft.fftshift(freq)    # IF in Hz
freq = freq / 1e6 + fMHz        # RF in MHz


# process all files in the idir specified above
# .h5 header files should have been produced
d0 = devs[0]    # start with one of the device
files = glob('%s/*.h5'%idir[d0])
files.sort()
nfiles = len(files)
if (nfiles < 1):
    sys.exit('no .h5 header files in %s'%idir[d0])

# split the data
nSplit = 1
iSplit = 0  # 0...nSplit-1
split = '%d_%d' % (nSplit, iSplit)

if (nfiles % nSplit == 0):
    nfiles1 = nfiles // nSplit
else:
    nfiles1 = nfiles // nSplit + 1
iStart = nfiles1*iSplit
iEnd = min(nfiles, iStart+nfiles1)

odir  = 'vis_plots_%.0fMHz-short' % fMHz
ofile = '%s/vis_split%s.h5' % (odir, split)
if (not os.path.isdir(odir)):
    call('mkdir -p %s'%odir, shell=True)

# misc setup
eye0 = np.identity(nAnt)
eye1 = eye0.copy()
eye1[nAnt-1, nAnt-1] = 0
eye2 = eye1.copy()
eye2[nAnt-2, nAnt-2] = 0


# one entry per file
sav_names = []
sav_epoch = []  # file open time in unix time
# n windows per file
sav_sec   = []  # window time in sec after file open
sav_off   = []  # device offset in clocks (offset>0 means dev0 needs to start later by this offset; dev0_C0=offset)
sav_peak  = []  # normalized correlation peak (>0.01 would be more reliable)
sav_mask  = []
sav_vis   = []  # auto_align on  
sav_vis2  = []  # auto_align off
sav_nullvis0  = []  # nulling-0
sav_nullvis1  = []  # nulling-1
sav_nullvis2  = []  # nulling-2
# start processing
for i in range(iStart, iEnd):
    dir0 = idir[devs[0]]
    fil0 = files[i]   #.h5 file
    print(fil0, '...', '(%d/%d)'%(i-iStart+1, nfiles1))
    dir1 = idir[devs[1]]
    fil1 = fil0.replace(dir0, dir1)
    fil1 = fil1.replace('_%s'%devs[0], '_%s'%devs[1])
    if (not os.path.isfile(fil1)):
        print('file not found: %s'%fil1, 'skipping ...')
        continue
    print('...', fil1)

    name = fil0.replace('%s/'%dir0, '')
    name = name.replace('_%s.h5'%devs[0], '')
    sav_names.append(name)

    tmp = name.split('_')
    dtstr = '_'.join(tmp[-2:])
    t0 = datetime.strptime(dtstr, '%y%m%d_%H%M%S')          # file open time
    topen = Time(t0, format='datetime').to_value('unix')    # in unix time
    sav_epoch.append(topen)

    attrs0 = getAttrs(fil0)
    rate0 = attrs0['rate'][0]

    # loop windows
    win_off  = []
    win_sec  = []
    win_peak = []
    win_mask = []
    win_vis  = []
    win_vis2 = []
    win_nullvis0 = []
    win_nullvis1 = []
    win_nullvis2 = []
    for j in range(nWin):
        ts = winSep * j + winOff        # window starting time in seconds
        win_sec.append(ts)
        cs = int(ts * rate0)            # starting tick
        winLen = int(winDur * rate0)    # window length in num of samples
        winLen -= winLen%2              # make it an even number
        print('...... starting time:', ts, 'sec')
        print('...... C0, nSamp:', cs, winLen)  # debugging
        midLen = winLen//2
        suptitle = '%s, win=%.2fsec'%(name, ts)

        # device lags
        #png1 = '%s/%s_win%d_corr.png' % (odir, name, j)
        #f1, s1 = plt.subplots(1,1,figsize=(8,6))
        # raw and masked spectra
        png2 = '%s/%s_win%d_spec.png' % (odir, name, j)
        f2, sub2 = plt.subplots(2,2,figsize=(12,9), sharex=True, sharey=True)
        s2 = sub2.flatten()
        # auto-align on vis plots
        png3 = '%s/%s_win%d_vispha.png' % (odir, name, j)
        f3, sub3 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        png4 = '%s/%s_win%d_visamp.png' % (odir, name, j)
        f4, sub4 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        # auto-align off vis plots
        png5 = '%s/%s_win%d_vis2pha.png' % (odir, name, j)
        f5, sub5 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        png6 = '%s/%s_win%d_vis2amp.png' % (odir, name, j)
        f6, sub6 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        # norma-spec and eigenmodes and nulling
        png7 = '%s/%s_win%d_eigenmode.png' % (odir, name, j)
        f7, sub7 = plt.subplots(2,1,figsize=(10,9), sharex=True, sharey=True)
        png8 = '%s/%s_win%d_nulling.png' % (odir, name, j)
        f8, sub8 = plt.subplots(2,2,figsize=(12,9), sharex=True, sharey=True)
        s8 = sub8.flatten()
        # eigenmode removal vis plots
        png9 = '%s/%s_win%d_vis3pha.png' % (odir, name, j)
        f9, sub9 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        png10 = '%s/%s_win%d_vis3amp.png' % (odir, name, j)
        f10, sub10 = plt.subplots(3,3,figsize=(12,9), sharex=True, sharey=True)
        for ai in range(3):
            for aj in range(3):
                if (ai>aj):
                    sub3[ai,aj].remove()
                    sub4[ai,aj].remove()
                    sub5[ai,aj].remove()
                    sub6[ai,aj].remove()
                    sub9[ai,aj].remove()
                    sub10[ai,aj].remove()
        s3 = sub3.flatten()
        s4 = sub4.flatten()
        s5 = sub5.flatten()
        s6 = sub6.flatten()
        s9 = sub6.flatten()
        s10 = sub5.flatten()

        data0, clock0 = load_cli_hdr(fil0, C0=cs, nSamp=winLen, verbose=verbose)
        #print('main: data0', data0.shape, data0.mask.shape)
        data1, clock1 = load_cli_hdr(fil1, C0=cs, nSamp=winLen, verbose=verbose)
        #print('main: data1', data1.shape, data1.mask.shape)
        data = np.ma.array([data0[:,0], data0[:,1], data1[:,0], data1[:,1]], shrink=False)
        data.mask = [data0.mask[:,0], data0.mask[:,1], data1.mask[:,0], data1.mask[:,1]]
        data.fill_value = 0.    # for correlation
        #norm = np.array([x.var() for x in data])
        #print('main: data', data.shape, data.mask.shape)
        del data0, data1, clock0, clock1
        gc.collect()

        # filter out strong RFI before determine offset
        # assume no offset
        FTdata = []
        frames = []
        for ch in range(4):
            chFT = maskedFFT(data[ch])  # default to 1024 ch
            frames.append(chFT.shape[0])
            FTdata.append(chFT)
        nframe = min(frames)
        for ch in range(4):
            FTdata[ch] = FTdata[ch][:nframe]
        FTdata = np.ma.array(FTdata)
        del chFT
        gc.collect()

        # covariance matrix and eigenmodes
        print('process eigen modes ...')
        Cov, norm = makeCov(FTdata, scale=True)
        W, V = Cov2Eig(Cov)
        nFTdata = FTdata / (norm.reshape((nAnt,1,nchan)))
        del FTdata
        gc.collect()

        # norm-spec and eigenmodes
        for ch in range(nAnt):
            cc = 'C%d'%ch
            sub7[0].plot(freq, 20*np.log10(np.abs(nFTdata[ch]).mean(axis=0)), color=cc, label='Ant%d'%ch)
            sub7[1].plot(freq, 10*np.log10(W[:,ch]), color=cc, label='W%d'%ch)
        sub7[0].set_title('normalized spec.')
        sub7[0].set_ylabel('power (dB)')
        sub7[0].legend()
        sub7[1].set_title('eigenmodes')
        sub7[1].set_ylabel('power (dB)')
        sub7[1].set_xlabel('freq (MHz)')
        sub7[1].legend()
        f7.tight_layout(rect=[0,0.03,1,0.95])
        f7.subplots_adjust(hspace=0.05)
        f7.savefig(png7)
        plt.close(f7)


        # diagnostic plot, masked spec
        rfi_mask = []
        for ch in range(4):
            var = nFTdata[ch].var(axis=0)    # real numbers
            clip_var = sigma_clip(var, sigma=10)
            rfi_mask.append(clip_var.mask)
        rfi_mask = np.any(rfi_mask, axis=0)
        win_mask.append(rfi_mask)

        for ch in range(4):
            ax = s2[ch]
            ax.plot(chan, 20.*np.log10(np.abs(nFTdata[ch]).mean(axis=0)), label='original')
            ax.plot(chan[~rfi_mask], 20.*np.log10(np.abs(nFTdata[ch]).mean(axis=0)[~rfi_mask]), label='masked')
            if (ch==0):
                ax.legend()
            ax.set_title('ch%d'%ch)
            ax.set_xlabel('chan')
            ax.set_ylabel('power (dB)')
        f2.tight_layout(rect=[0,0.03,1,0.95])
        f2.suptitle(suptitle)
        f2.savefig(png2)
        plt.close(f2)


        # nulling and diagnostic plots
        #null0 = np.zeros_like(nFTdata)   # nulling of nothing (sanity check)
        null1 = np.zeros_like(nFTdata)   # nulling of V3
        null2 = np.zeros_like(nFTdata)   # nulling of V3+V2
        for ii in range(nchan):
            Vi = V[ii]
            Vih = Vi.conjugate().T
            sane = np.dot(Vi, np.dot(eye0, Vih))
            rem1 = np.dot(Vi, np.dot(eye1, Vih))
            rem2 = np.dot(Vi, np.dot(eye2, Vih))
            #null0[:,:,ii] = np.tensordot(sane, nFTdata[:,:,ii], axes=(1,0))
            null1[:,:,ii] = np.tensordot(rem1, nFTdata[:,:,ii], axes=(1,0))
            null2[:,:,ii] = np.tensordot(rem2, nFTdata[:,:,ii], axes=(1,0))
        #null0.mask = nFTdata.mask
        null1.mask = nFTdata.mask
        null2.mask = nFTdata.mask




        # auto-aling off visibilities
        snap_vis = []   # raw vis
        #null_vis0 = []  # null0
        null_vis1 = []  # null1
        null_vis2 = []  # null2
        b = -1
        for ai in range(3):
            for aj in range(ai+1,4):
                b += 1
                corr = '%d-%d' % (ai,aj)

                # raw
                blvis = nFTdata[ai]*nFTdata[aj].conjugate()
                avgblvis = blvis.mean(axis=0)
                snap_vis.append(avgblvis)
                ax = sub5[ai,aj-1]
                ax.plot(chan, np.angle(avgblvis))
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.set_xlabel('chan')
                    ax.set_ylabel('phase (rad)')
                ax = sub6[ai,aj-1]
                ax.plot(chan[~rfi_mask], 10.*np.log10(np.abs(avgblvis[~rfi_mask]))) # note, 10*log10(vis), not 20*log10(vis)
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.set_xlabel('chan')
                    ax.set_ylabel('power (dB)')

                # nulling
                #blvis0 = null0[ai]*null0[aj].conjugate()
                #avgblvis0 = blvis0.mean(axis=0)
                #null_vis0.append(avgblvis0)
                blvis1 = null1[ai]*null1[aj].conjugate()
                avgblvis1 = blvis1.mean(axis=0)
                null_vis1.append(avgblvis1)
                blvis2 = null2[ai]*null2[aj].conjugate()
                avgblvis2 = blvis2.mean(axis=0)
                null_vis2.append(avgblvis2)
                ax = sub9[ai,aj-1]
                #ax.plot(chan, np.angle(avgblvis0), label='null0')
                ax.plot(chan, np.angle(avgblvis1), label='null1')
                ax.plot(chan, np.angle(avgblvis2), label='null2')
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.legend()
                    ax.set_xlabel('chan')
                    ax.set_ylabel('phase (rad)')
                ax = sub10[ai,aj-1]
                #ax.plot(chan, 10.*np.log10(np.abs(avgblvis0)), label='null0') # note, 10*log10(vis), not 20*log10(vis)
                ax.plot(chan, 10.*np.log10(np.abs(avgblvis1)), label='null1') # note, 10*log10(vis), not 20*log10(vis)
                ax.plot(chan, 10.*np.log10(np.abs(avgblvis2)), label='null2') # note, 10*log10(vis), not 20*log10(vis)
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.legend()
                    ax.set_xlabel('chan')
                    ax.set_ylabel('power (dB)')

                #del blvis, blvis0, blvis1, blvis2
                del blvis, blvis1, blvis2
                gc.collect()

        win_vis2.append(snap_vis)
        #win_nullvis0.append(null_vis0)
        win_nullvis1.append(null_vis1)
        win_nullvis2.append(null_vis2)

        f5.tight_layout(rect=[0,0.03,1,0.95])
        f5.suptitle('%s, auto_align=off'%suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f5.savefig(png5)
        plt.close(f5)
        f6.tight_layout(rect=[0,0.03,1,0.95])
        f6.suptitle('%s, auto_align=off'%suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f6.savefig(png6)
        plt.close(f6)

        f9.tight_layout(rect=[0,0.03,1,0.95])
        f9.suptitle('%s, nulling_vis'%suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f9.savefig(png9)
        plt.close(f9)
        f10.tight_layout(rect=[0,0.03,1,0.95])
        f10.suptitle('%s, nulling_vis'%suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f10.savefig(png10)
        plt.close(f10)



        ## skip auto-align in this version
        continue


        # filter out the strong RFI            
        filtered = []
        for ch in range(4):
            tmp = streamFilter(data[ch], rfi_mask)
            filtered.append(tmp)
        rdata = np.ma.array(filtered)
        rdata.fill_value = 0.    # for correlation
        norm2 = np.array([x.var() for x in rdata])
        winLen2 = rdata.shape[1]
        midLen2 = winLen2 // 2
        print('winLen2, midLen2:', winLen2, midLen2)
        del tmp, filtered
        gc.collect()


        # find correlation lags after filtering
        # loop pairs
        pair_off = []
        pair_peak = []
        for ipa, pa in enumerate(pairs):
            ai, aj = pa
            print('......... pair:', ai, aj)
            corr = correlate(rdata[ai].filled(), rdata[aj].filled(), mode='same') # complex, raw
            ncorr = np.abs(corr) / np.sqrt(norm2[ai]*norm2[aj]) / winLen2  # real, normalized
            peak = ncorr.max()
            off0 = np.argmax(ncorr)
            off = off0 - midLen2
            print('......... peak, off', peak, off)
            pair_off.append(off)
            pair_peak.append(peak)

            xoff = np.arange(len(corr)) - midLen2
            s1.plot(xoff, ncorr, label='[%d,%d] off=%d, peak=%.4f'%(ai,aj,off,peak))
        del rdata
        gc.collect()

        # diagnostic, correlation lags
        s1.legend()
        s1.set_xlabel('corr_lag (clock)')
        s1.set_ylabel('corr. coefficient')
        s1.set_title(suptitle)
        f1.savefig(png1)
        plt.close(f1)

        win_off.append(pair_off)
        win_peak.append(pair_peak)


        ## calculate visibilities
        med_off = int(np.mean(pair_off))
        shifted = []
        if (med_off > 0):
            shifted.append(data[0, med_off:])
            shifted.append(data[1, med_off:])
            shifted.append(data[2, :-med_off])
            shifted.append(data[3, :-med_off])
        elif (med_off < 0):
            shifted.append(data[0, :med_off])
            shifted.append(data[1, :med_off])
            shifted.append(data[2, -med_off:])
            shifted.append(data[3, -med_off:])
        else:
            shifted = data
        data = np.ma.array(shifted)

        FTdata = []
        for ch in range(4):
            FTdata.append(maskedFFT(data[ch]))
        FTdata = np.ma.array(FTdata)

        snap_vis = []
        b = -1
        for ai in range(3):
            for aj in range(ai+1,4):
                b += 1
                corr = '%d-%d' % (ai,aj)

                blvis = FTdata[ai]*FTdata[aj].conjugate()
                avgblvis = blvis.mean(axis=0)
                snap_vis.append(avgblvis)

                ax = sub3[ai,aj-1]
                ax.plot(chan, np.angle(avgblvis))
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.set_xlabel('chan')
                    ax.set_ylabel('phase (rad)')
                ax = sub4[ai,aj-1]
                ax.plot(chan[~rfi_mask], 10.*np.log10(np.abs(avgblvis[~rfi_mask]))) # note, 10*log10(vis), not 20*log10(vis)
                ax.text(0.02,0.95,corr,transform=ax.transAxes)
                if (b==0):
                    ax.set_xlabel('chan')
                    ax.set_ylabel('power (dB)')

        del FTdata
        gc.collect()

        win_vis.append(snap_vis)

        f3.tight_layout(rect=[0,0.03,1,0.95])
        f3.suptitle('%s, auto_align=on'%suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f3.savefig(png3)
        plt.close(f3)
        f4.tight_layout(rect=[0,0.03,1,0.95])
        f4.suptitle(suptitle)
        plt.subplots_adjust(wspace=0, hspace=0)
        f4.savefig(png4)
        plt.close(f4)




    sav_sec.append(win_sec)
    #sav_off.append(win_off)
    #sav_peak.append(win_peak)
    sav_mask.append(win_mask)
    #sav_vis.append(win_vis)
    sav_vis2.append(win_vis2)
    #sav_nullvis0.append(win_nullvis0)
    sav_nullvis1.append(win_nullvis1)
    sav_nullvis2.append(win_nullvis2)



    #if (i>=2):
    #    break


sav_names = [x.encode() for x in sav_names] # h5py does not support saving the unicode strings, need to encode first
sav_names = np.array(sav_names)
sav_epoch = np.array(sav_epoch)
sav_sec   = np.array(sav_sec)
#sav_off   = np.array(sav_off)
#sav_peak  = np.array(sav_peak)
sav_mask  = np.array(sav_mask)
#sav_vis   = np.ma.array(sav_vis)    # (nfiles, nwin, nbl, nchan)
sav_vis2  = np.ma.array(sav_vis2)   # (nfiles, nwin, nbl, nchan)
#sav_nullvis0  = np.ma.array(sav_nullvis0)   # (nfiles, nwin, nbl, nchan)
sav_nullvis1  = np.ma.array(sav_nullvis1)   # (nfiles, nwin, nbl, nchan)
sav_nullvis2  = np.ma.array(sav_nullvis2)   # (nfiles, nwin, nbl, nchan)

# save to file
adoneh5(ofile, sav_names, 'names')
adoneh5(ofile, sav_epoch, 'unixtime')
adoneh5(ofile, sav_sec,   'win_sec')
#adoneh5(ofile, sav_off,   'win_off')
#adoneh5(ofile, sav_peak,  'win_peak')
adoneh5(ofile, sav_mask,  'win_mask')
#adoneh5(ofile, sav_vis,   'win_vis')    # auto_align on
adoneh5(ofile, sav_vis2,  'win_vis2')   # auto_align off
#adoneh5(ofile, sav_nullvis0,  'win_nullvis0')   # nulling
adoneh5(ofile, sav_nullvis1,  'win_nullvis1')   # nulling
adoneh5(ofile, sav_nullvis2,  'win_nullvis2')   # nulling

attrs = {}
attrs['nWin']   = nWin
attrs['winDur'] = winDur
attrs['winOff'] = winOff
attrs['pairs']  = pairs
attrs['devs']   = devs
attrs['idir']   = [idir[x] for x in devs]
putAttrs(ofile, attrs)


# the final tally
t2 = time.time()
print('process files:', nfiles1, 'time:', t2-t1, 'secs')


