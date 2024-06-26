#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.time import Time
from astropy.stats import sigma_clip

from loadh5 import *
from util_func import *
from pyplanet import *
from delay_func2 import *
from array_func import blsubplots

DB = loadDB()

inp = sys.argv[0:]
pg  = inp.pop(0)

#nChan = 1024   # get from .vish5 file
ifmin   = 0
ifmax   = 0
rfmin   = 0
rfmax   = 0
aconf   = None
ai      = 0
aj      = 1
srcflux = 0.
f410    = 0.
f610    = 0.
do_fitClkOff = False
do_check = False
blfit   = 0         # default disable fitting baslines
anglim  = 30.
angsep  = 3.
poslim  = 0.5
possep  = 0.1
gant    = 6.0       # antenna max gain in dB, black LPDA
Ehphw   = 30.       # E-plane half-power half width in deg
Hhphw   = 60.       # H-plane half-power half width in deg
dellim  = 150.
delsep  = 0.5
body    = 'sun'


def atten(x, hphw):
    '''
    Gaussian attenuation (linear) given offset x and the half-power half width
    '''
    sig = hphw / np.sqrt(2.*np.log(2.))
    return np.exp(-x**2/2./sig**2)


usage = '''
plot all windows of a given .vish5 file
multiple files can be given
each file has a corresponding plot

syntax:
    %s <.vish5_files> [options]

options are:
    -b BODY         # specify the target name for fringe estimate
                    # (default: %s)

    --flux JY       # specify the source flux in Jy for SEFD calculation
                    # manually input the flux appropriate for the frequency
                    # NEW: use --interp is preferred

    --interp f410 f610
                    # alternatively, specify the 410MHz and 610MHz flux in Jy from the website
                    # flux is interpolated to the frequency under analysis

    --if MIN MAX    # IF range (in Hz) to use for time plot
                    # (both default to 0, for using the full range)

    --aconf ANT_CONT # specify an antenna configuration file
                    # default is none

    --gant GAIN     # antenna max gain in dB (%.1f dB)
                    # do not change --gant unless you are using a different antenna

    ## to select only a subset of antennas for processing
    --ai AI         # specify the antenna numbers (only needed if --aconf is used)
    --aj AJ         # (AI, AJ are default to 0, 1)

    ## for fitting clock offset
    --fit           # fit and remove the Clock Offset
    --plot          # plot the Clock Offset phase correction

    --del 'LIMIT SEP' # set the +/-LIMIT and SEP for delay (in samples)
                    # default: +/-150, 0.5

    ## for fitting antenna positions
    --init 'PARAMs' # assuming a linear array:
                    # param 1: direction of the array (in degree)
                    # param 2+: distance from Ant0 (in meter)
    --ang 'LIMIT SEP' # set the +/-LIMIT and SEP for angle (in degree)
                    # default: +/-30, 3
    --pos 'LIMIT SEP' # set the +/-LIMIT and SEP for antenna positions (in meter)
                    # default: +/-0.5, 0.1

''' % (pg, body, gant)

if (len(inp) < 1):
    sys.exit(usage)

files = []
while (inp):
    k = inp.pop(0)
    if (k == '--if'):
        ifmin = float(inp.pop(0))
        ifmax = float(inp.pop(0))
    elif (k == '-b'):
        body = inp.pop(0)
    elif (k == '--aconf'):
        aconf = inp.pop(0)
    elif (k == '--ai'):
        ai = int(inp.pop(0))
    elif (k == '--aj'):
        aj = int(inp.pop(0))
    elif (k == '--flux'):
        srcflux = float(inp.pop(0))
    elif (k == '--interp'):
        f410 = float(inp.pop(0))
        f610 = float(inp.pop(0))
        print(f410, f610)
    elif (k == '--gant'):
        gant = float(inp.pop(0))
    elif (k == '--hphw'):
        hphw = float(inp.pop(0))
    elif (k == '--fit'):
        do_fitClkOff = True
    elif (k == '--plot'):
        do_check = True
    elif (k == '--init'):
        tmp = inp.pop(0).split()
        ang0 = float(tmp.pop(0))
        pos0 = [float(x) for x in tmp]
        blfit = 1   # enable fitting basline
    elif (k == '--ang'):
        tmp = inp.pop(0).split()
        anglim = float(tmp[0])
        angsep = float(tmp[1])
    elif (k == '--pos'):
        tmp = inp.pop(0).split()
        poslim = float(tmp[0])
        possep = float(tmp[1])
    elif (k == '--del'):
        tmp = inp.pop(0).split()
        dellim = float(tmp[0])
        delsep = float(tmp[1])
    elif (k.startswith('-')):
        sys.exit('unknown option: %s'%k)
    else:
        files.append(k)
print('file:', files)



if (aconf is None):
    do_model = False
else:
    do_model = True
    tmp = os.path.basename(aconf)
    tmp2 = tmp.split('_')
    site = tmp2[0].lower()
    antXYZ = np.loadtxt(aconf, usecols=(1,2,3))

G0 = 10.**(gant/10.)

for fvis in files:
    if (not os.path.isfile(fvis)):
        print('file not found: %s --> skip' % fvis)
        continue
    print('processing:', fvis)

    cdir = fvis.replace('.vish5', '.check')
    if (not os.path.isdir(cdir)):
        call('mkdir -p %s'%cdir, shell=True)
    odir = fvis.replace('.vish5', '.output')
    if (not os.path.isdir(odir)):
        call('mkdir -p %s'%odir, shell=True)


    attrs = getAttrs(fvis)
    fcen = attrs['freq'][0]         # ch0 central freq in Hz
    rate = attrs['rate'][0]         # ch0 sampling rate in Hz

    rfreq = getData(fvis, 'freq')   # RF in MHz
    nChan = len(rfreq)

    if (f410>0 and f610>0):         # interpolate the 
        xs = (rfreq - 410)/(610-410)
        srcflux2 = f410 + (f610-f410)*xs
        # --interp takes precedence over --flux
        print('fcen=%.0fMHz'%(fcen/1e6), 'flux=%.2eJy'%srcflux)
        srcflux = srcflux2.mean()
    else:
        srcflux2 = np.array([srcflux])


    # file-open time
    unix_utc_open = attrs.get('unix_utc_open')
    if (unix_utc_open is None):
        fbase = os.path.basename(fvis)
        ii = fbase.find('_dev')
        dtstr = fbase[ii-13:ii]
        #print(fvis, ii, dtstr)
        topen = datetime.strptime(dtstr, '%y%m%d_%H%M%S')   # local time
        atopen = Time(topen, format='datetime')
        atopen -= 8/24  # convert to UTC
    else:
        atopen = Time(unix_utc_open, format='unix')
    print(atopen)

    ## unix_utc_open / atopen was converted from 'local' time to UTC, assuming UTC+8
    ## will be corrected if doing model
    tz_corr = 0.
    if (do_model):
        if (site == 'tnro'):    # UTC+7
            #tz_corr = 1.        # add 1 hr back to the UTC
            tz_corr = 0.        # turns out, system time was still Taipei time, so no correction
        atopen += tz_corr/24.
    print(atopen)



    nft  = getData(fvis, 'winNFT')  # shape: nWin, nInp, nChan
    spec1 = getData(fvis, 'winSpec') # shape: nWin, nBl, nChan
    sec  = getData(fvis, 'winSec')  # shape: nWin
    var  = getData(fvis, 'winVar')  # shape: nWin, nBl, nChan
    coeff1 = getData(fvis, 'winCoeff') # shape: nWin, nBl, nChan
    nWin = len(sec)

    # window time
    tmp = atopen + sec/86400        # convert to days
    tWin = tmp.to_datetime()


    if (len(spec1.shape) == 2): # backward-compatible with single_dev_vis2.py
        spec1  = spec1.reshape((nWin, 1, nChan))
        coeff1 = coeff1.reshape((nWin, 1, nChan))
        var    = var.reshape((nWin, 1, nChan))

    nInp = nft.shape[1]
    nBl = int(nInp*(nInp-1)/2)
    pairs = []
    for ai in range(nInp-1):
        for aj in range(ai+1,nInp):
            pairs.append([ai,aj])

    freq  = np.fft.fftfreq(nChan, d=1./rate)
    freq  = np.fft.fftshift(freq)
    freq  /= 1e6    # IF in MHz
    if (ifmin == 0):
        ifmin = freq.min()
    if (ifmax == 0):
        ifmax = freq.max()

    rfmin = ifmin + fcen
    rfmax = ifmax + fcen


    if (do_fitClkOff):
        ## find and remove clock offset
        print('removing clock offset ...')
        savClkOff = []
        spec3 = spec1.copy()
        coeff3 = coeff1.copy()
        params = np.arange(-10,10,0.03)
        #for i in range(1,2):   # testing
        for i in range(nWin):
            tag = 't%08.1f' % sec[i]
            ClkOff = []
            for b in range(nBl):
                vispha = np.ma.angle(spec1[i,b])
                viserr = np.ones(nChan) * 0.3
                like = lnlike(params, vispha, viserr, nChan)
                maxlike = like.max()
                maxidx  = like.argmax()
                maxpara = params[maxidx]
                #print(b, maxpara)
                phimod = phimodel(maxpara, nChan)
                spec3[i,b] = spec1[i,b] / np.exp(1.j*phimod)
                coeff3[i,b] = coeff1[i,b] / np.exp(1.j*phimod)
                ClkOff.append(maxpara)
            savClkOff.append(ClkOff)

            if (do_check):
                png = '%s/correlate_pha2_%s.png' % (cdir, tag)
                fig, sub = blsubplots(na=nInp, figsize=(12,9), squeeze=False)
                b = -1
                for ai in range(nInp-1):
                    for aj in range(ai+1,nInp):
                        b += 1
                        ax = sub[ai,aj-1]
                        ax.plot(freq, np.ma.angle(spec1[i,b]), label='original')
                        ax.plot(freq, np.ma.angle(spec3[i,b]), label='unwrap')
                        ax.set_xlabel('freq (MHz)')
                        ax.set_ylabel('phase (rad)')
                        ax.legend()
                fig.savefig(png)
                plt.close(fig)

        savClkOff = np.array(savClkOff)
        fig, sub = blsubplots(na=nInp, figsize=(12,9), squeeze=False)
        b = -1
        for ai in range(nInp-1):
            for aj in range(ai+1, nInp):
                b += 1
                ax = sub[ai,aj-1]
                ax.plot(sec, savClkOff[:,b], marker='.')
                ax.set_xlabel('time (sec)')
                ax.set_ylabel('clock off (sample)')
        fig.suptitle('%s'%(fvis,))
        fig.savefig('%s/clock_off.png' % cdir)
        plt.close(fig)

        adoneh5(fvis, savClkOff, 'winClkOff')
        adoneh5(fvis, spec3, 'winSpecClkCorr')
        adoneh5(fvis, coeff3, 'winCoeffClkCorr')

        spec2 = spec3 - spec3.mean(axis=0)
        coeff2 = coeff3 - coeff3.mean(axis=0)
        print('... clock offset done.')

    else:
        spec3 = getData(fvis, 'winSpecClkCorr')
        if (spec3 is None):
            spec2 = spec1 - spec1.mean(axis=0)
        else:
            spec2 = spec3 - spec3.mean(axis=0)
        coeff3 = getData(fvis, 'winCoeffClkCorr')
        if (coeff3 is None):
            coeff2 = coeff1 - coeff1.mean(axis=0)
        else:
            coeff2 = coeff3 - coeff3.mean(axis=0)
    adoneh5(fvis, spec2, 'winSpecSub')
    adoneh5(fvis, coeff2, 'winCoeffSub')
        

    xx = sec
    yy = rfreq
    X,Y = np.meshgrid(xx, yy, indexing='xy')
    tt = tWin
    T,Y = np.meshgrid(tt, yy, indexing='xy')
    ## plot waterfall of real/imag/amp/pha for spec1, spec2 and coeff1, coeff2
    if (do_check):

        for b, p in enumerate(pairs):
            ai, aj = p
            bname = '%d-%d'%(ai,aj)

            for pp in range(4):
                if (pp == 0):
                    png = '%s/waterfall_bl%s_spec_raw.png' % (cdir,bname)
                    pp_spec = spec1
                elif (pp == 1):
                    png = '%s/waterfall_bl%s_spec_sub.png' % (cdir,bname)
                    pp_spec = spec2
                elif (pp == 2):
                    png = '%s/waterfall_bl%s_coeff_raw.png' % (cdir,bname)
                    pp_spec = coeff1
                elif (pp == 3):
                    png = '%s/waterfall_bl%s_coeff_sub.png' % (cdir,bname)
                    pp_spec = coeff2

                fig, sub = plt.subplots(2,2,figsize=(15,10))

                # real
                ax = sub[0,0]
                if (pp > 1):
                    s=ax.pcolormesh(X,Y,pp_spec.real[:,b,:].T,shading='auto')
                else:
                    s=ax.pcolormesh(X,Y,pp_spec.real[:,b,:].T,shading='auto', norm=colors.SymLogNorm(linthresh=1e3))
                cb=plt.colorbar(s, ax=ax)
                cb.set_label('linear')
                ax.set_xlabel('time (sec)')
                ax.set_ylabel('freq (MHz)')
                ax.set_title('Real')
                # imag
                ax = sub[0,1]
                if (pp > 1):
                    s=ax.pcolormesh(X,Y,pp_spec.imag[:,b,:].T,shading='auto')
                else:
                    s=ax.pcolormesh(X,Y,pp_spec.imag[:,b,:].T,shading='auto', norm=colors.SymLogNorm(linthresh=1e3))
                cb=plt.colorbar(s, ax=ax)
                cb.set_label('linear')
                ax.set_xlabel('time (sec)')
                ax.set_ylabel('freq (MHz)')
                ax.set_title('Imag')
                # amp
                ax = sub[1,0]
                s=ax.pcolormesh(X,Y,10*np.ma.log10(np.ma.abs(pp_spec[:,b,:])).T,shading='auto')
                cb=plt.colorbar(s, ax=ax)
                cb.set_label('10*log10(abs())')
                ax.set_xlabel('time (sec)')
                ax.set_ylabel('freq (MHz)')
                ax.set_title('Power')
                # pha
                ax = sub[1,1]
                s=ax.pcolormesh(X,Y,np.ma.angle(pp_spec[:,b,:]).T,shading='auto')
                cb=plt.colorbar(s, ax=ax)
                cb.set_label('phase (rad)')
                ax.set_xlabel('time (sec)')
                ax.set_ylabel('freq (MHz)')
                ax.set_title('Phase')

                fig.tight_layout(rect=[0,0.03,1,0.95])
                fig.suptitle('%s'%png)
                fig.savefig(png)
                plt.close(fig)
        


    # spectral avg
    fa_r = np.ma.median(spec2.real, axis=2)
    fa_i = np.ma.median(spec2.imag, axis=2)
    fa_spec2 = fa_r + 1.j*fa_i  # shape (nWin, nBl)


    # window time
    tmp = atopen + sec/86400        # convert to days
    tWin = tmp.to_datetime()

    fsel = np.logical_and(freq>=ifmin, freq<=ifmax)
    avgIF = freq[fsel].mean()
    avgRF = avgIF + fcen/1e6        # effective central freq in MHz
    lam = 2.998e8 / (avgRF * 1e6)   # wavelength in meters
    flam = 2.998e8 / (rfreq * 1e6)
    flam0 = 2.998e8 / (400e6)       # a reference lambda at 400MHz

    chBW = rate / nChan             # channel bandwidth in Hz
    tInt = 1./rate * nChan          # integration time per spectrum in sec
    scale = np.sqrt(1.*chBW*tInt)   # conversion between variance and SEFD

    #G0 = 10.                        # on-axis gain (10dB --> 10)
    # fiducial num
    Tsys0 = 50.
    lam0  = 0.5
    Aeff0 = lam0**2 * G0 / (4*np.pi)
    SEFD0 = 2.*1.38e-23*1e26*Tsys0/Aeff0/scale  # converted to Jy
    #print('fiducial SEFD (150K, 0.5m):', SEFD0/1e6, 'MJy')

    Aeff = lam**2 * G0 / (4.*np.pi)
    scale2 = Aeff * scale / (2*1.38e-23) * 1e-26 # Tsys = (S/Jy / SNR) * scale2
    SEFD2 = Tsys0 / scale2
    print('fiducial SEFD (%.0fK, %.2fm):'%(Tsys0,lam), '%.1fMJy'%(SEFD2/1e6))
    fscale2 = flam**2 * G0 / (4.*np.pi) * scale / (2*1.38e-23) * 1e-26
    

    ## SEFD plot
    if (srcflux > 0):
        medAmp = np.ma.median(np.ma.abs(spec2)) # median level
        calFac = medAmp / srcflux
        sefd = np.ma.sqrt(var) / calFac * scale # SEFD in Jy
        sefd.fill_value = 0.
        par = np.percentile(sefd.filled(), [0, 15, 85, 100], axis=0)
    else:
        par = np.zeros((4, nChan))

    if (srcflux > 0 and False):
        f3, s3 = plt.subplots(1,1)
        png3 = '%s.SEFD.png' % fvis
        suptitle = 'file: %s, %.0fMHz, nWin=%d' % (fvis, fcen/1e6, nWin)

        ax = s3
        ax.fill_between(freq, par[0], par[3], alpha=0.2, color='b', label='min/max')
        ax.fill_between(freq, par[1], par[2], alpha=0.6, color='b', label='15%/85% pct')
        ax.grid(True, which='both', axis='y')
        ax.grid(True, which='major', axis='x')
        ax.legend()
        ax.set_yscale('log')
        ax.set_xlabel('freq (MHz)')
        ax.set_ylabel('SEFD (Jy)')
        ax.set_title(suptitle)
        f3.tight_layout()
        f3.savefig(png3)
        plt.close(f3)



    if (do_model):
        ut0 = tWin[0]
        b, obs = obsBody(body, time=ut0, site=site, retOBS=True, DB=DB)

        az = []
        el = []
        for ti in tWin:
            obs.date = ti
            b.compute(obs)
            az.append(b.az)
            el.append(b.alt)
        az = np.array(az)
        el = np.array(el)
        phi    = np.pi/2. - az
        theta  = np.pi/2. - el
        adoneh5(fvis, az, 'az')
        adoneh5(fvis, el, 'el')

        ## convert the az,el to HA,DEC but assuming latitude=0
        ## if the dipole is in N-S direction, then
        ## DEC is the E-plane angle
        ## HA is the H-plane angle

        ## method 1
        phi1 = 0
        DEC = np.arcsin(np.sin(el)*np.sin(phi1)+np.cos(el)*np.cos(phi1)*np.cos(az))
        HA  = np.arcsin((-1)*np.sin(az)*np.cos(el)/np.cos(DEC))

        ## method 2
        za = np.pi/2 - el
        pntr = np.sin(za)
        pntz = np.cos(za)
        pntx = pntr * np.sin(az)
        pnty = pntr * np.cos(az)
        EWoff = np.arctan2(pntx,pntz)
        NSoff = np.arctan2(pnty,pntz)
        # overrides method 1 now
        HA = -EWoff
        DEC = NSoff

        Eatt = atten(DEC/np.pi*180., Ehphw)
        Hatt = atten(HA/np.pi*180., Hhphw)
        att0 = Eatt*Hatt
        print('max atten:', att0.max())


        unitVec  = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)], ndmin=2).T

        if (blfit == 0):    # no fitting
            xyzparams = [antXYZ] 
            fitparams = [[0,0,0,0]]
        elif (blfit == 1):  # 4 antennas in a straight line
            angrange = np.arange(-anglim, anglim, angsep) + ang0
            a1range = np.arange(-lenlim, lenlim, lensep) + pos0[0]
            a2range = np.arange(-lenlim, lenlim, lensep) + pos0[1]
            a3range = np.arange(-lenlim, lenlim, lensep) + pos0[2]
            xyzparams = []
            fitparams = []
            for ang in angrange:
                th = ang/180.*np.pi
                for a1 in a1range:
                    for a2 in a2range:
                        for a3 in a3range:
                            tmp = np.array([0., a1, a2, a3])
                            tmpx = tmp * np.cos(th)
                            tmpy = tmp * np.sin(th)
                            tmpz = np.zeros_like(tmp)
                            xyzparams.append(np.array([tmpx, tmpy, tmpz]).T)
                            fitparams.append([ang, a1, a2, a3])
            xyzparams = np.array(xyzparams)
        nParams = len(xyzparams)
        print('fitting xyz, nParams =', nParams)


        xyzlike = []
        xyzmodComplex = []
        for fiti in range(nParams):
            xyz = xyzparams[fiti]
            tmpmodComplex = []
            b = -1
            like = 0.
            for ai in range(nInp-1):
                for aj in range(ai+1, nInp):
                    b += 1
                    BVec = xyz[aj] - xyz[ai]
                    
                    Bproj  = np.zeros_like(theta)
                    for j in range(len(unitVec)):
                        Bproj[j]  = np.dot(BVec, unitVec[j])  # in meters
                    modComplex = np.exp(-2.j*np.pi*Bproj/lam)
                    tmpmodComplex.append(modComplex)
                    resPhase = np.angle(fa_spec2[:,b]/modComplex)
                    resPhase -= np.median(resPhase)
                    like += -0.5*(resPhase**2/0.3**2).sum()  # fixed phase error
            xyzlike.append(like)
            xyzmodComplex.append(tmpmodComplex)

        xyzlike = np.array(xyzlike)
        maxidx = xyzlike.argmax()
        maxxyz = xyzparams[maxidx]
        print('best xyz:', maxxyz)
        print('best fit param:', fitparams[maxidx])
        blmodComplex = xyzmodComplex[maxidx]
        blmodPhase = np.angle(blmodComplex)
        adoneh5(fvis, maxxyz, 'ant_xyz')



    ## plotting for each baseline
    savTsys = []
    savTsysmask = []
    savTsys2D = []
    savTsys2Dmask = []
    savSEFD2D = []
    savSEFD2Dmask = []
    for b, pa in enumerate(pairs):
        ai, aj = pa
        corr = '%d-%d'%(ai,aj)

        if (do_model):
            modComplex = blmodComplex[b]
            modPhase = blmodPhase[b]

        ## vs. IF plot
        pwd  = os.getcwd()
        pwd1 = os.path.basename(pwd)
        for pp in range(2):
            if (pp == 0):
                spec = spec1[:,b]
                coeff = coeff1[:,b]
                #png1 = '%s.bl%s.all_win.png' % (fvis, corr)
                png1 = '%s/bl%s.all_win.png' % (odir, corr)
                suptitle = 'file: %s, bl:%s' % (fvis, corr)
            else:
                spec = spec2[:,b]
                coeff = coeff2[:,b]
                #png1 = '%s.bl%s.all_win.sub.png' % (fvis, corr)
                png1 = '%s/bl%s.all_win.sub.png' % (odir, corr)
                suptitle = 'file: %s, bl:%s, sub_t-mean' % (fvis, corr)

            # spectral avg
            fr = np.ma.median(coeff.real, axis=1)
            fi = np.ma.median(coeff.imag, axis=1)
            fa_coeff = fr + 1.j*fi
            if (do_model):
                # fring-stopping
                fs_coeff = fa_coeff / modComplex
                cr = np.ma.median(fs_coeff.real)
                ci = np.ma.median(fs_coeff.imag)
                ap_cal = cr + 1.j*ci    # cal both amp and phase
                p_cal = ap_cal / np.abs(ap_cal) # cal only phase
                calPhase = np.ma.angle(p_cal)
                # median phase calibrate
                cfs_coeff = fs_coeff / p_cal
                modPhase2 = np.angle(modComplex * p_cal)


            #f1, s1 = plt.subplots(2,2,figsize=(12,9),sharex=True)
            f1, s1 = plt.subplots(3,2,figsize=(12,12))

            ax = s1[1,0]
            ax.set_title('cross01')
            ax.set_xlabel('freq (MHz)')
            #ax.set_ylabel('vis.power (dBm)')
            ax.set_ylabel('xcorr coeff')
            ax.set_yscale('log')
            ax.set_ylim([2e-3, 2])
            ax = s1[1,1]
            ax.set_title('cross01')
            ax.set_xlabel('freq (MHz)')
            ax.set_ylabel('vis.phase (rad)')
            ax.set_ylim([-3.2,3.2])
            for j in range(2):
                ax = s1[0,j]
                ax.set_title('ch%d'%j)
                ax.set_xlabel('freq (MHz)')
                ax.set_ylabel('power (dBm)')


            for i in range(nWin):
                #print('   window: (%d/%d)'%(i+1,nWin))
                cc = 'C%d'%(i%10)

                ax = s1[1,0]
                #ax.plot(freq, toPower(np.ma.abs(spec[i]), mode='vis'), color=cc)
                ax.plot(freq, np.ma.abs(coeff[i]), color=cc)

                ax = s1[1,1]
                ax.plot(freq, np.ma.angle(spec[i]), color=cc, linestyle='none', marker='.')

                for j in range(2):
                    ax = s1[0,j]
                    if (j==1):
                        lab = '%04.0fsec'%sec[i]
                    else:
                        lab = ''
                    ax.plot(freq, toPower(nft[i,j]), color=cc, label=lab)

            med_coeff = np.ma.median(np.ma.abs(coeff), axis=1)
            par2 = np.percentile(med_coeff, [15,50,85])
            med_med_coeff = par2[1]
            std_med_coeff = (par2[2]-par2[0])/2.
            ax = s1[1,0]

            ax = s1[2,0]
            if (do_model):
                ax.plot(sec, sigma_clip(cfs_coeff.real), 'b-', label='real')
            ax.plot(sec, sigma_clip(med_coeff), 'b:', label='abs')
            ax.legend()
            ax.set_xlabel('window time (sec)')
            ax.set_ylabel('med_coefficient')
            ax.text(0.05, 0.90, 'med:%.3g+/-%.3g'%(med_med_coeff, std_med_coeff), transform=ax.transAxes)

            med_phase = np.ma.median(np.ma.angle(spec), axis=1)
            ax = s1[2,1]
            ax.plot(sec, med_phase)
            if (do_model):
                ax.plot(sec, modPhase2)
            ax.set_ylim([-3.2,3.2])
            ax.set_xlabel('window time (sec)')
            ax.set_ylabel('med_phase')


            #s1[0,1].legend()
            #s1[0,0].get_shared_y_axes().join(s1[0,0], s1[0,1], s1[1,0])
            s1[0,0].get_shared_y_axes().join(s1[0,0], s1[0,1])
            s1[0,0].get_shared_x_axes().join(s1[0,0], s1[0,1], s1[1,0], s1[1,1])
            s1[2,0].get_shared_x_axes().join(s1[2,0], s1[2,1])

            f1.tight_layout(rect=[0,0.03,1,0.95])
            f1.suptitle(suptitle)
            f1.savefig(png1)
            plt.close(f1)


            ## vs. time plot
            if (pp == 0):
                #png2 = '%s.bl%s.vs_time.png' % (fvis, corr)
                png2 = '%s/bl%s.vs_time.png' % (odir, corr)
                suptitle = 'file: %s, bl:%s\nIF=[%.1f,%.1f]'%(fvis, corr, ifmin, ifmax)
            else:
                #png2 = '%s.bl%s.vs_time.sub.png' % (fvis, corr)
                png2 = '%s/bl%s.vs_time.sub.png' % (odir, corr)
                suptitle = 'file: %s, bl:%s, sub_t-mean\nIF=[%.1f,%.1f]'%(fvis, corr, ifmin, ifmax)

            f2, s2 = plt.subplots(3,1, figsize=(8,9), sharex=True)

            # freq avg
            faspec = spec[:,fsel].mean(axis=1)
            faphase = np.ma.median(np.ma.angle(spec), axis=1)

            # vs sec or vs date
            #tt = sec
            tt = tWin
            
            # phase plot
            ax = s2[0]
            #ax.plot(tt, np.ma.angle(faspec), marker='.')
            ax.plot(tt, faphase, marker='.')
            if (do_model):
                ax.plot(tt, modPhase2)
            ax.set_ylim([-3.3, 3.3])
            ax.set_xlabel('time')
            ax.set_ylabel('vis.phase (rad)')

            # power plot
            ax = s2[1]
            ax.set_xlabel('time')
            sel = 1
            if (sel==0):
                if (do_model):
                    ax.plot(tt, sigma_clip(cfs_coeff.real), 'b-', label='real')
                ax.plot(tt, sigma_clip(med_coeff), 'b:', label='abs')
                ax.legend()
                ax.set_ylabel('corr strength')
            elif (sel==1):
                #if (do_model):
                    #Tsys1 = srcflux / cfs_coeff.real * scale2
                    #ax.plot(tt, sigma_clip(Tsys1), 'b-', label='from real')
                Tsys2 = srcflux / med_coeff * scale2
                # correct for Sun contribution to noise
                if (do_model):
                    Tsys2 *= att0
                Tsys2 *= (1-med_coeff)
                cTsys2 = sigma_clip(Tsys2)
                if (pp==1):
                    savTsys.append(cTsys2)
                    savTsysmask.append(cTsys2.mask)
                ax.plot(tt, cTsys2, 'b:', label='from abs')
                ax.text(0.05, 0.60, 'Tsys_min=%.0fK'%cTsys2.min(), transform=ax.transAxes)
                ax.text(0.05, 0.53, 'Tsys_med=%.0fK'%np.median(cTsys2), transform=ax.transAxes)
                ax.axhline(np.median(cTsys2), linestyle='--', color='k')
                ax.legend()
                ax.set_ylabel('Tsys (K)')

            # elevation plot
            ax = s2[2]
            if (do_model):
                ax.plot(tt, el/np.pi*180, marker='.')
            ax.set_xlabel('time')
            ax.set_ylabel('elevation (deg)')

            f2.tight_layout(rect=[0,0.03,1,0.95])
            f2.suptitle(suptitle)
            f2.savefig(png2)
            plt.close(f2)

        ## 2D Tsys, SEFD plot
        if (do_model):
            Tsys2D = srcflux2.reshape((1,-1)) / np.ma.abs(coeff) * fscale2.reshape((1,-1))
            Tsys2D *= att0.reshape((-1,1))
            # correct Tsys bias due to Solar flux
            Tsys2D *= (1.-np.ma.abs(coeff))
            # Tsys.shape = (nWin, nChan)
            savTsys2D.append(Tsys2D)
            savTsys2Dmask.append(Tsys2D.mask)
            Tsys_time = np.ma.median(Tsys2D, axis=1)
            Tsys_freq = np.ma.median(Tsys2D, axis=0)
            Tsys_time2 = Tsys2D.min(axis=1)
            Tsys_freq2 = Tsys2D.min(axis=0)
            #vmin = Tsys2D.min()
            #vmax = Tsys2D.max()
            vmin = 0
            vmax = 300

            SEFD2D = Tsys2D / fscale2.reshape((1,-1))
            SEFD2D *= (flam/flam0)**2
            savSEFD2D.append(SEFD2D)
            savSEFD2Dmask.append(SEFD2D.mask)
            SEFD2D_time = np.ma.median(SEFD2D, axis=1)
            SEFD2D_freq = np.ma.median(SEFD2D, axis=0)
            SEFD2D_time2 = SEFD2D.min(axis=1)
            SEFD2D_freq2 = SEFD2D.min(axis=0)

            fig = plt.figure(figsize=(12,9))
            gs = fig.add_gridspec(2,2,width_ratios=(3,1),height_ratios=(1,3), wspace=0, hspace=0, left=0.08, right=0.95, bottom=0.08, top=0.92)
            ax2D = fig.add_subplot(gs[1,0])
            axtp = fig.add_subplot(gs[0,0])
            axrt = fig.add_subplot(gs[1,1])

            s = ax2D.pcolormesh(T,Y,Tsys2D.T, vmin=vmin, vmax=vmax)
            ax2D.set_xlabel('time (UT)')
            ax2D.set_ylabel('freq (MHz)')
            #cb = plt.colorbar(s, ax=ax2D)
            #cb.set_label('Tsys (K)')

            axtp.plot(sec, Tsys_time, label='med')
            axtp.plot(sec, Tsys_time2, linestyle=':', label='min')
            axtp.legend()
            axtp.set_ylabel('Tsys (K)')
            axtp.set_ylim(vmin, vmax)
            axtp.grid(axis='both')
            axtp.set_xticks([])

            if (np.count_nonzero(Tsys_freq.mask) < nChan):
                axrt.plot(Tsys_freq, freq, label='med')
            if (np.count_nonzero(Tsys_freq2.mask) < nChan):
                axrt.plot(Tsys_freq2, freq, linestyle=':', label='min')
            axrt.legend()
            axrt.set_xlabel('Tsys (K)')
            axrt.set_xlim(vmin, vmax)
            axtp.grid(axis='both')
            axrt.set_yticks([])

            fig.text(0.76, 0.80, 'med(Tsys_freq.med) = %.0f K'%np.ma.median(Tsys_freq))
            fig.text(0.76, 0.77, 'med(SEFD_freq.med) = %.2f MJy'%(np.ma.median(SEFD2D_freq)/1e6))
            fig.suptitle(suptitle)
            png4 = '%s/bl%s.rf%d_%d.Tsys2D.png' % (odir, corr, rfmin, rfmax)
            fig.savefig(png4)
            plt.close(fig)


    if (do_model):
        savTsys = np.ma.array(savTsys, mask=savTsysmask)
        adoneh5(fvis, savTsys.T, 'winTsys')   # shape (nWin, nBl), after the transpose
        savTsys2D = np.ma.array(savTsys2D, mask=savTsys2Dmask)
        adoneh5(fvis, savTsys2D.transpose((1,0,2)), 'winTsys2D') # shape (nWin, nBl, nChan)
        savSEFD2D = np.ma.array(savSEFD2D, mask=savSEFD2Dmask)
        adoneh5(fvis, savSEFD2D.transpose((1,0,2)), 'winSEFD2D') # shape (nWin, nBl, nChan)



