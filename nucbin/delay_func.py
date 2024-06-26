#!/usr/bin/env python

from loadh5 import *


## the antenna-delay functions ##
def lnlike(params, freq, y, yerr):
    # data:
    # spec = complex spectra, shape [nsb, nb, nch]
    # svar = float variance,  shape [nsb, nb, nch]
    # data may have been smoothed in spactral domain

    ymod = phimodel(params, freq)

    # log likelihood
    like = -0.5*np.sum(pwrap2(y-ymod)**2/yerr**2)
    return like


def phimodel(params, freq):
    # model:
    # params = antenna-delays (sec) and phase offsets (rad)
    delays = params[:na]
    offset = params[na:]

    # data:
    # freq = RF frequency in Hz, shape [nsb, nch]
    s    = 2. * np.pi * (freq * 1.) * 1.      # phi = s * delay_ps

    ymod = np.zeros((nsb, nb, nch))
    b = -1
    for ai in range(na-1):
        for aj in range(ai+1, na):
            b += 1

            dtau = delays[aj] - delays[ai]
            phi0 = offset[aj] - offset[ai]
            phi = s * dtau + phi0
            
            ymod[:,b,:] = phi
    
    return ymod


def lnprob(params, ranages, freq, y, yerr):
    return lnflat(params, ranges) + lnlike(params, freq, y, yerr)


## the baselin-delay functions ##
def lnlike2(params, freq, y, yerr):
    # data: phase, phase_err
    # spec = complex spectra, shape [nsb, nch]
    # svar = float variance,  shape [nsb, nch]
    # data may have been smoothed in spactral domain

    ymod = phimodel2(params, freq)

    # log likelihood
    like = -0.5*np.ma.sum(pwrap2(y-ymod)**2/yerr**2)
    return like


def phimodel2(params, freq):
    # model:
    # params = antenna-delays (sec) and phase offsets (rad)
    delays = params[0]
    offset = params[1]

    # data:
    # freq = RF frequency in GHz, shape [nsb, nch]
    #s    = 2. * np.pi * (freq * 1.e9) * 1.e-12      # phi = s * delay_ps
    # make it more general, freq in Hz, delay in sec
    s    = 2. * np.pi * (freq * 1.) * 1.      # phi = s * delay

    ## the output (ymod) will match the shape of freq
    # if freq is flattened, then ymod is flattened as well
    #ymod = np.zeros((nsb, nch))
    ymod = np.zeros_like(freq)

    dtau = delays
    phi0 = offset
    phi = s * dtau + phi0
    
    ymod = phi
    
    return ymod


def fitfun(xdata, a, b):
    params = [a, b]
    return phimodel2(params, xdata)


def lnprob2(params, ranages, freq, y, yerr):
    return lnflat(params, ranges) + lnlike2(params, freq, y, yerr)


def scanLH(ranges, freq, y, yerr):
    ngrid = 30

    npar = len(ranges)
    #print('npar:', npar, 'ranges:', ranges)
    parlist = np.zeros((npar, ngrid))

    for i in range(npar):
        x = np.linspace(ranges[i,0], ranges[i,1], ngrid+1)
        parlist[i] = x[:ngrid]

    ### for 2 parameters only (npar == 2)
    LH = np.zeros((ngrid, ngrid))

    for i in range(ngrid):
        for j in range(ngrid):
            par = [parlist[0,i], parlist[1,j]]
            LH[i,j] = lnlike2(par, freq, y, yerr)

    maxLH = np.ma.max(LH)
    maxPos = np.unravel_index(np.ma.argmax(LH), LH.shape)
    maxpar = [parlist[0,maxPos[0]], parlist[1,maxPos[1]]]

    return LH, maxLH, maxpar


## general  function ##
def lnflat(params, ranges):
    # flat prior
    # input:
    #   params (N)
    #   ranges (N, 2) --> the lower and upper limits for each param
    # return log probability
    prob = 0.
    for p, r in zip(params, ranges):
        if (p>r.max() or p<=r.min()):
            prob = -np.inf
    return prob


def pwrap2(phi, lim=np.pi):
    '''
    wrap input phase between +/- lim

    input:
        phi     array of phase, any shape
        lim     e.g. np.pi or 180 (deg)
    output:
        phi     same shape as the input
    '''
    phi2 = phi.copy()
    while(np.any(phi2 >= lim)):
        phi2[phi2 >= lim] -= 2.*lim
    while(np.any(phi2 < -lim)):
        phi2[phi2 < -lim] += 2.*lim

    return phi2

