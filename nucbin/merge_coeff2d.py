#!/usr/bin/env python

from loadh5 import *
import matplotlib.pyplot as plt
from astropy.time import Time

freq_list = [300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740, 780]
#freq_list = [300, 340, 460, 500, 540, 660, 700]
nFreq = len(freq_list)

png = 'all_freq_coeff2d.png'
fig, sub = plt.subplots(2,2,figsize=(15,12), sharex=True, sharey=True)
for ii in range(2):
    sub[ii,0].set_ylabel('freq (MHz)')
    sub[1,ii].set_xlabel('time (UT)')

sub[0,0].set_title('real')
sub[0,1].set_title('imag')
sub[1,0].set_title('abs()')
sub[1,1].set_title('angle()')


for fi in range(nFreq):
    fcen = freq_list[fi]
    fvis = 'analyze_%.0f.vish5' % fcen

    print('plotting ', fvis, ' ...')

    attrs = getAttrs(fvis)
    epoch0 = attrs['unix_utc_open']

    coeff  = getData(fvis, 'winCoeff')
    if (len(coeff.shape) == 2):     # backward compatibility for single_dev_vis.py
        nTime, nChan = coeff.shape
        coeff = coeff.reshape((nTime, 1, nChan))
    coeff2 = getData(fvis, 'winCoeffSub')
    freq = getData(fvis, 'freq')
    tsec = getData(fvis, 'winSec')
    tWin = Time(epoch0+tsec, format='unix').to_datetime()

    #X, Y = np.meshgrid(tsec, freq, indexing='xy')
    X, Y = np.meshgrid(tWin, freq, indexing='xy')

    sub[0,0].pcolormesh(X,Y,coeff2.real[:,0,:].T)
    sub[0,1].pcolormesh(X,Y,coeff2.real[:,0,:].T)
    sub[1,0].pcolormesh(X,Y,np.ma.abs(coeff2)[:,0,:].T)
    sub[1,1].pcolormesh(X,Y,np.ma.angle(coeff2)[:,0,:].T)

    print('coeff min/max:', np.ma.abs(coeff2.min()), np.ma.abs(coeff2.max()))


fig.autofmt_xdate()
fig.tight_layout(rect=[0,0.03,1,0.95])
fig.subplots_adjust(wspace=0, hspace=0.06)
fig.savefig(png)
plt.close(fig)


