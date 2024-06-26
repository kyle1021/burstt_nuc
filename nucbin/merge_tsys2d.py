#!/usr/bin/env python

from loadh5 import *
import matplotlib.pyplot as plt
from astropy.time import Time

freq_list = [300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740, 780]
#freq_list = [300, 340, 460, 500, 540, 660, 700]
nFreq = len(freq_list)

zmin = 0
zmax = 150
tmin = 10800
tmax = 21600    
rfmin = 0
rfmax = 0

png = 'all_freq_tsys2d.png'
fig, sub = plt.subplots(2,2,figsize=(15,12), height_ratios=[1,3], width_ratios=[3,1])
sub[1,0].set_ylabel('freq (MHz)')
sub[1,0].set_xlabel('time (UT)')
sub[0,0].set_ylabel('Tsys (K)')
#sub[0,0].set_xlabel([])
sub[0,0].grid()
sub[1,1].set_xlabel('Tsys (K)')
#sub[1,1].set_yticks([])
sub[1,1].grid()
sub[0,1].remove()

sub[1,0].sharex(sub[0,0])
sub[1,0].sharey(sub[1,1])


tsys_all = []
tsys_freq2_all = []
for fi in range(nFreq):
    fcen = freq_list[fi]
    fvis = 'analyze_%.0f.vish5' % fcen

    print('plotting ', fvis, ' ...')

    attrs = getAttrs(fvis)
    epoch0 = attrs['unix_utc_open']

    tsys2d  = getData(fvis, 'winTsys2D')
    if (len(tsys2d.shape) == 2):     # backward compatibility for single_dev_vis.py
        nTime, nChan = tsys2d.shape
        tsys2d = tsys2d.reshape((nTime, 1, nChan))
    freq = getData(fvis, 'freq')
    tsec = getData(fvis, 'winSec')
    tWin = Time(epoch0+tsec, format='unix').to_datetime()

    tsel = np.logical_and(tsec>=tmin, tsec<=tmax)
    sel_tsys2d = tsys2d[tsel]
    if (rfmin == 0):
        rfmin = freq.min()
    else:
        rfmin = min(rfmin, freq.min())
    if (rfmax == 0):
        rfmax = freq.max()
    else:
        rfmax = max(rfmax, freq.max())

    tt1 = tWin[np.argmin(np.abs(tsec-tmin))]
    tt2 = tWin[np.argmin(np.abs(tsec-tmax))]


    #X, Y = np.meshgrid(tsec, freq, indexing='xy')
    X, Y = np.meshgrid(tWin, freq, indexing='xy')

    sub[1,0].pcolormesh(X,Y,tsys2d[:,0,:].T, vmin=zmin, vmax=zmax)
    sub[1,0].axvline(tt1, color='r', linestyle='--')
    sub[1,0].axvline(tt2, color='r', linestyle='--')

    tsys_freq  = sel_tsys2d.min(axis=(0,1))
    tsys_freq2 = np.ma.median(sel_tsys2d, axis=(0,1))
    tsys_freq2_all.append(tsys_freq2)
    if (fi==0):
        sub[1,1].plot(tsys_freq, freq, color='C1', label='min')
        sub[1,1].plot(tsys_freq2, freq, color='C0', label='med')
    else:
        sub[1,1].plot(tsys_freq, freq, color='C1')
        sub[1,1].plot(tsys_freq2, freq, color='C0')

    tsys_all.append(tsys2d)

tsys_all = np.asarray(tsys_all)
tsys_time  = tsys_all.min(axis=(0,2,3))
tsys_time2 = np.ma.median(tsys_all, axis=(0,2,3))

tsys_freq2_all = np.asarray(tsys_freq2_all)
tsys_freq2_med = np.median(tsys_freq2_all, axis=(0,1))

sub[0,0].plot(tWin, tsys_time, color='C1', label='min')
sub[0,0].plot(tWin, tsys_time2, color='C0', label='med')
sub[0,0].legend()
sub[0,0].set_ylim(zmin, zmax)
#sub[0,0].set_xlim(tWin[0], tWin[-1])

sub[1,1].legend()
print(rfmin, rfmax)
sub[1,1].set_ylim(rfmin, rfmax)
sub[1,1].set_xlim(zmin, zmax)
sub[1,1].axvline(tsys_freq2_med, color='k', linestyle='--')
fig.text(0.77, 0.80, 'med(Tsys): %.0f K'%tsys_freq2_med)

fig.tight_layout(rect=[0,0.03,1,0.95])
fig.subplots_adjust(wspace=0, hspace=0)
fig.savefig(png)
plt.close(fig)


