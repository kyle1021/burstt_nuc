#!/usr/bin/env python

from read_cli2 import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from glob import glob


## spectrum analyzer data, avg-100
fSA = 'Trace_1way-39k.csv'

SAdata = np.loadtxt(fSA, skiprows=44, delimiter=',', unpack=True)
# SAdata[0]: freq in Hz
# SAdata[1]: power in dB
#print(SAdata.shape)
#print(SAdata[:,:10])
funSA = interp1d(SAdata[0], SAdata[1], fill_value='extrapolate')  # interpolation function


dev = 'dev3'
files = glob('1to1/NS20dB*%s.log'%dev)
files.sort()

nchan = 1024
rate = 40e6
Ifreq = np.fft.fftfreq(nchan, d=1/rate) # Hz
Ifreq = np.fft.fftshift(Ifreq)

f1, s1 = plt.subplots(2,1,figsize=(8,8), sharey=True)

for flog in files:
    attrs = logParse(flog)
    fcen = attrs['freq'][0]             # Hz
    #rate = attrs['rate'][0]

    fmeta = flog.replace('.log', '.meta')
    data, clock = load_cli_data2(fmeta, nseg=1000, nskip=0)
    ftdata = maskedFFT(data[:,0], 1024)
    pow0 = 20*np.ma.log10(np.ma.abs(ftdata).mean(axis=0))
    ftdata = maskedFFT(data[:,1], 1024)
    pow1 = 20*np.ma.log10(np.ma.abs(ftdata).mean(axis=0))

    Rfreq = Ifreq + fcen                # Hz

    refpow = funSA(Rfreq)
    bp0 = pow0 - refpow
    bp1 = pow1 - refpow

    s1[0].plot(Ifreq[5:nchan-5]/1e6, bp0[5:nchan-5], label='%.0fMHz'%(fcen/1e6))
    s1[1].plot(Ifreq[5:nchan-5]/1e6, bp1[5:nchan-5], label='%.0fMHz'%(fcen/1e6))

for i in range(2):
    ax = s1[i]
    ax.legend()
    ax.set_xlabel('IF (MHz)')
    ax.set_ylabel('BP (dB)')
    ax.set_title('%s-ch%d'%(dev,i))
    ax.set_xlim([-22, 32])

f1.tight_layout(rect=[0,0.03,1,0.95])
f1.suptitle('1to1, blade-gain:20dB')
f1.savefig('bp20dB_%s.png' % dev)
plt.close(f1)


