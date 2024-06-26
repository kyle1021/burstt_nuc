#!/usr/bin/env python

from loadh5 import *
from glob import glob
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# a single LNA connected to dev4-ch0 (RX1)
# ch1 (RX2) was open
# 4 sets of data:
#   LNAon  --> gain: 20dB, 40dB
#   LNAoff --> gain: 20dB, 40dB
#
# save to a single file

fout = 'blade_terminate.h5'
png = 'all_terminate.png'
f1, s1 = plt.subplots(1,1,figsize=(15,5))


files = glob('LNA*/*.spec.h5')
files.sort()

for i in range(4):
    f = files[i]
    if (i == 0):
        name = 'LNAoff'
        gain = 20
    elif (i == 1):
        name = 'LNAoff'
        gain = 40
    elif (i == 2):
        name = 'LNAon'
        gain = 20
    elif (i == 3):
        name = 'LNAon'
        gain = 40

    dname = '%s_%ddB' % (name, gain)

    spec = getData(f, 'specdB')
    adoneh5(fout, spec[:,0], dname)

    if (i == 0):
        freq = getData(f, 'freqMHz')
        adoneh5(fout, freq, 'freqMHz')

    cc = 'C%d'%i
    for j in range(15):
        if (j == 0):
            s1.plot(freq[j], spec[j,0], color=cc, label=dname)
        else:
            s1.plot(freq[j], spec[j,0], color=cc)

s1.set_xlabel('freq (MHz)')
s1.set_ylabel('power (dB)')
s1.legend()
s1.set_ylim([-120,-55])
s1.set_title('LNA termination results')
f1.tight_layout(rect=[0,0.03,1,0.95])
f1.savefig(png)
plt.close(f1)


