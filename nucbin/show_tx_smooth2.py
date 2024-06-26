#!/usr/bin/env python

from read_cli import *
import matplotlib.pyplot as plt
import scipy.signal as sg
import sys, os.path
from subprocess import call


def load_csv(fname, samprate, nsamp):
    csv = np.loadtxt(fname, delimiter=',', max_rows=nsamp, unpack=True)
    data0 = csv[0] + 1j*csv[1]
    data1 = csv[2] + 1j*csv[3]
    #data0 = 1j*csv[0] + csv[1]
    #data1 = 1j*csv[2] + csv[3]
    time1 = np.arange(len(data0), dtype=float) / samprate
    return data0, data1, time1


winlen    = 50     # num of samples to avg
odir      = 'figs'
do_smooth = True
do_cross  = True
use_csv   = False
nseg      = 3000     # num of seg to read, each seg is 4096 bytes (1024 samples, for two channels)
samprate  = 20e6    # Hz


inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
plot a few diagnostics:
    - smoothed data showing the TX signal in the beginning
    - averaged spectrum
    - auto- and cross- correlation of data with +/-2 sample delay

%s <files> [options]

options are:
-o DIR      save outputs in the DIR
--avg NUM   set the number of samples to smooth (%d)

''' % (pg, winlen)


if (len(inp) < 1):
    sys.exit(usage)

files = []
while(inp):
    k = inp.pop(0)
    if (k == '-o'):
        odir = inp.pop(0)
    elif (k == '--avg'):
        winlen = int(inp.pop(0))
    else:
        files.append(k)

if (not odir == '.'):
    call('mkdir -p %s' % odir, shell=True)

win = np.zeros(winlen*3, dtype=float)
win[winlen:winlen*2] = 1
win /= winlen


data = []
for fname in files:
    #fname = 'cli1-rxtx-mimo.bin'
    if (not os.path.isfile(fname)):
        print('skip. file not found:', fname)
        next

    print('plotting', fname)
    nsamp = nseg*512    # 512 is the num of samples per 4096bytes seg in dual channel mode
    if (fname.endswith('csv')):
        use_csv = True
        data0, data1, time1 = load_csv(fname, samprate, nsamp)
    else:
        data0, data1, time1 = load_cli_data(fname, samprate, nseg=nseg, nskip=0, dual=True, verbose=0)
        data0 = data0.flatten()
        data1 = data1.flatten()
    data.extend([data0, data1])

    if (do_smooth):
        fig, ax = plt.subplots()
        #x = time1[0:]
        x = range(time1.size)
        xmax = 30*512
        smooth0 = sg.convolve(data0, win, mode='same')
        smooth1 = sg.convolve(data1, win, mode='same')
        ax.plot(x[:xmax], smooth0.real[:xmax], label='ch0.real')
        ax.plot(x[:xmax], smooth0.imag[:xmax], label='ch0.imag')
        ax.plot(x[:xmax], smooth1.real[:xmax], label='ch1.real')
        ax.plot(x[:xmax], smooth1.imag[:xmax], label='ch1.imag')
        fig.savefig('%s/%s.smooth.png' % (odir, fname))
        plt.close(fig)


nInp = len(data)
if (do_cross):
    nplt = 5
    nchan = 1024
    nsamp = len(time1)
    nspec = nsamp // nchan
    for chi in range(nInp):
        tmp0 = data[chi][2:(nspec-1)*nchan+2]
        ft0 = np.fft.fft(tmp0.reshape((-1,nchan)), axis=1)
        ft0 = np.fft.fftshift(ft0, axes=1)

        fid = chi // 2
        fname = files[fid]
        png1 = '%s/%s.spec.png' % (odir, fname)
        f1, sub1 = plt.subplots(2,2,figsize=(8,6))
        ax = sub1[0,0]
        ax.plot(range(nchan), np.abs(ft0).mean(axis=0))
        ax.set_xlabel('chan_idx')
        ax.set_ylabel('abs(data)')
        ax.set_title('ch0')
        ax = sub1[0,1]
        ax.plot(range(nchan), np.angle(ft0).mean(axis=0))
        ax.set_xlabel('chan_idx')
        ax.set_ylabel('angle(data)')
        for chj in range(chi, nInp):
            print('inputs:', chi, chj)
            png = '%s/%s.cc%d%d.png' % (odir, fname, chi, chj)
            fig, sub = plt.subplots(nplt,4,figsize=(12,nplt*2.0), squeeze=False)
            for i in range(-2,3):
                iy = i+2
                tmp1 = data[chj][2+i:(nspec-1)*nchan+2+i]
                ft1 = np.fft.fft(tmp1.reshape((-1,nchan)), axis=1)
                ft1 = np.fft.fftshift(ft1, axes=1)
                cross = ft0 * ft1.conjugate()

                if (i == 0 and chj==chi+1):
                    ax = sub1[1,0]
                    ax.plot(range(nchan), np.abs(ft1).mean(axis=0))
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('abs(data)')
                    ax.set_title('ch1')
                    ax = sub1[1,1]
                    ax.plot(range(nchan), np.angle(ft1).mean(axis=0))
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('angle(data)')

                    f1.tight_layout()
                    f1.savefig(png1)
                    plt.close(f1)

                    f2, sub2 = plt.subplots(2,2,figsize=(8,6))
                    png2 = '%s/%s.cross.png' % (odir, fname)
                    #scross = sg.convolve(cross, win.reshape((1,-1)), mode='same')  # spectral smoothing
                    scross = sg.convolve(cross, win.reshape((-1,1)), mode='same')   # temporal smoothing
                    ax = sub2[0,0]
                    for j in range(10):
                        ax.plot(range(1024), np.abs(scross[j*50,:]), color='C%d'%j)
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('abs(cross)')
                    ax = sub2[0,1]
                    for j in range(10):
                        ax.plot(range(1024), np.angle(scross[j*50,:]), color='C%d'%j)
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('angle(cross)')
                    ax = sub2[1,0]
                    for j in range(10):
                        ax.plot(range(1024), scross.real[j*50,:], color='C%d'%j)
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('real(cross)')
                    ax = sub2[1,1]
                    for j in range(10):
                        ax.plot(range(1024), scross.imag[j*50,:], color='C%d'%j)
                    ax.set_xlabel('chan_idx')
                    ax.set_ylabel('imag(cross)')

                    f2.tight_layout()
                    f2.savefig(png2)
                    plt.close(f2)

                ax = sub[iy, 0]
                s = ax.imshow(np.abs(cross).T, origin='lower', aspect='auto')
                cb = plt.colorbar(s, ax=ax)
                cb.set_label('abs(cross)')
                ax.set_xlabel('time_idx')
                ax.set_ylabel('chan_idx')

                ax = sub[iy,1]
                #ax.plot(range(nchan), np.abs(cross.mean(axis=0)))
                ax.plot(range(nchan), np.abs(cross).mean(axis=0))
                ax.set_xlabel('chan_idx')
                ax.set_ylabel('abs(cross)')

                ax = sub[iy, 2]
                s = ax.imshow(np.angle(cross).T, origin='lower', aspect='auto')
                cb = plt.colorbar(s, ax=ax)
                cb.set_label('angle(cross)')
                ax.set_xlabel('time_idx')
                ax.set_ylabel('chan_idx')

                ax = sub[iy,3]
                #ax.plot(range(nchan), np.angle(cross.mean(axis=0)))
                ax.plot(range(nchan), np.angle(cross).mean(axis=0))
                ax.set_xlabel('chan_idx')
                ax.set_ylabel('angle(cross)')

            fig.tight_layout()
            fig.savefig(png)
            plt.close(fig)

