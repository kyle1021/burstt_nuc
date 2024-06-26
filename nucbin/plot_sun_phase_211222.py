#!/usr/bin/env python

from pyplanet import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import timedelta
from read_meta3 import *


useData = False
#samprate = 5e6      # sampling rate in Hz
samprate = 10e6      # sampling rate in Hz  --> test, fixing slow fringe rate
seg_byte = 8000000  # taroge data config
size_byte = 4       # taroge data config
nchan = 1           # flatten the input data
nchan2 = 1024       # FFT channel number
tdur  = 1.          # length of data to read, in sec
tpseg = seg_byte/(size_byte*2)/samprate # duration per segment
ratt = 0.5  # reflected signal attanuation factor

site = 'taroge4'
t0 = datetime(2021,12,22,22,20,0)
b, obs = obsBody('sun', time=t0, site=site, retOBS=True)
print('== starting point ==')
print('Time(UTC) = %s' % obs.date)
print('RA2000  = %s' % b.ra)
print('DEC2000 = %s' % b.dec)
print('AZ = %s' % b.az)
print('EL = %s' % b.alt)

theta = np.pi/2. - b.alt
phi = np.pi/2. - b.az
unitVec = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]

antXYZ = np.loadtxt('TAROGE4.config', usecols=(1,2,3))  # (X,Y,Z)
BVec = antXYZ[2] - antXYZ[0]
print('Ant1, Ant3, BL1-3 (m):', antXYZ[0], antXYZ[2], BVec)
Bproj = np.dot(BVec, unitVec)
print('delay = %.3f m' % Bproj)
print('===================\n')


flist = 'sun_data.UT211222'
names = np.loadtxt(flist, usecols=(0,), dtype=str)
freq, dt1 = np.loadtxt(flist, usecols=(1,2), unpack=True)
# freq --> central frequency in MHz
# dt1 --> file beginning time in seconds relative to t0
nFile = len(names)
lamb = 3.e8 / (freq * 1.e6) # wavelength in meters



if (useData):
    fig, sub = plt.subplots(4,1,figsize=(10,12), sharex=True)
else:
    fig, sub = plt.subplots(3,1,figsize=(10,9), sharex=True)

#for i in range(nFile):
for i in range(4,21):
    name = names[i]
    d1 = timedelta(seconds=dt1[i])
    t1 = t0 + d1    # starting datetime of the file
    lam = lamb[i]

    dates = []
    az = []
    el = []
    obsPha = []
    obsAmp = []
    # calculate for every minute
    #for dmin in np.arange(0,2,1,dtype=float):
    for dmin in np.arange(0,10,1,dtype=float):
        di = timedelta(minutes=dmin)
        ti = t1 + di
        dates.append(ti)

        obs.date = ti
        b.compute(obs)

        az.append(b.az)
        el.append(b.alt)

        if (useData):
            tskip = dmin * 60.                          # seconds to skip
            nskip = np.floor(tskip/tpseg).astype(int)   # num of segments to skip
            nseg  = np.ceil(tdur/tpseg).astype(int)     # num of segments to read
            #print('nskip, nseg:', nskip, nseg)
            print('reading', name, 'nseg=', nseg, 'nskip=', nskip, 'tskip=%.1f'%tskip)
            fname0 = '../4TB/%s-2ant-ch0.meta' % name
            data0, timestamp = load_blade_data(fname0, nskip=nskip, nseg=nseg, nchan=nchan, size_byte=size_byte, seg_byte=seg_byte, meta=True, vector=False, verbose=0)
            fname1 = '../4TB/%s-2ant-ch1.meta' % name
            data1, timestamp = load_blade_data(fname1, nskip=nskip, nseg=nseg, nchan=nchan, size_byte=size_byte, seg_byte=seg_byte, meta=True, vector=False, verbose=0)
            dual = [data0, data1]
            for ch in range(2):
                data = dual[ch].flatten()
                nspec = len(data)//nchan2
                data = data[:nspec*nchan2].reshape((-1,nchan2))
                data = np.fft.fft(data, axis=1)
                data = np.fft.fftshift(data, axes=1)
                dual[ch] = data
            cross = dual[0] * dual[1].conjugate()
            #cross /= (np.abs(dual[0])*np.abs(dual[1]))
            crossSpec = cross.mean(axis=0)
            crossInt  = crossSpec.sum()
            obsPha.append(np.angle(crossInt))
            if (np.abs(crossInt) < 500):
                obsAmp.append(np.abs(crossInt))
            else:
                obsAmp.append(np.nan)
        else:
            obsPha.append(np.nan)
            obsAmp.append(np.nan)
        

    dates = np.array(dates)
    az = np.array(az)
    el = np.array(el)
    obsPha = np.array(obsPha)
    print(name, 'el:', el[0]/np.pi*180, el[-1]/np.pi*180)


    phi    = np.pi/2. - az
    theta  = np.pi/2. - el
    theta2 = np.pi/2. + el  # reflection off the see
    unitVec  = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]).T
    unitVec2 = np.array([np.sin(theta2)*np.cos(phi), np.sin(theta2)*np.sin(phi), np.cos(theta2)]).T
    Bproj  = np.zeros_like(theta)
    Bproj2 = np.zeros_like(theta2)
    for j in range(len(unitVec)):
        Bproj[j]  = np.dot(BVec, unitVec[j])  # in meters
        Bproj2[j] = np.dot(BVec, unitVec2[j]) # in meters
    phase  = np.angle(np.exp(2.j*np.pi*Bproj/lam))
    phase2 = np.angle(np.exp(2.j*np.pi*Bproj/lam) + ratt*np.exp(2.j*np.pi*Bproj2/lam))    # reflection attenuation: ratt
    phase  *= -1.
    phase2 *= -1.


    # az,el plot
    ax = sub[0]
    l1, = ax.plot(dates, az/np.pi*180., 'b-', label='az')
    l2, = ax.plot(dates, el/np.pi*180., 'g-', label='el')
    if (i==0):
        ax.legend(handles=[l1,l2])
    ax.set_ylabel('AZ,EL (deg)')
    ax.set_title('Sun trajectory')
    ax.set_aspect('auto')
    # delay plot
    ax = sub[1]
    l1, = ax.plot(dates, Bproj,  'b-', label='direct')
    l2, = ax.plot(dates, Bproj2, 'g-', label='reflect')
    if (i==0):
        ax.legend(handles=[l1,l2])
    ax.set_ylabel('delay1-3 (m)')
    ax.set_title('delay between Ant1-Ant3')
    ax.set_aspect('auto')

    cc = 'C%d'%(i%10)
    if (useData):
        ax = sub[3]
        ax.scatter(dates, obsAmp, color=cc)
        ax.set_ylabel('vis.ampld')
        ax.set_aspect('auto')

    # variable freq phase plot
    ax = sub[2]
    #ax.plot(dates, phase1var, label='direct')
    fMHz = freq[i]
    ax.plot(dates, phase, color=cc, label='%.0fMHz'%fMHz)
    ax.plot(dates, phase2, color=cc, linestyle=':')
    ax.scatter(dates, obsPha, color=cc)
    ax.legend()
    ax.set_ylabel('vis.phase (rad)')
    ax.set_title('solid: direct phase. dotted: direct + %.2f*reflect' % ratt)
    ax.set_aspect('auto')

    sub[-1].set_xlabel('Time (UTC)')


fig.autofmt_xdate()
fig.tight_layout(rect=[0,0.03,1,0.95])
fig.subplots_adjust(hspace=0.1)
fig.suptitle('Sun trajectory, phase for Baseline 1-3, 2021/12/22')
fig.savefig('sun_trajectory_211222.png')
plt.close(fig)


