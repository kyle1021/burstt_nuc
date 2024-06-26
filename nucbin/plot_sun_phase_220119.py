#!/usr/bin/env python

from pyplanet import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import timedelta
from read_meta3 import *


ratt = 0.5  # reflected signal attanuation factor
useData = False
if (useData):
    samprate = 5e6      # sampling rate in Hz
    seg_byte = 8000000  # taroge data config
    size_byte = 4       # taroge data config
    nchan = 1           # flatten the input data
    nchan2 = 1024       # FFT channel number
    tdur  = 1.          # length of data to read, in sec
    tpseg = seg_byte/(size_byte*2)/samprate # duration per segment

site = 'taroge4'
t0 = datetime(2022,1,19,22,40,0)
b, obs = obsBody('sun', time=t0, site=site, retOBS=True)
print('== starting point ==')
print('Time(UTC) = %s' % obs.date)
print('RA2000  = %s' % b.ra)
print('DEC2000 = %s' % b.dec)
print('AZ = %s' % b.az)
print('EL = %s' % b.alt)
print('===================\n')

theta = np.pi/2. - b.alt
phi = np.pi/2. - b.az
unitVec = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]

antXYZ = np.loadtxt('TAROGE4.config', usecols=(1,2,3))  # (X,Y,Z)
#BVec = antXYZ[2] - antXYZ[0]
#print('Ant1, Ant3, BL1-3 (m):', antXYZ[0], antXYZ[2], BVec)
#Bproj = np.dot(BVec, unitVec)
#print('delay = %.3f m' % Bproj)


UTstr = t0.strftime('%y%m%d')
flist = 'sun_data.UT%s' % UTstr
names = np.loadtxt(flist, usecols=(0,), dtype=str)
freq, dt1 = np.loadtxt(flist, usecols=(1,2), unpack=True)
# freq --> central frequency in MHz
# dt1 --> file beginning time in seconds relative to t0
nFile = len(names)
lamb = 3.e8 / (freq * 1.e6) # wavelength in meters


nAnt = 4
for ai in range(nAnt-1):
    for aj in range(ai+1, nAnt):
        BLname = '%d-%d' % (ai+1, aj+1)
        BVec = antXYZ[ai] - antXYZ[aj]

        png   = 'sun_trajectory_%s_%s.png' % (UTstr, BLname)

        fig, sub = plt.subplots(3,1,figsize=(10,9), sharex=True)

        #for i in range(2):
        for i in range(nFile):
            name = names[i]
            d1 = timedelta(seconds=dt1[i])
            t1 = t0 + d1    # starting datetime of the file
            lam = lamb[i]

            dates = []
            az = []
            el = []
            obsPha = []
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
                    cross /= (np.abs(dual[0])*np.abs(dual[1]))
                    crossSpec = cross.mean(axis=0)
                    crossInt  = crossSpec.sum()
                    obsPha.append(np.angle(crossInt))
                else:
                    obsPha.append(np.nan)
                

            dates = np.array(dates)
            az = np.array(az)
            el = np.array(el)
            obsPha = np.array(obsPha)


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


            # az,el plot
            ax = sub[0]
            l1, = ax.plot(dates, az/np.pi*180., 'b-', label='az')
            l2, = ax.plot(dates, el/np.pi*180., 'g-', label='el')
            if (i==0):
                ax.legend(handles=[l1,l2])
            ax.set_ylabel('AZ,EL (deg)')
            ax.set_title('Sun trajectory')
            # delay plot
            ax = sub[1]
            l1, = ax.plot(dates, Bproj,  'b-', label='direct')
            l2, = ax.plot(dates, Bproj2, 'g-', label='reflect')
            if (i==0):
                ax.legend(handles=[l1,l2])
            ax.set_ylabel('delay%s (m)' % BLname)
            ax.set_title('delay between Ant%d-Ant%d' % (ai+1, aj+1))
            # variable freq phase plot
            ax = sub[2]
            #ax.plot(dates, phase1var, label='direct')
            cc = 'C%d'%(i%10)
            fMHz = freq[i]
            ax.plot(dates, phase, color=cc, label='%.0fMHz'%fMHz)
            ax.plot(dates, phase2, color=cc, linestyle=':')
            ax.scatter(dates, obsPha, color=cc)
            ax.legend()
            ax.set_ylabel('vis.phase (rad)')
            ax.set_title('solid: direct phase. dotted: direct + %.2f*reflect' % ratt)
            ax.set_xlabel('Time (UTC)')


        fig.autofmt_xdate()
        fig.tight_layout(rect=[0,0.03,1,0.95])
        fig.suptitle('Sun trajectory, phase for Baseline %s, %s' % (BLname, UTstr))
        fig.savefig(png)
        plt.close(fig)


