#!/usr/bin/env python

from pyplanet import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import timedelta

site = 'taroge4'
b, obs = obsBody('sun', time=datetime(2021, 12, 22, 22, 20, 0), site=site, retOBS=True)
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
print(antXYZ[0], antXYZ[2], BVec)

Bproj = np.dot(BVec, unitVec)
print(Bproj)

t0 = datetime(2021,12,22,22,20,0)
dates = []
az = []
el = []
pbase = []
lamb = []
freq = []
for dmin in np.arange(0,210,1,dtype=float):
    d10 = timedelta(minutes=dmin)
    ti = t0 + d10
    dates.append(ti)
    fMHz = 260. + 5.*(dmin//10.)   # in MHz
    lam = 3.e8 / (fMHz * 1.e6)     # in meters
    #print(dmin, fMHz, lam)
    lamb.append(lam)
    freq.append(fMHz)

    obs.date = ti
    b.compute(obs)
    az.append(b.az)
    el.append(b.alt)

dates = np.array(dates)
az = np.array(az)
el = np.array(el)
lamb = np.array(lamb)
freq = np.array(freq)
#print(lamb)

lamb260 = 3.e8 / 260.e6 # for 260MHz, in meters

theta = np.pi/2. - el
phi   = np.pi/2. - az
theta2 = np.pi/2. + el  # reflection off the see
unitVec  = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]).T
unitVec2 = np.array([np.sin(theta2)*np.cos(phi), np.sin(theta2)*np.sin(phi), np.cos(theta2)]).T
Bproj = np.zeros_like(theta)
Bproj2 = np.zeros_like(theta2)
phase1 = np.zeros_like(Bproj)
phase2 = np.zeros_like(Bproj)
phase1var = np.zeros_like(Bproj)
for i in range(len(unitVec)):
    Bproj[i]  = np.dot(BVec, unitVec[i])  # in meters
    Bproj2[i] = np.dot(BVec, unitVec2[i]) # in meters
    phase1[i] = np.angle(np.exp(2.j*np.pi*Bproj[i]/lamb260))
    phase2[i] = np.angle(np.exp(2.j*np.pi*Bproj[i]/lamb260) + 0.5*np.exp(2.j*np.pi*lamb260/Bproj2[i]))    # attenuate reflection with an arbitrary factor (0.5)
    phase1var[i] = np.angle(np.exp(2.j*np.pi*Bproj[i]/lamb[i]))


fig, sub = plt.subplots(4,1,figsize=(10,12), sharex=True)
# az,el plot
ax = sub[0]
ax.plot(dates, az/np.pi*180., label='az')
ax.plot(dates, el/np.pi*180., label='el')
ax.set_ylabel('AZ,EL (deg)')
ax.set_title('Sun trajectory')
ax.legend()
# delay plot
ax = sub[1]
ax.plot(dates, Bproj)
ax.set_ylabel('delay1-3 (m)')
ax.set_title('delay between ant1-ant3')
# 260MHz phase plot
ax = sub[2]
ax.plot(dates, phase1, label='direct')
ax.plot(dates, phase2, label='direct+0.5reflect')
ax.set_ylabel('vis.phase (rad)')
ax.set_title('phase at 260MHz')
ax.legend()
# variable freq phase plot
ax = sub[3]
#ax.plot(dates, phase1var, label='direct')
ufreq = np.unique(freq)
for i in range(len(ufreq)):
    cc = 'C%d'%(i%10)
    fMHz = ufreq[i]
    w = freq == fMHz
    print(i, fMHz, freq[w])
    ax.plot(dates[w], phase1var[w], color=cc, label='%.0fMHz'%fMHz)
ax.legend()
ax.set_ylabel('vis.phase (rad)')
ax.set_title('phase from 260 to 360 MHz')
ax.set_xlabel('Time (UTC)')

fig.autofmt_xdate()
fig.suptitle('Sun trajectory, phase for Baseline 1-3, 2021/12/22')
fig.savefig('sun_trajectory.png')
plt.close(fig)

