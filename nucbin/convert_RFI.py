#!/usr/bin/env python

from loadh5 import *
import matplotlib.pyplot as plt
from astropy.time import Time
from datetime import timedelta



inp = sys.argv[0:]
pg  = inp.pop(0)

ntune   = 15
f0      = 220.    # MHz
nchan   = 1024
nInp    = 2
tsep    = 900  # wf time resolution in seconds

usage = '''
convert .spec.h5 files (of a given device) in SED_specs/
to a regularized format for subsequent plotting

syntax:
    %s <dev> [options]

    <dev> is 'dev1', 'dev2', ... etc.

    options are:
    --tsep TSEP             waterfall plot time resolution: seconds (%d)
    -m nTune                number of freq tunings (match loop_run.py setting)
    -f fcen0                starting central freq in MHz (match loop_run.py setting)

''' % (pg, tsep)

if (len(inp) < 1):
    sys.exit(usage)

while (inp):
    k = inp.pop(0)
    if (k == '--tsep'):
        tsep = float(inp.pop(0))
    elif (k == '-m'):
        ntune = int(inp.pop(0))
    elif (k == '-f'):
        f0 = float(inp.pop(0))
    elif (k.startswith('-')):
        sys.exit('unknown option: %s'%k)
    else:
        dev = k


files = glob('SED_specs/SED_data*_%s.spec.h5' % dev)
files.sort()
nfiles = len(files)

# file open times:
ftimes = []
for f in files:
    ii = f.find('_dev')
    tstr = f[ii-11:ii]
    #print(f, '-->', tstr)
    dt = datetime.strptime(tstr, '%y%m%d_%H%M')
    ftimes.append(dt)
ftimes = np.array(ftimes)
fmt = '%y%m%dH%H'
s1 = ftimes.min().strftime(fmt)
s2 = (ftimes.max()+timedelta(hours=1)).strftime(fmt)
t1 = Time(datetime.strptime(s1, fmt), format='datetime').to_value('unix')
t2 = Time(datetime.strptime(s2, fmt), format='datetime').to_value('unix')
print('bracket times:', s1, s2, t1, t2)
totsec = t2 - t1
nsamp = int(totsec/tsep)
tsamp = []
for i in range(nsamp):
    tdelt = tsep * (0.5+i)
    tt = t1 + tdelt
    tsamp.append(tt)
tsamp = np.array(tsamp)


#fout = 'regularized_%s.spec.h5' % dev
fout = 'grid_%s-%s_%s.spec.h5' % (dev, s1, s2)

fcen = np.zeros(ntune)
for i in range(ntune):
    fcen[i] = f0 + 40.*i


allspec = np.ma.zeros((nsamp, ntune, nInp, nchan))
allspec.mask = np.ones_like(allspec, dtype=bool)    # all masked by default
allfreq = np.ma.zeros((nsamp, ntune, nchan))
allfreq.mask = np.ones_like(allfreq, dtype=bool)    # all masked by default
alltime = np.ma.zeros((nsamp, ntune))
alltime.mask = np.ones_like(alltime, dtype=bool)    # all masked by default

print('nsamp:', nsamp)
for k in range(nfiles):
    f = files[k]
    print('reading ', f, '  ...', '(%d/%d)'%(k,nfiles))

    ispec = getData(f, 'specdB')
    ifreq = getData(f, 'freqMHz')
    itime = getData(f, 'unixtime')
    medtime = np.median(itime)
    i = np.argmin(np.abs(tsamp-medtime))
    print(medtime, '-->', i, tsamp[i], medtime-tsamp[i])
    if (ispec.shape == (ntune, nInp, nchan)):
        allspec[i] = ispec
        allspec.mask[i] = False
        allfreq[i] = ifreq
        allfreq.mask[i] = False
        alltime[i] = itime
        alltime.mask[i] = False
    else:
        nfq = len(ifreq)
        ifcen = ifreq.mean(axis=1)
        for j in range(nfq):
            jdx = np.argmin(np.abs(fcen-ifcen[j]))
            allspec[i,jdx] = ispec[j]
            allspec.mask[i,jdx] = False
            allfreq[i,jdx] = ifreq[j]
            allfreq.mask[i,jdx] = False
            alltime[i,jdx] = itime[j]
            alltime.mask[i,jdx] = False

print('saving file:', fout, '...')
adoneh5(fout, allspec, 'specdB')
adoneh5(fout, allfreq, 'freqMHz')
adoneh5(fout, alltime, 'unixtime')
adoneh5(fout, tsamp, 'gridtime')    # also in unix time

