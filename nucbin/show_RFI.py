#!/usr/bin/env python

from loadh5 import *
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.time import Time

home = os.getenv('HOME')

inp = sys.argv[0:]
pg  = inp.pop(0)

nInp  = 2
#ntune = 15
edge  = 4
#nchan = 1024
ylim  = [-130, -55]
tfmt  = '%y%m%d %H'
tfmt2 = '%y%m%dH%H'
t1str = ''
t2str = ''
vmin  = 999.
vmax  = 999.


usage = '''
given a regularized spec file, produce RFI and WF plots

syntax:
    %s <grid_file> [options]

    options are:
    --t1 'yymmdd HH'     starting plot time (default: start of grid time)
    --t2 'yymmdd HH'     ending plot time (default: end of grid time)
    --vlim vmin vmax     set color scale limits for waterfall plots


''' % (pg,)

if (len(inp) < 1):
    sys.exit(usage)

while (inp):
    k = inp.pop(0)
    if (k == '--t1'):
        tmp = inp.pop(0)
        t1 = datetime.strptime(tmp, tfmt)
        t1str = t1.strftime(tfmt2)
    elif (k == '--t2'):
        tmp = inp.pop(0)
        t2 = datetime.strptime(tmp, tfmt)
        t2str = t2.strftime(tfmt2)
    elif (k == '--vlim'):
        vmin = float(inp.pop(0))
        vmax = float(inp.pop(0))
    elif (k.startswith('-')):
        sys.exit('unknown option: %s'%k)
    else:
        fspec = k

#fspec = 'regularized_%s.spec.h5' % dev
ii = fspec.find('_dev')
tmp = fspec[ii+1:]
ii2 = tmp.find('-')
dev = tmp[:ii2]
print('dev:', dev)

allspec = getData(fspec, 'specdB')
## plot_sed2.py saved corrected power in dBm now, no need to adjust
#allspec[:,8] -= 10.   # 540MHz
#allspec[:,9] -= 10.   # 580MHz
#allspec[:,14] -= -20. # 780MHz
allfreq = getData(fspec, 'freqMHz')
allepoch = getData(fspec, 'unixtime')
alldatetime = Time(allepoch, format='unix').to_value('datetime')
alldatetime = np.ma.array(alldatetime, mask=allepoch.mask)
#day1 = alldatetime[0,0].strftime('%y%m%d')
#day2 = alldatetime[-1,-1].strftime('%y%m%d')
gridepoch = getData(fspec, 'gridtime')
griddatetime = Time(gridepoch, format='unix').to_value('datetime')
if (t1str == ''):
    t1 = griddatetime[0]
    t1str = t1.strftime(tfmt2)
if (t2str == ''):
    t2 = griddatetime[-1]
    t2str = t2.strftime(tfmt2)

nwin, ntune, nchan = allfreq.shape


#fref = '%s/cli-command/220422_cafe/hotload_%s.spec.h5' % (home, dev)
#reffreq = getData(fref, 'freqMHz')
#refspec = getData(fref, 'specdB')
fref = '%s/local/bin/blade_terminate.h5' % (home,)
reffreq = getData(fref, 'freqMHz')
#refspec = getData(fref, 'specdB')
refspec1 = getData(fref, 'LNAon_40dB')
#refspec1 -= 20. # scale to 20dB
refspec2 = getData(fref, 'LNAoff_40dB')
#refspec2 -= 20. # scale to 20dB

## fig 1: median spec and min/max
wtime = np.logical_and(griddatetime>=t1, griddatetime<=t2)
medspec = np.ma.median(allspec[wtime], axis=0)
minspec = allspec[wtime].min(axis=0)
maxspec = allspec[wtime].max(axis=0)
#avgfreq = np.ma.mean(allfreq, axis=0)
avgfreq = np.ma.median(allfreq, axis=0)

filledspec = np.ma.filled(allspec, np.nan)
p15spec = np.nanpercentile(filledspec[wtime], 15, axis=0)
p85spec = np.nanpercentile(filledspec[wtime], 85, axis=0)

f1, s1 = plt.subplots(2,1,figsize=(15,9), sharex=True, sharey=True)
chw   = range(edge, nchan-edge)
for ch in range(nInp):
    ax = s1[ch]
    lh = [] # handles for legend
    ll = [] # labels for legend
    for i in range(len(reffreq)):
        fl1, = ax.plot(reffreq[i,chw], refspec1[i,chw], color='gray', linestyle='--', label='TermLNAon')
        fl2, = ax.plot(reffreq[i,chw], refspec2[i,chw], color='gray', linestyle=':', label='TermLNAoff')
        if (i == 0):
            lh.extend([fl1, fl2])
            ll.extend(['TermLNAon', 'TermLNAoff'])

    for i in range(ntune):
        cc = 'C%d'%(i%10)
        #ax.plot(avgfreq[i,chw], medspec[i,ch,chw], color=cc)
        fl1 = ax.fill_between(avgfreq[i,chw], maxspec[i,ch,chw], minspec[i,ch,chw], color=cc, alpha=0.2, label='min/max')
        fl2 = ax.fill_between(avgfreq[i,chw], p85spec[i,ch,chw], p15spec[i,ch,chw], color=cc, alpha=0.6, label='15%/85% pct')
        if (i == 0):
            lh.extend([fl1, fl2])
            ll.extend(['min/max', '15%/85% pct'])

    ax.set_xlim(avgfreq.min(), avgfreq.max())
    ax.set_xlabel('freq (MHz)')
    ax.set_ylabel('ampld (dB)')
    ax.legend(lh, ll)
    ax.set_title('ch%d'%ch)
    ax.set_ylim(ylim)

f1.tight_layout(rect=[0,0.03,1,0.95])
f1.suptitle('%s--%s RFI, %s' % (t1str, t2str, dev))
f1.savefig('RFI_%s-%s_%s.png' % (t1str, t2str, dev))
plt.close(f1)


## fig 2: waterfall of relative changes
f2, s2 = plt.subplots(2,1,figsize=(15,9))
relspec = allspec - medspec.reshape((1,ntune,nInp,nchan))
#onedatetime = alldatetime[:,0]
#rowepoch = allepoch.mean(axis=1)
#onedatetime = Time(rowepoch, format='unix').to_value('datetime')
for ch in range(nInp):
    ax = s2[ch]
    if (vmin == 999.):
        vmin = relspec[:,:,ch,chw].min()/2
    if (vmax == 999.):
        vmax = relspec[:,:,ch,chw].max()/2
    for j in range(ntune):
        #m=ax.pcolormesh(avgfreq[j,chw], griddatetime, relspec[:,j,ch,chw], vmin=vmin, vmax=vmax)
        #m=ax.pcolormesh(avgfreq[j,chw], onedatetime, relspec[:,j,ch,chw], vmin=vmin, vmax=vmax)
        #print('j=',j,'min/max=',relspec[:,j,ch,chw].min(), relspec[:,j,ch,chw].max())
        X, Y = np.meshgrid(griddatetime, avgfreq[j,chw], indexing='xy')
        #X, Y = np.meshgrid(gridepoch, avgfreq[j,chw], indexing='xy')
        m=ax.pcolormesh(X, Y, relspec[:,j,ch,chw].T, vmin=vmin, vmax=vmax, shading='nearest')
    ax.set_ylabel('freq (MHz)')
    #ax.set_xlabel('Unix Epoch (sec)')
    ax.set_ylabel('Datetime')
    #ax.set_xlabel('freq (MHz)')
    #ax.set_ylabel('Datetime')
    ax.set_title('ch%d'%ch)
    cb = plt.colorbar(m, ax=ax)
    cb.set_label('power/med (dB)')
    ax.set_xlim([t1, t2])

f2.tight_layout(rect=[0,0.03,1,0.95])
f2.suptitle('%s--%s RFI, %s' % (t1str, t2str, dev))
f2.savefig('WF_%s-%s_%s.png' % (t1str, t2str, dev))
plt.close(f2)
