#!/usr/bin/env python

import sys, os.path
from subprocess import call, run, PIPE
from datetime import datetime, timedelta
import time

home = os.getenv('HOME')

## with single device, dual channel, the max bandwidth is 40MHz
#fconf = 'uni-dev4.conf'
#bw    = 40      # bandwidth in MHz

## if both devices are used (4 channels in total), the max bandwith is 20MHz
fconf = 'user.conf'
bw    = 40

fcen0 = 220     # starting central freq in MHz
ntune = 15      # 15 tunings, 220--780MHz

gain0 = 40.     # starting bladeRF gain in dB
dur   = 10      # recording duration
loop_dur = 900  # time between loops
odir0 = '/pool-frb/cli-command/2022-dongsha'
#dur   = 1       # recording duration
#odir0 = '/pool-frb/cli-command/test'
#odir0 = '.'
name0 = 'data'
show  = False
do_check = False
do_sed = True
dev = 'dev3'    # dongsha setting
pdir = '%s/data_SED' % home
do_test = False
nLoop = 0


inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
run N hours of recording

each loop contains 15 frequency tunings with %dsec recording at each tuning
(next loop will start after %dsec)

N hours will be translated to nLoop

syntax:
    %s N [options]

    for --test, N is nLoop
    otherwise, N is nHour (use --loop N to specify nLoop outside of --test)

    --loop N        specify nLoop directly
    --test          run a short test
                    each tuning is 1sec, next loop to start after 60sec instead
    --check         produce quick_look2 plots (histograms...)
    --no-disp       skip showing the png files

''' % (dur, loop_dur, pg)

if (len(inp) < 1):
    sys.exit(usage)

while (inp):
    k = inp.pop(0)
    if (k == '--no-disp'):
        show = False
    elif (k == '--loop'):
        nLoop = int(inp.pop(0))
    elif (k == '--test'):
        dur = 1
        loop_dur = 60   # next loop to start immediately if plotting took too long
        odir0 = '/pool-frb/cli-command/test'
        pdir = '%s/test_SED' % home
        show = True
        do_test = True
    elif (k == '--check'):
        do_check = True
    elif (k.startswith('-')):
        sys.exit('unknown option: %s'%k)
    else:
        nVar = int(k)

if (do_test):
    nLoop = nVar
    print('test: %d loops' % nLoop)
else:
    if (nLoop > 0):
        print('data: %d loops' % nLoop)
    else:
        nHour = nVar
        nLoop = int(nHour*3600/loop_dur)
        print('data: %d hours --> %d loops' % (nHour, nLoop))

if (not os.path.isfile(fconf)):
    sys.exit('config file not found: %s'%fconf)

cmd = 'lsusb|grep "Nuand bladeRF 2.0"'
#cmd = 'bladeRF-meta -e exit'
result = run(cmd, shell=True, stdout=PIPE)
cmdout = result.stdout.decode().strip()
if (cmdout == ''):
    sys.exit('bladeRF not connected?')



call('mkdir -p %s'%odir0, shell=True)
if (not os.path.isdir(odir0)):
    sys.exit('failed to create output dir: %s'%odir0)

if (not os.path.isdir(pdir)):
    call('mkdir -p %s'%pdir, shell=True)


for j in range(nLoop):
    gain = gain0
    now = datetime.now()    # local time
    tstr = now.strftime('%y%m%d_%H%M')
    odir = '%s/%s_%s' % (odir0, name0, tstr)
    call('mkdir -p %s'%odir, shell=True)

    name = '%s/%s' % (odir, name0)

    for i in range(ntune):
        fcen = fcen0 + i*bw
        print('loop:%d, tuning:%d, fcen%.0fMHz'%(j,i,fcen))
        cmd = 'run_blade.pl %.0f -r %.0f -f %.0f -g %.1f --name %s --conf %s'%(dur, bw, fcen, gain, name, fconf)
        call(cmd, shell=True)
        #print(cmd)

    if (do_check):
        cmd2 = 'quick_look2.py %s*.meta' % (name,)
        if (show):
            cmd2 += ' --disp &'
        call(cmd2, shell=True)

    if (do_sed):
        cmd2 = 'plot_sed2.py %s --dev %s' % (odir,dev)
        call(cmd2, shell=True)
        cmd2 = 'cp -f %s/SED* %s/' % (odir, pdir)
        call(cmd2, shell=True)
        if (show):
            call('display %s/SED*.png &'%odir, shell=True)

    newloop = now + timedelta(seconds=loop_dur)
    print('next loop to start at:', newloop)
    now2 = datetime.now()
    ii = 0
    while (now2 < newloop):
        time.sleep(1)   # wait 1 sec till specified time reached
        ii += 1
        now2 = datetime.now()
        if (ii%60==0):
            print('...now:', now2)

