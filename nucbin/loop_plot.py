#!/usr/bin/env python

import sys, os, os.path
from subprocess import call
from datetime import datetime, timedelta
import numpy as np
import time
#from mp_func import *

home = os.getenv('HOME')

def getParam(fconf, keyword):
    '''
    get a list of params with the keyword
    defined in the .conf file
    example:
        keyword = '@server'
    '''
    fh = open(fconf, 'r')
    tmp = fh.readlines()
    while (tmp):
        k = tmp.pop(0)
        if (k.startswith('#')):
            continue

        if (keyword in k):
            sub = k[k.find('(')+1:k.find(')')]
            param = sub.split(',')
    param = [x.strip().replace("'",'') for x in param]
    return param
        

## with single device, dual channel, the max bandwidth is 40MHz
#fconf = 'uni-dev4.conf'
#bw    = 40      # bandwidth in MHz

## if both devices are used (4 channels in total) on a single PC, the max bandwith is 20MHz
## however, for two devices on two PCs, the max bandwidth is still 40MHz
fconf = 'user.conf'
bw    = 40

fcen0 = 220     # starting central freq in MHz
ntune = 15      # 15 tunings, 220--780MHz

gain0 = 40.     # starting bladeRF gain in dB
glist = None
gentry = {}
name0 = 'data'
do_sed = True
sopt  = ''
do_header = False

# for testing
#dur   = 1       # recording duration
#odir0 = '.'
#odir0 = '/pool-frb/cli-command/test2'
# for 1h data taking
#dur   = 235
# for 15min data taking
dur   = 20      # reduce data size: 20sec*15tunings = 300sec
odir0 = '.'
loop_dur = 900  # wait to begin next loop 15min
bw20  = False



inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
measure a wide bandwidth with multiplt tunings
can loop the tunings multiple times
optionally chaning the gain in each loop

syntax:
    %s N [options]

    N is the number of loops to repeat the whole sweep

options are:
-m nTune        # of tunings in each loop (%d)
-d duration     how long to record at each tuning (%.0f sec)
--ld loop_duration
                time btween loop in seconds (%d sec)

-f fcen0        starting central freq in MHz (%.0f MHz)
-r rate         sample rate in Msps, equal to bandwidth in MHz
                this is also the separation between tunings
                (%.0f Msps)

-g gain         gain for the bladeRF (%.0f dB)
                i.e. glist = np.ones(ntune)*gain
--glist "g0 g1 g2 ... gN"
                explicitly set the gain of each freq
                number of elements (N) must be the same as nTune
                the list needs to be quoted
--gentry i G    modify the i-th entry of the glist to G (dB)
                can be combined with -g or --glist
                ** note: i starts from 0

-o odir         the parent directory where data will be saved
                actual data is saved under %s_yymmdd_HHMM/
                (%s/)
--name name     string used in output folder and filename
                (%s)

--test          set name=test and dur=1 (sec)

--header        use 'cli2header.py' after each loop
''' % (pg, ntune, dur, loop_dur, fcen0, bw, gain0, name0, odir0, name0)

if (len(inp) < 1):
    sys.exit(usage)

while (inp):
    k = inp.pop(0)
    if (k == '-m'):
        ntune = int(inp.pop(0))
    elif (k == '-d'):
        dur = int(inp.pop(0))
    elif (k == '-r'):
        rate = int(inp.pop(0))
        bw = rate
        if (np.isclose(rate, 20.)):
            bw20 = True
            print('set bw20 True')
    elif (k == '-f'):
        fcen0 = int(inp.pop(0))
    elif (k == '-g'):
        gain0 = int(inp.pop(0))
    elif (k == '-o'):
        odir0 = inp.pop(0)
    elif (k == '--name'):
        name0 = inp.pop(0)
    elif (k == '--test'):
        name0 = 'test'
        dur = 1
        #loop_dur = 300
        loop_dur = 30
        sopt = '--name test'
    elif (k == '--glist'):
        tmp = inp.pop(0)
        glist = [float(x) for x in tmp.split()]
    elif (k == '--gentry'):
        istr = inp.pop(0)   # i as a str
        gentry[istr] = float(inp.pop(0))
    elif (k == '--ld'):
        loop_dur = int(inp.pop(0))
    elif (k == '--header'):
        do_header = True
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        nloop = int(k)


fqlist = np.arange(ntune)*bw + fcen0

if (glist is None):
    glist = np.ones(ntune) * gain0
else:
    if (len(glist) != ntune):
        sys.exit('error setting glist! num of element inconsistent with nTune.')
for k in gentry.keys():
    ientry = int(k)
    glist[ientry] = gentry[k]

if (False): # debugging
    for i in range(ntune):
        print(i, fqlist[i], glist[i])
    sys.exit()

server = getParam(fconf, '@server')
device = getParam(fconf, '@active')
cwd = os.getcwd()
base_cwd = os.path.basename(cwd)
dir_cwd  = os.path.dirname(cwd)
parent = getParam(fconf, '@parent')
for i in range(len(parent)):
    if parent[i] == '$ldir':
        parent[i] = dir_cwd
#print(parent)

for j in range(nloop):
    #gain = gain0 - j*10     # 20, 10, 0, -10
    #gain = gain0 * 1
    now = datetime.now()
    tstr = now.strftime('%y%m%d_%H%M')
    name = '%s_%s/%s' % (name0, tstr, name0)
    if (odir0 != '.'):
        name = '%s/%s' % (odir0, name)
        odir = '%s/%s_%s' % (odir0, name0, tstr)
    else:
        odir = '%s_%s' % (name0, tstr)

    if (not os.path.isdir(odir)):
        call('mkdir -p %s'%odir, shell=True)
    #os.chdir(odir0)

    for i in range(ntune):
        fcen = fcen0 + i*bw
        print('fcen=%.0fMHz'%fcen)
        #if (np.isclose(fcen, 540.) or np.isclose(fcen, 580.)):
        #    gain = gain0 - 10.
        #elif (np.isclose(fcen, 780.)):
        #    gain = gain0 - 20.
        #else:
        #    gain = gain0 - 0.
        #gain = gain0 - 0.
        gain = glist[i]
        cmd = 'run_blade.pl %.0f -r %.0f -f %.0f -g %.1f --name %s --conf %s'%(dur, bw, fcen, gain, name, fconf)
        call(cmd, shell=True)
        #print(cmd)

    # produce SED for each loop
    if (do_sed):
        pdir = 'SED_plots'
        call('mkdir -p %s'%pdir, shell=True)
        sdir = 'SED_specs'
        call('mkdir -p %s'%sdir, shell=True)

        cmd0 = 'plot_sed2.py %s %s'%(odir, sopt)
        children = []
        for si, srv in enumerate(server):
            pid = os.fork()
            if (pid > 0):   # the parend
                children.append(pid)
                print(srv, pid)
            else:           # the child
                dev = device[si]
                cwd = '%s/%s' % (parent[si], base_cwd)
                #fref = '%s/cli-command/220422_cafe/hotload_%s.spec.h5' % (home, dev)
                fref = '%s/local/bin/blade_terminate.h5' % (home,)
                cmd = cmd0 + ' --ref %s --ylim -130 -55' % (fref,) 
                cmd += ' --dev %s'%dev
                if (bw20):
                    cmd += ' --bw20'
                #print('common cmd:')
                #print(cmd)

                if (srv == 'local'):
                    print('plotting on: local...')
                    #print(cmd)
                    #os._exit(0)
                    call(cmd, shell=True)
                    call('cp %s/SED*.png %s'%(odir, pdir), shell=True)
                    call('cp %s/SED*.spec.h5 %s'%(odir, sdir), shell=True)
                    if (do_header):
                        cmd = 'cli2header.py %s/*.meta' % odir
                        call(cmd, shell=True)
                else:
                    cmd1 = 'ssh %s -T "cd %s; ~/miniconda3/bin/python ~/local/bin/%s"' % (srv, cwd, cmd)
                    print('plotting on: %s ...' % srv)
                    #print(cmd1)
                    #os._exit(0)
                    call(cmd1, shell=True)
                    cmd2 = 'scp %s:%s/%s/SED*.png %s'%(srv, cwd, odir, pdir)
                    call(cmd2, shell=True)
                    cmd2 = 'scp %s:%s/%s/SED*.spec.h5 %s'%(srv, cwd, odir, sdir)
                    call(cmd2, shell=True)
                    if (do_header):
                        cmd = 'ssh %s -T "cd %s; ~/miniconda3/bin/python ~/local/bin/cli2header.py %s/*.meta"' % (srv, cwd, odir)
                        call(cmd, shell=True)
                os._exit(0) # important to finish the child here

        # wait all children to finish
        for pid in children:
            print('waiting', pid)
            os.waitpid(pid, 0)
        print('parallel plotting done.')

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

#cmd2 = 'quick_look.py --disp %s*.meta' % (name,)
#call(cmd2, shell=True)

