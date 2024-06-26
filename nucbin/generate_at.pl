#!/usr/bin/env python

import sys, os.path
from subprocess import call
from datetime import datetime, timedelta

# the bladeRF command to use
blade = 'run_blade_meta-nuc-a4.pl'

# default starting time
now = datetime.now()    # in local time
start = now + timedelta(minutes=1)
tstr = start.strftime('%y%m%d %H%M')

# other defaults
name  = 'test'  # at-test.sh, test_yymmdd_HHMMSS_dev?.meta/.log
rate  = 40      # Msps
fcen0 = 490     # MHz
nTune = 10      # num of tunings/at commands
delMin = 6      # separation of the at commands in minutes
opt   = ''



inp = sys.argv[0:]
pg  = inp.pop(0)

usage = '''
generate shell scripts that executes the bladeRF command using "at" system command

usage:
    %s <NAME> [options]

    <NAME> is used for script name: "at-NAME.sh"
            and also data file names: "NAME_yymmdd_HHMMSS_dev?.meta/.log"

options are:
    -r RATE     sample rate in Msps (%.1f)
    -f FCEN     starting central freq in MHz (%.1f)
    -n N        N commands (%d)
                each command is separated in central freq by the sample rate (bandwidth)
    -d DMIN     separation of the commands in minutes (%.1f)
    -t "yymmdd HHMM"    starting time of the commands in local time
                        (%s)
    -o OPT      common option to be passed to the run_blade command

''' % (pg, rate, fcen0, nTune, delMin, tstr)

if (len(inp)<1):
    sys.exit(usage)


while (inp):
    k = inp.pop(0)
    if (k == '-r'):
        rate = float(inp.pop(0))
    elif (k == '-f'):
        fcen0 = float(inp.pop(0))
    elif (k == '-n'):
        nTune = int(inp.pop(0))
    elif (k == '-d'):
        delMin = float(inp.pop(0))
    elif (k == '-t'):
        tstr = inp.pop(0)
    elif (k == '-o'):
        opt = inp.pop(0)
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        name = k

durSec = delMin*60 - 5  # recording duration for each tuning in seconds

ofile = 'at-%s.sh' % name
OUT = open(ofile, 'w')

t0 = datetime.strptime(tstr, '%y%m%d %H%M')
for i in range(nTune):
    fcen = fcen0 + i*rate
    t0 += timedelta(minutes=delMin)
    #atstr = t0.strftime('%H:%M %y%m%d')
    atstr = t0.strftime('%H:%M %m%d%y')

    bcmd = "%s %.0f -r %.1f -f %.1f --name %s " % (blade, durSec, rate, fcen, name)
    bcmd += opt
    cmd = "echo '%s'" % bcmd
    cmd += ' | at %s' % atstr

    print(cmd, file=OUT)

OUT.close()
print('script saved:', ofile)

