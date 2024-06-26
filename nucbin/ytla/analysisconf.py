#!/usr/bin/env python

import numpy as np
import pandas as pd
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import sys, os.path
from datetime import datetime
from subprocess import call
from glob import glob
import platform
#print(platform.python_version())

from warnings import simplefilter
simplefilter('ignore')


## correlator related
na      = 7
nb      = int(na * (na-1) / 2)
nsb     = 2
nch     = 1024

bw0     = 2240.         # default bandwidth if not in the attrs
param_ver0  = '2010v2'  # default deformation model
ant_conf0   = '7x1.4'   # default antenna configuration


## OT related
OTremote1       = 'tcs:~/Oimage'
OTremote2       = 'macmini:~/Oimage'
OTremote3       = 'macmini2:~/Oimage'

## local setup
host    = platform.node()
#print(host)
#if (host == 'monitorytla'):
if (host == 'monitorYTLA2'):
    bindir      = '/Reduction/analysis/py38bin'
    obindir     = '/Reduction/analysis/py38obin'
    OTlocal     = '/Reduction/Oimage'
elif (host == 'amiba_corr'):
    #bindir     = '/home/corr/kylin/bin'
    bindir      = '/home/corr/kylin/testbin'
elif (host.startswith('Linde-MacBook-Pro')):
    bindir      = '/Users/kylin/analysis/py38bin'
else:
    bindir      = '.'   #??    

#print(bindir)
