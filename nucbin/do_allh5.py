#!/usr/bin/env python

from analysisconf import *

dirs = glob('data_*')
dirs.sort()
ndir = len(dirs)

cmd = 'cli2header.py *.meta'

for di, d in enumerate(dirs):
    print(d, '(%d/%d)'%(di+1,ndir))
    os.chdir(d)
    call(cmd, shell=True)
    os.chdir('..')

    #if (di > 1):
    #    break


