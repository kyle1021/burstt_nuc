#!/usr/bin/env python


import os, sys, os.path
from subprocess import call


host = os.getenv('HOSTNAME')
pwd  = os.getcwd()
print(host, pwd)
if (not host=='coma18' or not pwd=='/coma/Dropbox/BRST_project/NUC_local/bin'):
    sys.exit('this script is intended to distribute from coma18 to all NUCs')

dest = [
#    'nuc4',
#    'dongsha',
    'cafe2',
    'cafe3',
#    'cafe5',
    'cafe6',
    'nantou7',
    'nantou8',
    'lyudao9', #'nuc9',
    'nuc11',
    'nuc12',
    'lngabul'
    ]


for server in dest:
    print('server:', server)
    cmd = 'rsync -avu . -e ssh %s:~/local/bin/' % server
    #print(cmd)
    call(cmd, shell=True)


