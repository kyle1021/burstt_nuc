#!/usr/bin/env python

import struct
import sys, os.path
import time


inp = sys.argv[0:]
pg  = inp.pop(0)

#-- defaults --
seg_byte    = 131072    # number of bytes in a segment
samp_byte   = 4         # number of bytes of a sample (I+Q)
seg_samp    = seg_byte // samp_byte # number of samples per segment


usage = '''
this script convert a dual-channel binary file recorded by bladeRF (no meta)
into two single-channel binary files

syntax:
    %s <file(s)>

    multiple files can be supplied. each will be processed sequentially.

''' % (pg,)

if (len(inp) < 1):
    sys.exit(usage)


files = []
while (inp):
    k = inp.pop(0)
    files.append(k)

print('files:', files)
nFile = len(files)

for i in range(nFile):
    fname = files[i]
    if (os.path.isfile(fname)):
        if (fname.endswith('.bin')):
            print('process:', fname)
        else:
            print('warning, file type may be wrong:', fname)
    else:
        print('.. skipped (no file):', fname)
        continue

    fsize = os.path.getsize(fname)

    fout1 = fname + '.rx1'
    fout2 = fname + '.rx2'

    t1 = time.time()
    fh  = open(fname, 'rb')
    fh1 = open(fout1, 'wb')
    fh2 = open(fout2, 'wb')
    while True:
        buf = fh.read(seg_byte)
        if (len(buf) > 0):
            nbyte = len(buf)
            nsamp = nbyte // samp_byte
            tmp = struct.unpack('<%di'%nsamp, buf)
            nhalf = nsamp // 2
            # even sample, Rx1
            tmp1 = tmp[0::2]
            rx1buf = struct.pack('<%di'%nhalf, *tmp1)
            fh1.write(rx1buf)
            # odd sample, Rx2
            tmp1 = tmp[1::2]
            rx2buf = struct.pack('<%di'%nhalf, *tmp1)
            fh2.write(rx2buf)
        else:   # reached EOF
            break

    fh.close()
    fh1.close()
    fh2.close()
    t2 = time.time()

    print('.. file size (bytes):', fsize)
    print('.. finished in %.2fsec' % (t2-t1))



