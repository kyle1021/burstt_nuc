#!/usr/bin/env python

from read_cli2 import *
from loadh5 import *

inp = sys.argv[0:]
pg  = inp.pop(0)

seg_byte = 131072   # bytes per segment
acount0  = 32768    # fiducial num of samples per segment
use_meta = True
redo     = False    # skip existing header file
verbose  = 0        # turn off info display

usage = '''
generate a header file for the data saved with bladeRF-cli or bladeRF-meta

usage:
%s <files> [options]

options are:
--no-meta       data is binary without meta
                (default is with meta)

--redo          re-generate header file even if it exists

--verbose       turn on some info display
-v              (same as --verbose)

''' % (pg,)

if (len(inp) < 1):
    sys.exit(usage)


files = []
while (inp):
    k = inp.pop(0)
    if (k == '--no-meta'):
        use_meta = False
    elif (k == '--redo'):
        redo = True
    elif (k == '--verbose' or k == '-v'):
        verbose = 1
    elif (k.startswith('-')):
        sys.exit('unknown option: %s' % k)
    else:
        files.append(k)


while (files):
    f = files.pop(0)
    if (f.endswith('.log')):
        name = f.replace('.log','')
    elif (f.endswith('.meta')):
        name = f.replace('.meta','')
    elif (f.endswith('.bin')):
        name = f.replace('.bin','')
    elif (f.endswith('.h5')):
        name = f.replace('.h5','')

    fhdr = '%s.h5' % name
    flog = '%s.log' % name
    if (use_meta):
        fbin = '%s.meta' % name
    else:
        fbin = '%s.bin' % name

    if (flog in files):
        files.remove(flog)
    if (fbin in files):
        files.remove(fbin)


    if (os.path.isfile(fhdr)):
        if (not redo):
            print('skip existing header: %s' % fhdr)
            continue

    fs = os.path.getsize(fbin)
    if (fs == 0):
        print('skip 0 size: %s' % fbin)
        continue


    print('processing: %s' % name)
    attrs = logParse(flog)
    nInp = attrs['nInp']
    if (use_meta):
        tick, acount, overrun = extractMeta(fbin, verbose=verbose)
    else:
        nseg = fs//seg_byte     # num of segments in the file
        nsamp = acount0 // nInp # num of clocks per segment
        tick = np.arange(nseg, dtype=np.uint64)*nsamp
        acount = np.ones(nseg, dtype=np.uint16)*acount0
        overrun = np.zeros(nseg, dtype=uint8)

    adoneh5(fhdr, tick, 'fpga_tick')
    adoneh5(fhdr, acount, 'actual_count')
    adoneh5(fhdr, overrun, 'overrun')
    putAttrs(fhdr, attrs)


