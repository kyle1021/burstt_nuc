#!/usr/bin/env python
from analysisconf import *

import h5py
#import astropy.stats.biweight as bwt
from astropy.stats import sigma_clip
from deformation.ytlaplatform import *
from pyplanet import *

from scipy import signal

#print('loadh5: test version')


def selectVis(mask, select, count=False):
    '''
    given a list of channel selections
    mask (True) the un-selected channels

    input mask is the visibility mask
    channel is axis=2

    select is a list of strings
    e.g. ['-510:-10', '377:379', 128]
    will select channels 10 to 510 in LSB
    channels 377 to 379 in USB
    and also channel 128 in USB
    '''
    if (count): # ignore input mask (which is required, so set to None)
        sh = (nsb,1,nch,1)
    else:
        sh = mask.shape

    in_mask = np.ones(sh, dtype=bool)   # init = all masked

    if (isinstance(select, (int,str))):
        if (isinstance(select, str)):
            sel = int(select)
        else:
            sel = select
        if (sel > 0):
            in_mask[1,:,sel,:] = False
        elif (sel < 0):
            in_mask[0,:,-sel,:] = False

    elif (isinstance(select, list)):
        for item in select:
            part = item.split(':')
            if (len(part) == 1):
                sel = int(part[0])
                if (sel > 0):
                    in_mask[1,:,sel,:] = False
                elif (sel < 0):
                    in_mask[0,:,-sel,:] = False
            elif (len(part) == 2):
                sel = [int(part[0]), int(part[1])]
                sel.sort()
                if (sel[0]>0 and sel[1]>0):
                    in_mask[1,:,sel[0]:sel[1]+1,:] = False
                elif (sel[0]<0 and sel[1]<0):
                    in_mask[0,:,-sel[1]:-sel[0]+1,:] = False
                else:   # sel[0]<=0 and sel[1]>=0
                    in_mask[1,:,:sel[1]+1,:] = False
                    in_mask[0,:,:-sel[0]+1,:] = False

    if (count):
        n_sel = np.count_nonzero(~in_mask)
        return n_sel
    else:
        return in_mask


def bldict():
    '''
    return a dictionary
    bl['ij'] = b
    '''
    bl = {}
    b = -1
    for ai in range(na-1):
        for aj in range(ai+1, na):
            b += 1
            corr = '%d%d' % (ai, aj)
            bl[corr] = b
    return bl


def ldcorr2(fbase, na=na):
    '''
    load the new RAW.h5 data
    Note: the phase of USB is flipped to be consistent with the LSB

    return: time, auto, cross
    '''
    dname = os.path.dirname(fbase)
    if (dname == ''):
        dname = '.'
    part = os.path.basename(fbase).split('.')
    fname = dname + '/%s.RAW.h5' % '.'.join(part[:2])

    try:
        f = h5py.File(fname, 'r')
        time = f.get('timestampEpoch')[()]
        ndata = len(time)
        #print ndata, time
    except:
        print('error accessing .RAW.h5 file:', fname)
        sys.exit()


    sbs = ['lsb', 'usb']
    auto  = np.zeros((nsb, na, nch, ndata))
    cross = np.zeros((nsb, nb, nch, ndata), dtype=complex)

    for s in range(nsb):
        # auto
        for ai in range(na):
            dname = '%s/auto%d%d' % (sbs[s], ai, ai)
            auto[s,ai] = f.get(dname)
    
        # cross
        b = -1
        for ai in range(na-1):
            for aj in range(ai+1, na):
                b += 1

                cross[s,b].real = f.get('%s/cross%d%d/real' % (sbs[s], ai, aj))[()]
                cross[s,b].imag = f.get('%s/cross%d%d/imag' % (sbs[s], ai, aj))[()]

    cross[1,:].imag *= -1.

    f.close()

    return time, auto, cross


def ldcorr(fname, na=na):
    '''
    load the raw data from correlator (dataout.sh)
    Note: the phase of USB is flipped to be consistent with the LSB

    return: time, auto, cross
    '''


    if (fname.endswith('.')):
        fname = fname.rstrip('.')
    elif (fname.endswith('.timestamp')):
        fname = fname.rstrip('.timestamp')

    ftime = fname + '.timestamp'

    if (os.path.isfile(ftime)):
        print('use existing .timestamp file')
    else:
        print('no existing .timestamp file')
        print('will use ad hoc time info (1 sec/pt)')

    dname = os.path.dirname(fname)
    if (dname == ''):
        dname = '.'
    fbase = fname


    print('... ', datetime.now().isoformat())


    #-- determine data length --
    aname = dname + '/' + fbase + '.%s.%s.h5' % ('lsb', 'auto')
    ha = h5py.File(aname, 'r')
    bname = 'auto%d%d' % (0, 0)
    arr = np.array(ha.get(bname))
    (nch, ndata) = arr.shape
    ha.close()

    #-- load data --
    if (os.path.isfile(ftime)):
        time = np.loadtxt(ftime)
    else:
        time = np.arange(ndata)
    

    w = ['lsb', 'usb']
    auto  = np.zeros((nsb, na, nch, ndata))
    cross = np.zeros((nsb, nb, nch, ndata), dtype=complex)
    for s in range(nsb):
        b = -1

        xtype = 'auto'
        aname = dname + '/' + fbase + '.%s.%s.h5' % (w[s], xtype)
        print(aname, '--> ', os.path.isfile(aname))
        ha = h5py.File(aname, 'r')
        xtype = 'cross'
        cname = dname + '/' + fbase + '.%s.%s.h5' % (w[s], xtype)
        print(cname, '--> ', os.path.isfile(cname))
        hc = h5py.File(cname, 'r')

        for i in range(na):
            bname = 'auto%d%d' % (i, i)
            arr = np.array(ha.get(bname))
            auto[s, i] = arr
                
                    
            for j in range(i+1, na):
                b += 1
                bname = '/cross%d%d' % (i, j)
                arr = np.array(hc.get('%s/real' % bname))
                cross[s, b].real = arr
                arr = np.array(hc.get('%s/imag' % bname))
                if (s == 1):    # flip the phase or imaginary part for USB
                    arr *= -1.
                cross[s, b].imag = arr

        ha.close()
        hc.close()

    print('... ', datetime.now().isoformat())
    return time, auto, cross


#def ldBandpass():


def ldoneh5(h5name):
    # load the dataset from All-in-one h5 file
    print('... ', datetime.now().isoformat())

    with h5py.File(h5name, 'r') as f:
        #na     = f.attrs['na']
        #nb     = f.attrs['nb']
        #nsb    = f.attrs['nsb']
        #nch    = f.attrs['nch']
        #ndata  = f.attrs['ndata']

        print('...  timestamp')
        time = np.array(f.get('timestamp'))
        print('...  auto-corr (real)')
        auto  = np.array(f.get('auto'))
        print('...  cross-corr (complex)')
        cross = np.array(f.get('cross'))

    print('... ', datetime.now().isoformat())
    return time, auto, cross


def modCalTable(fname, op='now', ctable=None):
    '''
    modify the cal_table of the H5 file
    op = one of the following: 'set', 'add', 'rem', 'now', 'all'
    ctable = a string or a list of tables
    
        set --> set the cal_table to ctable
        add --> add ctable to the existing cal_table
        rem --> remove ctable from the existing cal_table
        now --> show the existing cal_table, ignore ctable
        all --> show all available tables in the file, ignore ctable
    '''
    gname = 'calibration'   # default group name
    try:
        F = h5py.File(fname, 'a')
        pri_table = F.attrs['cal_table']
        g = F[gname]
        all_table = list(g.keys())
        if (isinstance(pri_table, (str,bytes))):
            pri_table = [pri_table]
        elif (isinstance(pri_table, np.ndarray)):
            pri_table = pri_table.tolist()
        # else it is assumed to be a list
    except:
        sys.exit('error accessing cal_table in: %s' % fname)

    if (op == 'now'):
        F.close()
        return pri_table

    if (op == 'all'):
        F.close()
        return all_table

    if (ctable is None):    # and op not: 'now', 'all'
        F.close()
        return None

    # regularize the tables as lists
    if (isinstance(ctable, (str,bytes))):
        ctable = [ctable]
    elif (isinstance(ctable, np.ndarray)):
        ctable = ctable.tolist()
    # else it is assumed to be a list


    # building the new list of tables
    if (op == 'add'):
        new_table = pri_table + ctable

    elif (op == 'rem'):
        new_table = pri_table[:]
        for t in ctable:
            if (t in new_table):
                new_table.remove(t)

    elif (op == 'set'):
        new_table = ctable[:]

    # checking duplicates and invalid tables
    fin_table = []
    all_table2 = all_table[:]
    all_table2.append('raw')    # make raw a valid table
    for t in new_table:
        if (t in all_table2 and not t in fin_table):
            fin_table.append(t)
    
    F.attrs['cal_table'] = fin_table
    F.close()

    return fin_table



def getAttrs(fname, dest='/'):
    '''
    retrieve the top-level attributes from the correlator data
    (representative file: .lsb.auto.h5)
    input:
        specific filename of the h5 file

        dest is a string for the destination group or dataset
        default is to get the root level attributes

    return: dict
    '''

    if (not os.path.isfile(fname)):
        print('error finding the source h5 file:', fname)
        sys.exit()

    f = h5py.File(fname, 'r')
    try:
        g = f.get(dest)
    except:
        #sys.exit('error accessing file: %s' % fname)
        print('destination not found: %s' % dest)
        return None

    attrs = {}
    for k in list(g.attrs.keys()):
        attrs[k] = g.attrs.get(k)
    f.close()

    return attrs


def putAttrs(fname, attrs, dest='/'):
    '''
    put the attrs (dict) as attributes in the oneh5 file
    dest is a string for the destination group or dataset
    default is to put in root level attributes
    '''
    try:
        f = h5py.File(fname, 'a')
        g = f.get(dest)

        for k in list(attrs.keys()):
            g.attrs[k] = attrs[k]
        f.close()
    except:
        print('error accessing oneh5 file:', fname)
        sys.exit()

    return 1


def getData(fname, dname):
    '''
    return the <dname> dataset array
    if <dname>.mask exist, the data is load as masked array
    '''
    #print('debug: in getData')
    try:
        f = h5py.File(fname, 'r')
    except:
        sys.exit('error opening file: %s' % fname)

    if (f.get(dname)):
        d = f[dname][()]
    else:
        #sys.exit('%s not found in %s' % (dname, fname))
        print('%s not found in %s' % (dname, fname))
        print('... return None')
        return None

    mname = '%s.mask' % dname
    if (f.get(mname)):
        m = f[mname][()]
        arr = np.ma.array(d, mask=m)
    else:
        arr = d

    f.close()

    return arr    


def chanAvg(avg_vis, avg_var, flag, nlen, chlim=[20,820], stdlim=5000.):
    '''
    channel average the visibility using inverse variance weighting
    input channels can be flagged if there are emission lines

    input:
        avg_vis[nsb, nb, nch, nunit], complex       time-averaged vis
        avg_var[nsb, nb, nch, nunit], complex       corresponding variance for real/imag parts
        flag[nsb, nb, nch, nunit], int              input flag: 0=good, >0=bad (bit-mask)
        nlen[nunit], int                            number of data points per unit
                                                    if nlen == None, the scatter will not be returned
                                                    also, the flag=8 will not be added
        chlim                                       channel range to average
        stdlim                                      a threshold for raw scatter, in Jy
                                                    5000 Jy is suitable for 5sec integration time

    output:
        int_vis[nsb, nb, nunit], complex            time- and channel- averaged vis
        int_var[nsb, nb, nunit], float              an estimate of variance in amplitude
        int_flag[nsb, nb, nunit], int               a new flag for baselines with excessive variance
        scatter[nsb, nb, nunit], float              estimate of raw data scatter (variance)
    '''
    avg_vis = avg_vis.reshape((nsb, nb, nch, -1))
    avg_var = avg_var.reshape((nsb, nb, nch, -1))
    flag    = flag.reshape((nsb, nb, nch, -1))
    sh = avg_vis.shape
    nunit = sh[3]

    if (isinstance(nlen, (list, np.ndarray, int))):
        do_scatter = True
        if (nunit == 1):
            nlen = np.array(nlen)
    else:
        # disable scatter
        do_scatter = False

    #chlim = np.array(chlim).reshape((-1,2))
    #if (chlim.shape == (1,2)):
    #   chlim = np.append(chlim, -chlim, axis=0)    # duplicate channel range for 

    ma_vis = np.ma.array(avg_vis, mask=flag)
    ma_var = np.ma.array(avg_var, mask=flag)

    ivw = 2. / (ma_var.real + ma_var.imag)      # inverse variance weighting
    #ivw = np.ma.array(2. / (ma_var.real + ma_var.imag), mask=flag)     # inverse variance weighting

    #sumwt = ivw[:,:,chlim[0]:chlim[1],:].sum(axis=2)
    #int_vis = np.sum((ma_vis*ivw)[:,:,chlim[0]:chlim[1],:], axis=2) / sumwt
    #int_var = 1. / sumwt
    sumwt   = np.ma.sum(ivw[:,:,chlim[0]:chlim[1],:], axis=2)
    int_vis = np.ma.sum((ma_vis*ivw)[:,:,chlim[0]:chlim[1],:], axis=2) / sumwt
    int_var = np.ma.array(1. / sumwt, mask=sumwt.mask)


    if (do_scatter):
        #int_flag = int_vis.mask.astype(int) * 8        # bit mask 8 for excessive variance after calibration
        #w = int_var > stdlim**2
        # raw scatter ~ int_var * nlen * (chlim[1]-chlim[0])
        scatter = int_var * (chlim[1]-chlim[0]+1)
        for ui in range(nunit):
            scatter[:,:,ui] *= nlen[ui]
        w = scatter > stdlim**2
        int_flag = np.zeros((nsb, nb, nunit), dtype=int)
        #print 'int_flag', int_flag.shape
        int_flag[w] += 8

        # scatter returns a flat variance for Jupiter on-source units (60sec, 13 points) and outlier units (120sec, 25 points)
        # scatter2 shows that outlier units has 3 times smaller variance than on-source units
        # scatter2 outlier variance is also about 2 times smaller thant scatter outlier units

        scatter2 = (np.ma.var(ma_vis[:,:,chlim[0]:chlim[1],:].real, axis=2) + np.ma.var(ma_vis[:,:,chlim[0]:chlim[1],:].imag, axis=2)) / 2.
        #scatter2 = np.ma.var(np.abs(ma_vis[:,:,chlim[0]:chlim[1],:]), axis=2)
        # since the purpose of scatter2 is to estimate the noise variance with outlier pointings (low to no signal),
        # we need to multiply the abs.variance by 1/0.42 (if we assume 0 SNR) to estimate the underlying noise variance
        # (assuming Gaussian random noise)
        # on the other hand, real.var and imag.var do not require this scaling
        #scatter2 /= 0.42
        for ui in range(nunit):
            scatter2[:,:,ui] *= nlen[ui]

        debug = False
        if (debug):
            print(avg_vis.shape)
            print(nunit, nlen.shape, nlen)
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            for dbi in range(3):
                dbpng = 'debug%d.png' % dbi
                x = list(range(nunit))
                fig, s2d = plt.subplots(na-1, na-1, sharex=True, sharey=False, figsize=(12,9))
                b = -1
                for ai in range(na-1):
                    for aj in range(ai+1, na):
                        b += 1
                        if (dbi == 0):
                            s2d[ai,aj-1].plot(x, scatter2[0,b,:]/scatter[0,b,:], color='C2', marker='v', linestyle='', label='lsb')
                            s2d[ai,aj-1].plot(x, scatter2[1,b,:]/scatter[1,b,:], color='C3', marker='^', linestyle='', label='usb')
                        elif (dbi == 1):
                            s2d[ai,aj-1].plot(x, scatter[0,b,:], color='C0', marker='+', linestyle='', label='scatter1')
                            s2d[ai,aj-1].plot(x, scatter2[0,b,:], color='C1', marker='x', linestyle='', label='scatter2')
                        elif (dbi == 2):
                            s2d[ai,aj-1].plot(x, scatter[1,b,:], color='C0', marker='+', linestyle='', label='scatter1')
                            s2d[ai,aj-1].plot(x, scatter2[1,b,:], color='C1', marker='x', linestyle='', label='scatter2')
                            
                s2d[na-2, na-2].legend()
                if (dbi==0):
                    fig.suptitle('scatter2/scatter1')
                elif (dbi==1):
                    fig.suptitle('scatter1 vs scatter2 lsb')
                elif (dbi==2):
                    fig.suptitle('scatter1 vs scatter2 usb')
                fig.tight_layout()
                fig.savefig(dbpng)
                plt.close(fig)
            sys.exit()
        # debug done

    if (do_scatter):
        return int_vis, int_var, int_flag, scatter
    else:
        return int_vis, int_var


def vecAvg(time, cross, pnt, cliptype=1):
    '''
    vector avg the visibility data in time according to the pnt units
    input:
        time[ndata]
        cross[nsb, nb, nch, ndata], complex
        pnt[ncol, nunit]

        pnt[0:] = sid, uid, t1, t2, sod, ra, dec, az, el, hp, op, sp, flag, target
        t1, t2: begin/end time of the unit in epoch time
        sod:    median time of the unit in second-of-day

    average and variance are estimated after a sigma_clipping to avoid spikes and sudden change in signal
    (e.g. caused by timestamp error)

    return:
        avg_time[nunit]                         averaged time (epoch time) of the unit
        avg_vis[nsb, nb, nch, nunit]            complex, integrated visibility
        avg_var[nsb, nb, nch, nunit]            complex, variance of real and imag parts, after the integration
        nlen[nunit]                             number of data points in this unit
        nkeep[nsb, nb, nch, nunit]              number of valid points in this unit
        nclip[nsb, nb, nch, nunit]              number of valid points in this unit
        flag[nsb, nb, nch, nunit]               vecAvg flagging; 0=good, 1=bad
    '''

    nunit = pnt[0].size
    pnt = pnt.reshape((-1, nunit))

    avg_time  = np.zeros(nunit)
    avg_vis   = np.ma.zeros((nsb, nb, nch, nunit), dtype=complex)
    avg_var   = np.ma.zeros((nsb, nb, nch, nunit), dtype=complex)

    nlen    = np.zeros(nunit, dtype=int)
    nkeep   = np.zeros((nsb, nb, nch, nunit), dtype=int)
    nclip   = np.zeros((nsb, nb, nch, nunit), dtype=int)

    flag    = np.zeros((nsb, nb, nch, nunit), dtype=int)    # 0=good, 1=bad

    for ui in range(nunit):
        tw = np.logical_and(time>=pnt[2][ui], time<=pnt[3][ui])
        nlen[ui] = np.count_nonzero(tw)

        avg_time[ui] = time[tw].mean()
        if (cliptype==1):
            cvr = sigma_clip(cross[:,:,:,tw].real, axis=3)
            cvi = sigma_clip(cross[:,:,:,tw].imag, axis=3)
            cmask = cvr.mask + cvi.mask
            nclip[:,:,:,ui] = np.count_nonzero(cmask, axis=3)
            nkeep[:,:,:,ui] = np.count_nonzero(~cmask, axis=3)

            cvr.mask = cmask
            cvi.mask = cmask
            avg_vis[:,:,:,ui] = cvr.mean(axis=3) + 1.j*cvi.mean(axis=3)
            avg_var[:,:,:,ui] = cvr.var(axis=3)/nkeep[:,:,:,ui] + 1.j*cvi.var(axis=3)/nkeep[:,:,:,ui]
            avg_var[:,:,:,ui].mask += (nkeep[:,:,:,ui] <= 1)
            avg_vis[:,:,:,ui].mask = avg_var[:,:,:,ui].mask
            flag[:,:,:,ui] = avg_var[:,:,:,ui].mask.astype(int) # 1=bad; 0=good

        elif (cliptype==2):
            for s in range(nsb):
                for b in range(nb):
                    for c in range(nch):
                        vr = cross[s,b,c,tw].real
                        vi = cross[s,b,c,tw].imag

                        cvr = sigma_clip(vr)
                        cvi = sigma_clip(vr)
                        cmask = cvr.mask + cvi.mask
                        nclip[s,b,c,ui] = np.count_nonzero(cmask)
                        nkeep[s,b,c,ui] = np.count_nonzero(~cmask)

                        cvr.mask = cmask
                        cvi.mask = cmask
                        v2  = cvr + 1j*cvi
                        if (nkeep[s,b,c,ui] > 1):
                            avg_vis[s,b,c,ui] = v2.mean()
                            # avg_var should reflect the variance after integration
                            avg_var[s,b,c,ui] = (np.var(v2.real) + 1.j * np.var(v2.imag))/ nkeep[s,b,c,ui]
                        else:
                            avg_vis[s,b,c,ui] = 0.
                            avg_var[s,b,c,ui] = 0.
                            flag[s,b,c,ui]    = 1


    return avg_time, avg_vis, avg_var, nlen, nkeep, nclip, flag


def unitAvg(time, auto, utime):
    '''
    average the auto-corr according to the time and pnt info
    input:
        time[n]
        auto[nsb,na,nch,n]
        utime[m,2]  <-- a list of (t_begin, t_end) to avg in between
                        e.g. utime = zip(pnt[2,:], pnt[3,:])

    output:
        avg_auto[nsb,na,nch,m]
    '''
    nunit = len(utime)
    sh0 = auto.shape    # input shape
    sh1 = (sh0[0], sh0[1], sh0[2], nunit)
    avg_auto = np.zeros(sh1)
    #print 'sh0, sh1', sh0, sh1

    for ui in range(nunit):
        tw = np.logical_and(time>=utime[ui][0], time<=utime[ui][1])
        #print 'debug: ui, uauto.shape', ui, auto[:,:,:,tw].shape
        avg_auto[:,:,:,ui] = auto[:,:,:,tw].mean(axis=3)

    return avg_auto


def auto2cal(avg_auto):
    '''
    take the time-avg auto-corr and return the geometric avg for each baseline

    input:
        avg_auto[nsb, na, nch, nunit]

    output:
        autocal[nsb, nb, nch, nunit]
    '''

    sh0 = avg_auto.shape
    nunit= sh0[-1]
    sh1 = (nsb, nb, nch, nunit)
    autocal = np.zeros(sh1)

    b = -1
    for ai in range(na-1):
        for aj in range(ai+1, na):
            b += 1
            autocal[:,b,:,:] = np.sqrt(avg_auto[:,ai,:,:] * avg_auto[:,aj,:,:])

    return autocal


def wtoneh5(h5name, time, auto, cross):
    # write out the full dataset into an All-in-One h5 file
    print('... ', datetime.now().isoformat())

    s = auto.shape
    #nsb        = s[0]  # already determined
    na  = s[1]
    # nch       = s[2]  # already determined
    ndata       = s[3]

    nb = int(na * (na-1) / 2)

    #with h5py.File(h5name, 'w') as f:
    f = h5py.File(h5name, 'w')

    f.attrs['na']       = na
    f.attrs['nb']       = nb
    f.attrs['nsb']      = nsb
    f.attrs['nch']      = nch
    f.attrs['ndata']    = ndata

    print('...  timestamp')
    f.create_dataset('timestamp', data = time)
    print('...  auto-corr (real)')
    f.create_dataset('auto',      data = auto)
    print('...  cross-corr (complex)')
    f.create_dataset('cross',     data = cross)

    f.close()

    print('... ', datetime.now().isoformat())


def ldraw(h5name):
    '''
    load the raw (uncorrected) correlator data

    return: time, auto, cross
    '''
    # load the dataset from All-in-one h5 file
    print('... ', datetime.now().isoformat())

    with h5py.File(h5name, 'r') as f:
            #na         = f.attrs['na']
            #nb         = f.attrs['nb']
            #nsb        = f.attrs['nsb']
            #nch        = f.attrs['nch']
            #ndata      = f.attrs['ndata']

            print('...  timestamp')
            time = np.array(f.get('timestamp'))
            print('...  auto-corr (real)')
            auto  = np.array(f.get('auto'))
            print('...  cross-corr (complex)')
            cross = np.array(f.get('cross'))

    print('... ', datetime.now().isoformat())
    return time, auto, cross


def wtrawvis(fname, avg_time, avg_cross, var_cross):
    '''
    write the raw, time_averaged visibility into a file
    '''

    try:
        f = h5py.File(fname, 'w')
        f.create_dataset('time',     data=avg_time)
        f.create_dataset('cross',    data=avg_cross)
        f.create_dataset('variance', data=var_cross)

        f.attrs['cal_table'] = 'raw'

        f.close()
    except:
        print('error writing raw vis into:', fname)
        return 0

    return 1


def ldmeanbp(pp, loghz, body='jupiter', bwmhz=2240.):
    '''
    return the pol mean bandpass (complex) at specified LO_GHz from body

    pol: str, 'X' or 'Y'
    loghz: float, LO freq in GHz
    bwmhz: float, BW in MHz
    '''
    try:
        LOstr = '%05.2fGHz' % loghz
    except:
        print('error reading LO_GHz (float):', loghz)
        sys.exit()

    if (bwmhz == 2240.):
        BWstr = ''
    else:
        BWstr = 'BW%.0f_' % bwmhz
        if (bwmhz == 1920.):
            body = 'mars'   # no jupiter for 1920MHz

    pol = pp.strip('ytla7') # may or may not contain it

    sbody = solbody(body)   # standardize the planet name
    dname = '%s%s_%s' % (BWstr, sbody, LOstr)
    mbfile = '%s/bandpass_mean.ytla7%s.h5' % (bindir, pol)

    proceed = False
    fh = h5py.File(mbfile, 'r')
    if (fh.get(dname)):
        proceed = True
        fh.close()
    else:
        print('no mean bandpass in file?', dname, mbfile)
        sys.exit()

    if (proceed):
        avgcal = getData(mbfile, dname)
        flag   = avgcal.mask.astype(int)
        return avgcal, flag


def ldrawvis(fname, use_acorr=False, rmflags=[], xrm=False):
    '''
    load raw visibility, the one after time averaging

    rmflags: a list of integer bit masks to be ignored in the flag array
             will be removed from the visibility masks

    xrm: F-engine xtalk removal (from raw vis)
         applicable to data before June/2020

    return: avg_time, avg_cross, var_cross, flag
    '''

    try:
        f = h5py.File(fname, 'r')
        avg_time  = f.get('time')[()]
        avg_cross = f.get('cross')[()]
        var_cross = f.get('variance')[()]
        flag      = f.get('flag')[()]

        f.close()
    except:
        print('error reading raw vis from:', fname)
        return None

    if (xrm):
        if ('ytla7X' in fname):
            pname = 'ytla7X'
        elif ('ytla7Y' in fname):
            pname = 'ytla7Y'
        fxtalk = '%s/xtalk_raw.%s.h5' % (bindir, pname)
        xtalk = getData(fxtalk, 'spec')
        xtalk = xtalk.reshape(nsb, nb, nch, 1)
        xsmooth16 = getData(fxtalk, 'smooth16')
        xsmooth16 = xsmooth16.reshape(nsb, nb, nch, 1)
        avg_cross -= xsmooth16

    #tm_flag = np.median(flag, axis=(0,1,2))
    #avg_time = np.ma.array(avg_time, mask=tm_flag)
    avg_time = np.ma.array(avg_time, mask=np.zeros_like(avg_time, dtype=bool))  # try not to flag anything in time


    nflag = flag.copy() # a new flag
    for bi in rmflags:
        nflag -= np.bitwise_and(nflag, int(bi))


    if (use_acorr):
        try:
            dname = 'calibration/acorr'
            f = h5py.File(fname, 'r')
            cal = f.get(dname)[()]
            f.close()
        except:
            print('error reading calibration/acorr')
            return None
        avg_cross /= cal
        var_cross /= cal**2

    avg_cross = np.ma.array(avg_cross, mask=nflag)
    var_cross = np.ma.array(var_cross, mask=nflag)

    return avg_time, avg_cross, var_cross, nflag


def ldcalvis(fname, table='primary', use_acorr=False, rmflags=[], xrm=False):
    '''
    load calibrated visibility

    table: select the calibration table to apply
        examples:
        raw             no correction, same as ldrawvis
        deformation     phase correction from platform deformation model
        primary         the designated primary calibration table (attrs['cal_table'])

        anything else will be loaded from the dataset 'calibration/<table>'

        if the table specified indicates it used acorr corrected data,
        use_acorr is set

    return: avg_time, avg_cross, var_cross, flag
    '''
    #print 'debug: in ldcalvis, table=', table

    try:
        f = h5py.File(fname, 'r')
        g = f.get('calibration')
    except:
        print('error accessing file:', fname)
        return None

    # obtain the primary table(s)
    ctable = f.attrs.get('cal_table')
    if (isinstance(ctable, (list,np.ndarray))):
        ptable = ctable[:]
        pncals = len(ptable)
    else:
        ptable = [ctable]
        pncals = 1
    

    for i in range(pncals):
        x = ptable[i]
        if (isinstance(x, bytes)):
            ptable[i] = x.decode()

    if (isinstance(table, (list,np.ndarray))):
        ncals = len(table)
        if ('primary' in table):
            table.remove('primary')
            table.extend(ptable)
            ncals += pncals - 1
    else:
        if (table.lower() == 'primary'):
            table = ptable
            ncals = pncals
        else:
            ctable = table.lower()
            table = [ctable]
            ncals = 1

    if (not use_acorr): # use_acorr == False
        for ctable in table:
            part = ctable.split('.')
            if (len(part) < 2):
                continue    # do not change use_acorr
            else:
                tail = part[-1]
                if (tail.startswith('a')):
                    use_acorr = True    # change if and only if this is true

    try:
        avg_time, avg_cross, var_cross, flag = ldrawvis(fname, use_acorr=use_acorr, rmflags=rmflags, xrm=xrm)
    except:
        print('error accessing the file:', fname)
        return None


    (nsb, nb, nch, nunit) = avg_cross.shape
    #calibration = np.ones_like(avg_cross, dtype=complex)
    calibration = np.ma.ones(avg_cross.shape, dtype=complex)
    calibration.mask = np.zeros(avg_cross.shape, dtype=bool)

    for ctable in table:
        #print 'debug: in ldcalvis, loop table=', ctable
        if (ctable == 'raw'):
            calibration *= 1.

        else:
            if (g.get(ctable)):
                #cal = getData(f, 'calibration/%s' % ctable)
                cal = g.get(ctable)[()]
                if (g.get('%s.mask' % ctable)):
                    m = g['%s.mask' % ctable][()]
                    cal = np.ma.array(cal, mask=m)
            else:
                print('error finding the cal_table:', ctable)
                return None
            cal = cal.reshape((nsb, nb, nch, -1))
            calibration *= cal

    f.close()


    cnum = calibration.shape[3]
    # reserved for time-varying calibrations
    # for example, in planet primary calibration, usually cnum = 1
    # and in pointing offset correction, cnum = nunit

    if (cnum == 1 or cnum == nunit):
        #avg_cross /= calibration
        #var_cross /= np.abs(calibration)**2
        avg_cross = np.ma.array(avg_cross.data/calibration, mask=avg_cross.mask)
        var_cross = np.ma.array(var_cross.data/np.abs(calibration)**2, mask=var_cross.mask)

    else:
        print('error applying calibration.')
        return None

    return avg_time, avg_cross, var_cross, flag


def manualFlag(antcfg='ant_flag.config', corrcfg='corr_flag.config'):
    '''
    read the manual flagging from antcfg or corrcfg if they exist
    return the manual flag as a dict of mask array: mask[pp].shape=(nsb,nb,nch)
    pp is one of ['X', 'Y']

    if either of antcfg or corrcfg is None, disable that manual flagging
    '''
    shape0 = (nsb, nb, nch)
    print(shape0)

    mask = {}
    mask['X'] = np.zeros(shape0, dtype=bool) # init a new mask array
    mask['Y'] = np.zeros(shape0, dtype=bool) # init a new mask array

    bld = bldict()

    if (antcfg is None or not os.path.isfile(antcfg)):
        print('manual flag: no ant_flag.')
    else:
        print('manual flag: using ant_flag.')
        goodAnt = np.loadtxt(antcfg, dtype=bool)
        print(goodAnt)

        b = -1
        for ai in range(na-1):
            for aj in range(ai+1, na):
                b += 1
                if (goodAnt[ai]==False or goodAnt[aj]==False):
                    #print 'man flag:', ai, aj
                    mask['X'][:,b,:] = True
                    mask['Y'][:,b,:] = True

    if (corrcfg is None or not os.path.isfile(corrcfg)):
        print('manual flag: no corr_flag.')
    else:
        print('manual flag: using corr_flag.')

        blflag = pd.read_csv(corrcfg, comment='#', sep='\s+', names=['PP', 'I', 'J'])
        for pp in ['X', 'Y']:
            conditions = blflag.PP == '%s%s' % (pp,pp)
            print('pp = ', pp)
            print(conditions.tolist())
            for i in range(blflag.PP.size):
                if (conditions[i] == False):
                    continue

                ai = blflag.I[i]
                if (blflag.J[i] == -1): # ant_j is any antenna
                    ajs = list(range(na))
                    ajs.remove(ai)      # except ant_i
                else:
                    ajs = [blflag.J[i]] # ant_j is specific

                for aj in ajs:
                    if (aj > ai):
                        ai2 = ai
                        aj2 = aj
                    else:
                        ai2 = aj
                        aj2 = ai
                    corr = '%d%d' % (ai2, aj2)
                    b = bld[corr]
                    print(blflag.PP[i], corr, '-->', b)
                    mask[pp][:,b,:] = True
                    mask[pp][:,b,:] = True

    return mask


def applyFlag(flag, manflag, pp=None, fname=None):
    '''
    apply the manualFlag to existing flag:
    - if flag is a bit mask, remove original manual flag (16) and apply the new one
    - if flag is a boolean mask (True = flagged), add the new flag
    - if flag is None, create a new flag with shape (nsb, nb, nch, 1)
    
    return updated flag and flag16
    '''
    if (pp is None and fname is None):
        sys.exit('error in applyFlag: need to specify either pp or fname')

    else:
        if (pp is None):    # need to determine pp from fname
            if ('ytla7X' in fname):
                pp = 'X'
            elif ('ytla7Y' in fname):
                pp = 'Y'
            else:
                sys.exit('error in applyFlag: can not determine pp from %s' % fname)
        elif (pp == 'X' or pp == 'Y'):  # the only valid pp
            pass
        else:
            sys.exit("error in applyFlag: pp can only be 'X' or 'Y'")

    # pp is known
    mmask = manflag[pp]
    # a boolean array with shape (nsb, nb, nch)

    # input flag
    if (flag is None):      # no input flag
        sh = (nsb, nb, nch, 1)
        flag16 = mmask.reshape(sh) * 16
        flag = flag16.copy()
        return flag, flag16

    else:                   # input is provided
        sh = flag.shape
        if (len(sh) != 4):
            print('inconsistent shape (nsb, nb, nch, nunit) <-->', sh)
            sys.exit('error in applyFlag: unexpected input flag shape')

        flag16 = np.ones(sh, dtype=int) * 16
        flag16 *= mmask.reshape((nsb, nb, nch, 1))

        if (flag.dtype == int):     # is a bit mask
            flag -= np.bitwise_and(flag, 16)    # remove old flag16
            flag += flag16

        elif (flag.dtype == bool):  # is a boolean mask
            nflag = flag.astype(int) + flag16
            flag = nflag.astype(bool)

        else:
            sys.exit('error in applyFlag: invalud flag type')

        return flag, flag16


def adoneh5(h5name, darray, dname, over=True):
    '''
    set over = False if you do not want to overwrite existing dataset
    '''
    #print datetime.now().isoformat()

    with h5py.File(h5name, 'a') as f:
        print('... ', dname)
        if (f.get(dname)):
            if (over):
                del f[dname]
            else:
                return      # do nothing if over==False
        f.create_dataset(dname, data = darray)

        if (isinstance(darray, np.ma.MaskedArray)):
            mname = '%s.mask' % dname
            if (f.get(mname)):
                del f[mname]
            f.create_dataset(mname, data = darray.mask)

    #print datetime.now().isoformat()


def newraw(rawh5, na):
    # when rawh5 (i.e. .raw.oneh5) does not exist, create a new one from .timestamp file
    # update: can create a new one from basename (without .timestamp)

    if (rawh5.endswith('.raw.oneh5')):
            base = rawh5.rstrip('.raw.oneh5')
    else:
            print('expecting rawh5 file ends with ".raw.oneh5"')
            return None

    (time, auto, cross) = ldcorr(base, na)

    print('output to %s' % rawh5)
    wtoneh5(rawh5, time, auto, cross)


def getLO(fname):
    '''return: LO_MHz'''
    lomhz = 0.

    try:
        f = h5py.File(fname, 'r')
        if (f.attrs.get('LO')):
            lomhz = f.attrs.get('LO')
        elif (f.attrs.get('LOFreq')):
            lomhz = f.attrs.get('LOFreq')
        else:
            print('LO info not found!')
            return None
        if (isinstance(lomhz, (list, np.ndarray))):
            lomhz = lomhz[0]
        #print 'LO =', lomhz, 'MHz'
        f.close()
    except:
        print('error accessing file:', fname)
        return None

    return lomhz


def getBW(fname):
    ''' retrieve the bw info from h5 file
    return default bw (defined in analysisconf.py) if not defined
    '''
    attrs = getAttrs(fname)
    if ('bw' in attrs):
        bw = attrs['bw']
        if (isinstance(bw, (list, np.ndarray))):
            bw = bw[0]
    else:
        bw = bw0

    return bw


def putLO(fname, lomhz):
    if (os.path.isfile(fname)):
        ha = h5py.File(fname, 'a')
        ha.attrs['LO'] = lomhz
        ha.close()
    else:
        print('putLO:: error opening ', fname)
        sys.exit()

    return None


def getPointing(fname):
    '''
    return: pnt_array, pnt_hdr
    '''
    try:
        f = h5py.File(fname, 'r')
        pnt     = f.get('pointing')[()]
        pnt_hdr = f.get('pnt_header')[()]

        f.close()
    except:
        print('getPointing error.')
        return None

    return pnt, pnt_hdr


def getCalInfo(fname, table='primary'):
    '''return: cal_az, cal_el, cal_hp, cal_sp'''
    try:
        f = h5py.File(fname, 'r')
        if (table == 'primary'):
            ctable = f.attrs.get('cal_table')
        else:
            ctable = table
        d = f.get('calibration/%s' % ctable)
        cal_az = d.attrs.get('cal_az')
        cal_el = d.attrs.get('cal_el')
        cal_hp = d.attrs.get('cal_hp')
        cal_sp = d.attrs.get('cal_sp')

        f.close()
    except:
        print('getCalInfo error.')
        return None

    return cal_az, cal_el, cal_hp, cal_sp


def getDeformation(fname):
    '''return: dz, tilt'''
    try:
        f = h5py.File(fname, 'r')
        dz   = f.get('deform_dz')[()]
        tilt = f.get('deform_tilt')[()]

        f.close()
    except:
        print('getDeformation error.')
        return None

    return dz, tilt


def putDeformation(fname, dz, tilt):

    try:
        f = h5py.File(fname, 'r')
        f.create_dataset('deform_dz',   data=dz)
        f.create_dataset('deform_tilt', data=tilt)

        f.close()
    except:
        print('putDeformation error.')
        return None

    return 1


def getIntLen(fname):
    '''return: intLen, numInt'''
    try:
        f = h5py.File(fname, 'r')
        intLen = f.attrs.get('intLen')[0]
        numInt = f.attrs.get('NumInt')[0]
        f.close()
    except:
        print('getIntLen error.')
        return None

    return intLen, numInt


def putPntDef(fname, pnt, param_ver=param_ver0, ant_conf=ant_conf0, phaseCorr=True):
    '''
    saves pointing info and the deformation info into the .oneh5 file
    deformation is rederived based on the pointing info
    also saves a calibration table for phase correction, if enabled

    <input>
    fname:  .oneh5 for data to be written
    pnt:    array of pointing info
            sid, uid, t1, t2, sod, ra, dec, az, el, hp, op, sp, osf
            where
            sid         schedule id (within the obs. date)
            uid         unit id     (within the schedule)
            t1, t2      begin and end time of the unit relative to the begin of uid=1
            sod         avg. time of the unit in second-of-day
            ra          avg. RA in hr
            dec         avg. DEC in deg
            az          avg. Az in deg
            el          avg. El in deg
            hp          avg. HexPol in deg
            op          avg. ObsPol in deg
            sp          avg. SkyPol in deg
            osf         on-source-flag (1=On; 0=Off)
            ha          avg. Hour Angle in hr
    <optional>
    param_ver:  deformation model version (default %s)
    ant_conf:   antenna config (default %s)

    <return>
    1:  success
    0:  fail    
    ''' % (param_ver0, ant_conf0)


    ncol = len(pnt)
    print('ncol=', ncol)
    if (ncol < 14):
        print('not enough info in .pointing file?')
        return 0

    pnt_hdr0 = ['sid', 'uid', 't1', 't2', 'sod', 'ra', 'dec', 'az', 'el', 'hp', 'op', 'sp', 'osf', 'ha']
    pnt_hdr = [x.encode() for x in pnt_hdr0]

    nunit = pnt[0].size

    antxy  = antpos(ant_conf)
    #print('antxy', antxy)
    pntarr = list(zip(pnt[7], pnt[8], pnt[9]))
    #print('pntarr', pntarr)
    dz, tilt = deform(pntarr, antxy, param_ver)
    #print('dz', dz, 'tilt', tilt)

    try:
        f = h5py.File(fname, 'a')
        dpnt = f.create_dataset('pointing', data = pnt[:14])
        #print('pointing' in f)
        dhdr = f.create_dataset('pnt_header', data = pnt_hdr)
        #print('pnt_header' in f)
        dpnt.attrs['pnt_header'] = pnt_hdr

        ddz  = f.create_dataset('deform_dz', data = dz)
        dtil = f.create_dataset('deform_tilt', data = tilt)
        ddz.attrs['param_ver']  = param_ver
        ddz.attrs['ant_conf']   = ant_conf
        ddz.attrs['antxy']      = antxy
        f.close()

    except:
        print('error writing pointing and deformation into h5 file', fname)
        return 0


    if (phaseCorr):
        lomhz = getLO(fname)
        ch = np.array(list(range(nch)))
        freq = np.zeros((nsb, nch))
        # channel frequency in MHz
        freq[0] = 84000. + lomhz - bw0 * (ch + 1.) / float(nch) # lsb
        freq[1] = 84000. + lomhz + bw0 * (ch + 1.) / float(nch) # usb

        # conversion from delay(mm) to phase(rad) 
        scale = 2. * np.pi / 2.998e5 * freq     # lambda_mm = 2.998e11 mms-1/ (freq * 1e6 s-1)

        calibration = np.zeros((nsb, nb, nch, nunit), dtype=complex)

        b = -1
        for ai in range(na-1):
            for aj in range(ai+1, na):
                b += 1

                for s in range(nsb):
                    for ui in range(nunit):
                        dp = (dz[ui,aj] - dz[ui,ai]) * scale[s]
                        calibration[s, b, :, ui] = np.exp(1.j * dp[ui])

        f = h5py.File(fname, 'a')
        dname = 'calibration/deformation'
        if (f.get(dname)):
            f[dname] = calibration
        else:
            f.create_dataset(dname, data=calibration)
        ## optional
        #f.attrs['cal_table'] = 'deformation'
        f.close()

        dname = 'calibration/mod.deform'
        adoneh5(fname, calibration.conjugate(), dname)
        #print('added mod.deform?')
        

    return 1


def putTarget(fname, target):
    '''
    write the target info for each unit
    '''
    try:
        f = h5py.File(fname, 'a')
        dtarg = f.create_dataset('target', data=target)
        f.close()
    except:
        print('error saving target info into file', fname)
        sys.exit()


def getTarget(fname):
    '''
    load the target info for each unit
    '''
    try:
        f = h5py.File(fname, 'r')
        if (f.get('target')):
            target = f['target'][()]
        else:
            target = []
        f.close()
    except:
        print('error getting target infor from file:', fname)
        sys.exit()

    return target


def checkLF(flag, ColdRx):
    '''
    flag out baselines based on ColdRx status (determined by loading factors)
    a string 14 digits (e.g. 11111111111111) showing the status of I and Q LF
    the order being I0 Q0 I1 Q1 I2 Q2 ...

    input: flag[nsb, nb, nch, nunit]

    return: modified flag[nsb, nb, nch, nunit]

    note: 0=good, >0=bad (bit mask)
        bit mask:   1   vecAvg, sigma_clip
                    2   ColdRx
    '''
    goodAnt = np.ones(na, dtype='bool')
    for ai in range(na):
        ii = 2*ai
        if (ColdRx[ii:ii+2] != '11'):
            goodAnt[ai] = False

    flagLF = np.zeros_like(flag)

    b = -1
    for ai in range(na-1):
        for aj in range(ai+1, na):
            b += 1

            if (not (goodAnt[ai] and goodAnt[aj])):
                flag[:,b,:,:] = np.bitwise_or(flag[:,b,:,:], 2)
                flagLF[:,b,:,:] += 2

    return flag, flagLF


def getFlag(fname):
    '''
    return the flag[nsb, nb, nch, nunit]

    note: 0=good, >0=bad (bit mask)
        bit mask:   1   vecAvg, sigma_clip
                    2   ColdRx
                    4   planet data flag
                    8   excessive variance after calibration
                    16  use file
    '''
    try:
        f = h5py.File(fname, 'r')
        flag = f.get('flag')
        f.close()
    except:
        print('getFlag error.')
        return None

    return flag

    
def getUnitLen(fname, full=True):
    '''
    from .avgh5, return the observing unit lengths
    '''
    try:
        fh = h5py.File(fname, 'r')
        nlen  = fh.get('num_unit_len')[()]
    except:
        print('error loading obs unit lengths')
        return None

    if (full):
        try:
            nkeep = fh.get('num_unit_keep')[()]
            nclip = fh.get('num_unit_clip')[()]
        except:
            print('error loading obs unit lengths')
            return None
        return nlen, nkeep, nclip
    else:
        return nlen


def vecSmooth(cspec, kern='g', gsigma=6., twidth=32, axis=2, mode=1):
    '''
    input:
        cspec       ndarray or a masked array
                    nominally, this is the complex visibility

        kern        'g' for Gaussian
                    't' for tophat

        gsigma      sigma of the Gaussian smoothing window

        twidth      width of the tophat

        axis        apply smoothing to which axis

        mode        0 for signal.convolve --> fast with n-D array, but fail for maskedArray
                    1 for np.convolve --> works with maskedArray, but 1-D at a time
    '''
    if (kern == 'g'):
        glen    = int(gsigma * 6.)
        win     = signal.gaussian(glen, gsigma) # 1-D array, length = glen
    elif (kern == 't'):
        glen    = int(twidth)
        win     = np.ones(glen)

    if (mode == 0):
        sh0         = np.ones_like(cspec.shape)
        try:
            sh0[axis] = glen
        except:
            print('error smoothing axis', axis)
            sys.exit()
        
        win     = win.reshape(tuple(sh0))
        #win    = win.reshape((1,1,glen))

        # spectrally smoothed quantities
        sspec   = signal.convolve(cspec, win/win.sum(), 'same')
        if (isinstance(cspec, np.ma.MaskedArray)):
            sspec   = np.ma.array(sspec, mask=cspec.mask)

    elif (mode == 1):
        if (axis !=2):
            print('vecSmooth mode=1 can only process array in the format spec[nsb, nb, nch]')
            sys.exit()
        sspec = cspec.copy()
        for sb in range(nsb):
            for b in range(nb):
                sspec[sb,b] = np.ma.convolve(cspec[sb,b], win/win.sum(), 'same')

    return sspec, win


## the module can be used to convert from correlator output format
## to the oneh5 format
if (__name__ == '__main__'):

    inp = sys.argv[0:]
    pg  = inp.pop(0)

    files   = []
    quick   = False         # if True, skip existing .avg file
    use_pnt = True
    tmax    = 0.
    tseg    = 0.


    usage = '''
    program needs one argument

    %s <file_base> [options]

        <file_base> = something like './data/2017_Oct_26_03_35_27.ytla7X'
        if the .h5 files are in the current directory, the path can be omitted

        <file_base> can also ends with '.timestamp', '.pointing', etc.
        the trailing part after the 2nd '.' will be ignored

        options are:

        --quick             # skip existing .avgh5 files
                            # default is to overwrite it

        --no-pointing       # testing data without pointing
        --tmax TMAX         # (with --no-pointing) set max time-range to use
        --tseg TSEG         # (with --no-pointing) set the segment time

    ''' % pg

    if (len(inp) < 1):
        print(usage)     
        sys.exit()

    while (inp):
        k = inp.pop(0)
        if (k == '--no-pointing'):
            use_pnt = False
        elif (k == '--quick'):
            quick = True
        elif (k == '--tmax'):
            tmax = inp.pop(0)
            try:
                tmax = float(tmax)
                print('tmax set:', tmax)
            except:
                print('error decoding TMAX.')
                sys.exit()
        elif (k == '--tseg'):
            tseg = inp.pop(0)
            try:
                tseg = float(tseg)
                print('tseg set:', tseg)
            except:
                print('error decoding TSEG.')
                sys.exit()
        elif (k.startswith('-')):
            sys.exit('unknown option: %s' % k)
        else:
            #fname = k
            files.append(k)


    for fname in files:
        fname = os.path.basename(fname)
        fpart = fname.split('.')
        fname = '.'.join(fpart[:2])
        dtobj = datetime.strptime(fpart[0], '%Y_%b_%d_%H_%M_%S')
        date1 = dtobj.strftime('%Y-%m-%d')
        date2 = dtobj.strftime('%y%m%d')

        ftime = fname + '.timestamp'
        fpnt  = fname + '.pointing'
        flah  = fname + '.lsb.auto.h5'
        fRAW  = fname + '.RAW.h5'


        fout = fname + '.avgh5'
        if (os.path.isfile(fout) and quick):
            print('%s exists, skip...' % fout)
            continue
        else:
            print('proceed: ', fRAW)


        if (date2 < '180605'):
            if (not os.path.isfile(ftime)):
                print('error finding timestamp file', ftime)
                sys.exit()
            if (not os.path.isfile(flah)):
                print('error finding correlator file', flah)     # a representative one
                sys.exit()
            print('reading data ...')
            time, auto, cross = ldcorr(fname, na)
            attrs = getAttrs(flah)
        else:
            if (not os.path.isfile(fRAW)):
                print('error finding .RAW.h5 file', fRAW)
                sys.exit()
            time, auto, cross = ldcorr2(fname)
            attrs = getAttrs(fRAW)


        # load or create pnt
        if (use_pnt):
            if (not os.path.isfile(fpnt)):
                print('error finding pointing file', fpnt)
                print('assume no-pointing instead')
                use_pnt = False
            else:
                pnt = np.loadtxt(fpnt, usecols=list(range(14)), unpack=True, ndmin=2)

        if (use_pnt == False):
            intLen, numInt = getIntLen(fRAW)
            tlast = time[0] + intLen * numInt
            if (tseg == 0.):
                if (tmax > 0.):
                    tmax += 1   # add 1sec buffer
                    tend = min(tlast, time[0]+tmax)
                else:
                    tend = tlast
                pnt = np.array([1, 1, time[0], tend]).reshape((-1,1))
                target = ['unknown'.encode()]
            else:       #(tseg > 0.)
                pnt = []
                tbeg = time[0]
                uid = 0
                while (tbeg < time[-1]):
                    tend = min(time[-1], tbeg+tseg)
                    uid += 1
                    pnt.append([1, uid, tbeg, tend])
                    #print uid, tbeg, tend
                    tbeg = tend
                pnt = np.array(pnt).T
                target = []
                for i in range(len(pnt[0])):
                    target.append('unknown'.encode())

        print('time-averaging ...')
        #avg_time, avg_cross, var_cross = vecAvg(time, cross, pnt)
        avg_time, avg_vis, avg_var, nlen, nkeep, nclip, flag = vecAvg(time, cross, pnt)
        avg_auto = unitAvg(time, auto, list(zip(pnt[2], pnt[3])))
        autocal = auto2cal(avg_auto)


        print('writing data ...')
        wtrawvis(fout, avg_time, avg_vis, avg_var)
        adoneh5(fout, nlen,  'num_unit_len')
        adoneh5(fout, nkeep, 'num_unit_keep')
        adoneh5(fout, nclip, 'num_unit_clip')
        adoneh5(fout, autocal, 'calibration/acorr')
        adoneh5(fout, avg_auto, 'auto')

        print('copying attributes ...')
        putAttrs(fout, attrs)

        print('basic flagging ...')
        if (attrs.get('ColdRx')):
            ColdRx = attrs['ColdRx'].decode()
            flag1, flagLF = checkLF(flag, ColdRx)
            #print('ColdRx', ColdRx)
            #print('flagLF', flagLF)
        else:
            flag1 = flag
            flagLF = np.zeros_like(flag)
        adoneh5(fout, flag1,  'flag')       # primary flag
        adoneh5(fout, flag,   'flag1_vecAvg')    # record of separate flags
        adoneh5(fout, flagLF, 'flag2_ColdRx')    # record of separate flags


        if (use_pnt):
            print('write pointing and deformation ...')
            putPntDef(fout, pnt)    # use default param_ver and ant_conf
            try:
                target = np.loadtxt(fpnt, usecols=(14,), dtype='S', ndmin=1)
                putTarget(fout, target)
            except:
                print('no target info in the pointing file.')
        else:
            adoneh5(fout, pnt, 'pointing')
            pnt_hdr0 = ['sid', 'uid', 't1', 't2']
            pnt_hdr = [x.encode() for x in pnt_hdr0]
            putAttrs(fout, {'pnt_header':pnt_hdr}, 'pointing')
            adoneh5(fout, pnt_hdr, 'pnt_header')
            adoneh5(fout, target, 'target')
            


