#!/usr/bin/env python

from analysisconf import *
import multiprocessing as mp
from itertools import repeat

def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)

def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)


if (__name__ == '__main__'):
    kwargs = {'shell':True}

    args_iter = []
    kwargs_iter = []
    for i in range(2):
        if (i==0):
            #iarg = ['echo $HOSTNAME']
            iarg = ['ls /pool-frb/cli-command/test/*.meta']
        elif (i==1):
            iarg = ['ssh nuc6 -T "ls /pool-frb/cli-command/test/*.meta"']
        args_iter.append(iarg)
        kwargs_iter.append(kwargs)

    print(args_iter, kwargs_iter)

    with mp.Pool(2) as p:
        #results = p.starmap(call, args)
        results = starmap_with_kwargs(p, call, args_iter, kwargs_iter)

    #for res in results:
    #    print(res)
