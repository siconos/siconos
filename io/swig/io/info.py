#!/usr/bin/env python

from siconos.io.mechanics_io import Hdf5
import sys, argparse
import numpy as np

fn = sys.argv[1]

def summarize(io):
    spos_data = io.static_data()
    dpos_data = io.dynamic_data()
    cf_data = io.contact_forces_data()
    solv_data = io.solver_data()
    t0 = dpos_data[:, 0].min()
    t1 = dpos_data[:, 0].max()
    times, counts = np.unique(dpos_data[:, 0], return_counts=True)
    print ('Time simulated: {0} to {1} = {2} steps'.format(t0, t1, len(times)))

    cf_times, cf_counts = np.unique(cf_data[:, 0], return_counts=True)
    min_cf = cf_counts.min()

    # Are there times where there are no contact forces?
    if len(np.setdiff1d(times, cf_times, assume_unique=True)) > 0:
        min_cf = 0

    print ('')
    print ('            {0:>10} {1:>10} {2:>10}'.format('Min','Avg','Max'))
    print ('            {0:->10} {1:->10} {2:->10}'.format('','',''))
    print ('Objects:    {0: >10} {1: >10} {2: >10}'
           .format(counts.min(), int(counts.mean()), counts.max()))
    print ('Contacts:   {0: >10} {1: >10} {2: >10}'
           .format(min_cf, int(cf_counts.mean()), cf_counts.max()))
    print ('Iterations: {0: >10} {1: >10} {2: >10}'
           .format(int(solv_data[:,1].min()), int(solv_data[:,1].mean()),
                   int(solv_data[:,1].max())))
    print ('Precision:  {0: >10.3g} {1: >10.3g} {2: >10.3g}'
           .format(solv_data[:,2].min(), solv_data[:,2].mean(),
                   solv_data[:,2].max()))
    print ('Loc. Prec.: {0: >10.3g} {1: >10.3g} {2: >10.3g}'
           .format(solv_data[:,3].min(), solv_data[:,3].mean(),
                   solv_data[:,3].max()))

from timeit import Timer
try:
    with Hdf5(mode='r',io_filename=fn) as io:
        if io.dynamic_data() is None or len(io.dynamic_data()) == 0:
            print 'Empty simulation found.'
        else:
            print ('')
            print ('Filename: "{0}"'.format(fn))
            summarize(io)
except IOError, e:
    print ('Error reading "{0}"'.format(fn))
    print (e)
