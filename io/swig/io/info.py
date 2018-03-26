#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Show information about a Siconos mechanics-IO HDF5 file.
"""

# Lighter imports before command line parsing
from __future__ import print_function
import sys, argparse

parser = argparse.ArgumentParser(
    description = __doc__)
parser.add_argument('file', metavar='input', type=str, nargs='+',
                    help = 'input file(s) (HDF5)')
parser.add_argument('-O','--list-objects', action = 'store_true',
                    help = 'List object names in the file')
parser.add_argument('-C','--list-contactors', action = 'store_true',
                    help = 'List contactor names in the file')
parser.add_argument('-V','--version', action='version',
                    version='@SICONOS_VERSION@')

if __name__=='__main__':
    args = parser.parse_args()

# Heavier imports after command line parsing
from siconos.io.mechanics_hdf5 import MechanicsHdf5
import numpy as np

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

def list_objects(io):
    print ('')
    print ('Objects:')
    print ('')
    print ('{0:>5} {1:>15} {2:>6} {3:>6}'.format(
        'Id','Name','Mass','ToB'))
    print ('{0:->5} {0:->15} {0:->6} {0:->6}'.format(''))
    for name, obj in io.instances().items():
        print ('{0:>5} {1:>15} {2:>6.4g} {3:>6.4g}'.format(
            obj.attrs['id'], name,
            obj.attrs['mass'],
            obj.attrs['time_of_birth']))

def list_contactors(io):
    print ('')
    print ('Contactors:')
    print ('')
    print ('{0:>5} {1:>15} {2:>9} {3:>9}'.format(
        'Id','Name','Type','Primitive'))
    print ('{0:->5} {0:->15} {0:->9} {0:->9}'.format(''))
    for name, obj in io.shapes().items():
        print ('{0:>5} {1:>15} {2:>9} {3:>9}'.format(
            obj.attrs['id'], name,
            obj.attrs['type'],
            obj.attrs['primitive'] if 'primitive' in obj.attrs else ''))

if __name__=='__main__':
    try:
        with MechanicsHdf5(mode='r', io_filename=args.file[0]) as io:
            if io.dynamic_data() is None or len(io.dynamic_data()) == 0:
                print ('Empty simulation found.')
            else:
                print ('')
                print ('Filename: "{0}"'.format(args.file[0]))
                summarize(io)
                if args.list_objects:
                    list_objects(io)
                if args.list_contactors:
                    list_contactors(io)
    except IOError as e:
        print ('Error reading "{0}"'.format(args.file[0]))
        print (e)
