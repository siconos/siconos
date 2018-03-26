#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Numerically compare the contents of two Siconos mechanics-IO HDF5 simulation files.

Sole output on standard output shall be the maximum difference between
the columns selected for comparison.  Return code 0 will be returned
if maximum difference is less than the threshold, otherwise 1 will be
returned if greater, or 2 if the files could not be opened or do not
contain the necessary tables.  Non-matching time indices shall be
considered as failing the comparison.
"""

# Lighter imports before command line parsing
from __future__ import print_function
import os, sys, argparse, re

parser = argparse.ArgumentParser(
    description = __doc__)
parser.add_argument('fns_in', metavar='<file>', type=str, nargs=2,
                    help = 'input files (HDF5)')
parser.add_argument('--columns','-c', metavar='<colspec>', type=str,
                    help = """
list of columns to include in the comparison (0-indexed,
comma-separated, may be prefixed with table name, and multiple lists
separated by semicolon, e.g., "dynamic:0,1,2;cf:0,1,2"; default
table="dynamic").  See Siconos HDF5 file format docs for meaning of
columns.""")
parser.add_argument('--start', metavar='<time>', type=float,
                    help = 'start time in seconds from beginning the comparison')
parser.add_argument('--end', metavar='<time>', type=float,
                    help = 'end time in seconds end the comparison')
parser.add_argument('--interval', metavar='<time>', type=float,
                    help = 'amount of time after start to compare')
parser.add_argument('--threshold', metavar='<float>', type=float, default=1e-14,
                    help = 'threshold for comparison')
parser.add_argument('-V','--version', action='version',
                    version='@SICONOS_VERSION@')

if __name__ == '__main__':
    args = parser.parse_args()
    columns = args.columns
    tablenames = set()
    if args.columns is None:
        columns = [('dynamic', None)]
        tablenames.update(['dynamic'])
    else:
        columns = []
        for spec in args.columns.split(';'):
            if ':' in spec:
                table, cols = spec.split(':')
            else:
                table = 'dynamic'
                cols = spec
            if ',' in cols:
                columns += [(table, int(col)) for col in cols.split(',')]
            else:
                columns += [(table, None)]
            tablenames.update([table])

# Heavier imports after command line parsing
import numpy as np
import h5py

def verify_tables(tablenames, columns, io1, io2):
    """Verify files contain tables and tables have same columns."""
    for table in tablenames:
        if not table in io1['data']:
            print('File "{}" does not have table "{}".'.format(
                args.fns_in[0], table), file=sys.stderr)
            sys.exit(2)
        if not table in io2['data']:
            print('File "{}" does not have table "{}".'.format(
                args.fns_in[1], table), file=sys.stderr)
            sys.exit(2)
        for c in columns:
            if c[1] is None:
                continue
            if c[0] == table and c[1] >= io1['data'][table].shape[1]:
                print('Table "{}" in file "{}" does not have specified column {}.'
                      .format(table, args.fns_in[0], c[1]), file=sys.stderr)
                sys.exit(2)
            if c[0] == table and c[1] >= io2['data'][table].shape[1]:
                print('Table "{}" in file "{}" does not have specified column {}.'
                      .format(table, args.fns_in[1], c[1]), file=sys.stderr)
                sys.exit(2)
        if io1['data'][table].shape[1] != io2['data'][table].shape[1]:
            print('Tables "{}" do not have same number of columns in each file.'
                  .format(table), file=sys.stderr)
            sys.exit(2)

def compare_tables(tablenames, columns, io1, io2):
    """Compare tables and columns given."""
    maxdiff = 0.0
    for table in tablenames:
        for c in columns:
            if c[0]!=table:
                continue
            t1 = io1['data'][table]
            t2 = io2['data'][table]

            # start/end/interval times
            S1, S2 = 0, 0
            E1, E2 = t1.shape[0]-1, t2.shape[0]-1
            if args.start is not None:
                S1 = np.searchsorted(t1[:,0], args.start)
                S2 = np.searchsorted(t2[:,0], args.start)
            if args.end is not None:
                E1 = np.searchsorted(t1[:,0], args.end)
                E2 = np.searchsorted(t2[:,0], args.end)
            elif args.interval is not None:
                E1 = np.searchsorted(t1[:,0], args.start + args.interval)
                E2 = np.searchsorted(t2[:,0], args.start + args.interval)
            for t,n,s,T,fn in [(S1,args.start,'Start',t1,args.fns_in[0]),
                               (S2,args.start,'Start',t2,args.fns_in[1]),
                               (E1,args.end,'End',t1,args.fns_in[0]),
                               (E2,args.end,'End',t2,args.fns_in[1])]:
                if t >= T.shape[0]:
                    print('{} time {} beyond the end of table "{}" for file "{}".'
                          .format(s, n, table, fn), file=sys.stderr)
                    sys.exit(2)

            # TODO: we assume same sampling rate for now, later,
            # support linear interpolation?
            if (E1-S1) != (E2-S2):
                print(('Tables "{}" do not have same number of rows in each file '
                       +'for requested range.').format(table), file=sys.stderr)
                sys.exit(2)

            # Compare one column at a time
            cols = [c[1]]
            if c[1] is None:
                cols = range(t1.shape[1])
            for j in cols:
                d = np.abs(t1[S1:E1,j] - t2[S2:E2,j]).max()
                if d > maxdiff:
                    maxdiff = d
    print(maxdiff)
    return int(maxdiff >= args.threshold)

if __name__ == '__main__':
    with h5py.File(args.fns_in[0], mode='r') as io1:
        with h5py.File(args.fns_in[1], mode='r') as io2:
            verify_tables(tablenames, columns, io1, io2)
            sys.exit(compare_tables(tablenames, columns, io1, io2))
