#!/usr/bin/env python

import numpy as np
import h5py
import os, sys, argparse

parser = argparse.ArgumentParser(
    description = 'Copy a Siconos HDF5 simulation file, filtering the contents.')
parser.add_argument('fns_in', metavar='input', type=str, nargs='+',
                    help = 'input file(s) (HDF5)')
parser.add_argument('fn_out', metavar='output', type=str, nargs=1,
                    help = 'output file (HDF5)')
parser.add_argument('--start', metavar='time', type=float,
                    help = 'time in seconds to cut the start of the simulation')
parser.add_argument('--end', metavar='time', type=float,
                    help = 'time in seconds to cut the end of the simulation')
parser.add_argument('--interval', metavar='time', type=float,
                    help = 'time between frames to preserve')
parser.add_argument('--gzip', action = 'store_true',
                    help = 'enable compression in copy')
parser.add_argument('--single', action = 'store_true',
                    help = 'use single-precision floats in copy')

class CopyVisitor(object):
    """The CopyVisitor is called for each group and dataset in the HDF5
       file, and is responsible for copying the structure to the new
       HDF5 file."""
    def __init__(self, time_filter = None):
        self.time_filter = None
        if time_filter is not None:
            self.time_filter = np.vectorize(time_filter)
        self.time_idx = None

    def visitor(self, path, obj):
        gr = io_out
        names = path.split('/')
        if len(names) > 1:
            for i,name in enumerate(names[:-1]):
                if name in gr:
                    gr = gr[name]
                else:
                    gr = gr.create_group(name)
                    gr_in = io_in['/'.join(names[:i+1])]
                    for a in gr_in.attrs:
                        gr.attrs[a] = gr_in.attrs[a]

        # Copy a dataset after creating the groups it is in
        if obj.__class__ == h5py.Dataset:
            chunks = obj.chunks
            if np.any(np.array(chunks) > np.array(obj.shape)):
                chunks = None

            shape = obj.shape
            time_idx = None
            if self.time_filter is not None:
                if path in ['data/cf', 'data/dynamic', 'data/velocities']:
                    if self.time_idx is None or self.times is None:
                        # Get indexes of filtered times in data/dynamic
                            dyn = obj.file['data/dynamic']
                            self.time_idx = self.time_filter(dyn[:,0]).nonzero()[0]
                            self.times = dyn[self.time_idx, 0]
                    # Get indexes of corresponding times in current dataset
                    if path is 'data/dynamic':
                        time_idx = self.time_idx
                    else:
                        time_idx = np.in1d(obj[:,0], self.times).nonzero()[0]
                    # Shape of filtered dataset
                    shape = (len(time_idx),) + tuple(shape[1:])

            # Create the dataset, supply compression and dtype
            # arguments, possibly overridden by command line arguments
            comp = ((obj.compression is True or args.gzip)
                    and chunks is not None)
            if comp:
                chunks = (4000,) + tuple(obj.shape[1:])
            ds = gr.create_dataset(obj.name,
                                   dtype = [obj.dtype,'f4'][args.single],
                                   shape = shape,
                                   maxshape = obj.maxshape,
                                   chunks = chunks,
                                   compression = [obj.compression, True][comp],
                                   compression_opts = [obj.compression_opts,9][comp],
                                   shuffle = comp,
                                   fletcher32 = obj.fletcher32)
            for a in obj.attrs:
                ds.attrs[a] = obj.attrs[a]

            # Copy the filtered or unfiltered dataset
            if time_idx is not None:
                ds[xrange(len(time_idx)),:] = obj[time_idx,:]
            else:
                ds[:] = obj

        # Some groups might be empty so we copy them anyway for their attributes
        elif obj.__class__ == h5py.Group:
            if path in io_out:
                gr = io_out[gr]
            else:
                gr = io_out.create_group(path)
            for a in obj.attrs:
                gr.attrs[a] = obj.attrs[a]

        else:
            print('Unknown type "{0}": {1}'.format(path, str(obj.__class__)))

if __name__ == '__main__':
    args = parser.parse_args()
    if os.path.exists(args.fn_out[0]):
        print('Output file "{0}" already exists!'.format(args.fn_out[0]))
        sys.exit(1)

    class TimeFilter(object):
        def __init__(self):
            self.marker = None
            self.last = None
        def __call__(self, t):
            res = True
            if args.end is not None:
                res = res and t <= args.end
            if args.start is not None:
                res = res and t >= args.start

            # Interval filter is true only for every arg.interval time index
            if res is True and self.marker is None:
                self.marker = t
            if args.interval is not None and self.marker is not None:
                if self.last is None:
                    if t >= self.marker:
                        self.last = t
                        self.marker += args.interval
                    else:
                        res = False
                else:
                    if t > self.last:
                        res = False
                        self.last = None
            return res

    with h5py.File(args.fns_in[0], mode='r') as io_in:
        with h5py.File(args.fn_out[0], mode='w') as io_out:
            io_in.visititems(CopyVisitor(
                time_filter = TimeFilter()
            ).visitor)
