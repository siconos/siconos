#!/usr/bin/env @PYTHON_EXECUTABLE@
"""
Description: Export a Siconos mechanics-IO HDF5 file to raw data .dat file
"""

from siconos.io.vview import VView, VRawDataExportOptions
from siconos.io.mechanics_hdf5 import MechanicsHdf5

if __name__=='__main__':
    # Parse command-line
    opts = VRawDataExportOptions()
    opts.parse()

    ## Options and config already loaded above
    with MechanicsHdf5(io_filename=opts.io_filename, mode='r') as io:
        vview = VView(io, opts)
        vview.export_raw_data()
