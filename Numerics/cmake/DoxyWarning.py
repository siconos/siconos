#!/usr/bin/env python

import sys
import os
from subprocess import check_call

def touch(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)
doxygen = sys.argv[2]
cdir = sys.argv[1]
output_file = sys.argv[-3]
ppfile = os.path.join(sys.argv[1],
                      os.path.splitext(os.path.basename(sys.argv[-1]))[0])
check_call([doxygen, '{0}.config'.format(ppfile)])
touch(output_file)
