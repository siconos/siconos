"""Update license headers

Usage:

- set correct header in instring (to be replaced)
- set correct header in outstring (expected result)

cd <siconos source dir>
python3 update_license_header.py

"""

import os
import fileinput
rootpath = './'
flist = []
for root, directories, files in os.walk(rootpath):
    for name in files:
        flist.append(os.path.join(root, name))

    exclude = ['.png', '.eps', '.hdf5', '.h5',
               '.jpg', '.gz', '.pdf', '.gif',
               '.npz', '.mat', '.zip', '.prt.1', '.dox',
               '.maple', '.sxd', '.sce', '.TXT',
               '.sci', '.stp']
for ext in exclude:
    flist[:] = [d for d in flist if not d.endswith(ext) and '.git' not in d]

for f in flist:
    with open(f) as ff:
        # try:
        print(f)
        ff.readline()
        # except:
        # print(f)


instring = 'Copyright 2020 INRIA'
# instring ='Copyright (C) 2005, 2018 by INRIA'
# instring ='Siconos-Numerics, Copyright INRIA 2005-2015'
# instring = 'Siconos, Copyright INRIA 2005-2016'

outstring = 'Copyright 2024 INRIA'

with fileinput.input(files=flist, inplace=True) as f:
    for line in f:
        print(line.replace(instring, outstring), end='')
