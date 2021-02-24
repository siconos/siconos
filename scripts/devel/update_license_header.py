
import os
import fileinput
rootpath = './sico-2021'
flist = []
for root, directories, files in os.walk(rootpath):
    for name in files:
        flist.append(os.path.join(root, name))

    exclude = ['.png', '.eps', '.hdf5', '.h5','.jpg', '.gz','.pdf', '.npz']
for ext in exclude:
    flist[:] = [d for d in flist if not d.endswith(ext) and '.git' not in d]

# for f in flist:
#     with open(f) as ff:
#         try:
#             ff.readline()
#         except:
#             print(f)


instring ='Copyright 2020 INRIA'
#instring ='Copyright (C) 2005, 2018 by INRIA'
#instring ='Siconos-Numerics, Copyright INRIA 2005-2015'
#instring = 'Siconos, Copyright INRIA 2005-2016'

outstring='Copyright 2021 INRIA'

with fileinput.input(files=flist, inplace=True) as f:
    for line in f:
        print(line.replace(instring, outstring), end='')
