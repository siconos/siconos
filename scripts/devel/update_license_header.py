import os
import sys
from pathlib import Path

root_dir_instance = Path('.')
siconos_files = [item.absolute() for item in root_dir_instance.glob("**/*") if not item.is_dir()]
siconos_files = [item for item in siconos_files if not '.git' in item.as_posix()]
siconos_files = [item for item in siconos_files if not 'update_license' in item.as_posix()]
#siconos_files = [os.path.join(item.parent.name, item.name) for item in root_dir_instance.glob("**/*") if not item.is_dir()]

search = ['Copyright 2016 INRIA', 'Copyright 2018 INRIA', 'Copyright 2019 INRIA']

replace = ['Copyright 2020 INRIA'] * 3

if (len(search) != len(replace)) :
    sys.exit("Error: search does not match with replace")

modified = []
count = 0
# Loop through files and replace "search" with "replace"
for filename in siconos_files:
    #filename = Path(filename).resolve()
    #print(filename)
    # Read in the file

    try:
        with open(filename, 'r') as file :
            filedata = file.read()

        fdata2 = filedata # copy
        for i in range(len(search)):
            #if filedata.find(search[i]) > -1:
            fdata2 = fdata2.replace(search[i], replace[i])

            #modified.append(filename)
               
            # Write the file out again
        #if filename in modified:
        # print(filename)
        if fdata2 == filedata:
            #with open(filename, 'w') as file:
            print("oualou", filename)
                #file.write(fdata2)
        else:
            with open(filename, 'w') as file:
                file.write(fdata2)
            #pass#print("ok")
    except UnicodeDecodeError:
        pass
 
