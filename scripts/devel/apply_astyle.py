"""Apply some astyle config to all siconos files.

Usage :

python scripts/devel/apply_astyle.py

in the siconos source dir.

"""

from pathlib import Path
import subprocess

astyle_cmd = ['astyle', '--style=ansi', '-U', '-v', '-s2']

currentdir = Path.cwd()

# Get a list of all C files
all_c_files = list(currentdir.glob('**/*.c'))
# Get a list of all CXX files
all_cxx_files = list(currentdir.glob('**/*.cpp'))


# Apply astyle command to all of them

for file in all_c_files:
    cmd = astyle_cmd + [file]
    subprocess.run(cmd)

for file in all_cxx_files:
    cmd = astyle_cmd + [file]
    subprocess.run(cmd)


# clean .orig files
orig_files = list(currentdir.glob('**/*.orig'))
for file in orig_files:
    file.unlink()
                
