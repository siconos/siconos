#!bin/bash
#
# Check if CI_PROJECT_DIR is set AND not empty

# -- Fclib download --
cd $HOME
git clone https://github.com/FrictionalContactLibrary/fclib.git

# Creates a directory to build libs.
mkdir -p $HOME/build
cd $HOME/build

# Configure and build fclib
cmake $HOME/fclib -DFCLIB_HEADER_ONLY=OFF
make
make install
