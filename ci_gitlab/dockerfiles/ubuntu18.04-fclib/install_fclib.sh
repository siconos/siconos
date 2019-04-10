#!bin/bash
#
# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}

# -- Fclib download --
cd $CI_PROJECT_DIR
git clone https://github.com/FrictionalContactLibrary/fclib.git

# Creates a directory to build libs.
mkdir -p $CI_PROJECT_DIR/build
cd $CI_PROJECT_DIR/build

# Configure and build fclib
cmake $CI_PROJECT_DIR/fclib -DFCLIB_HEADER_ONLY=OFF
make
make install
