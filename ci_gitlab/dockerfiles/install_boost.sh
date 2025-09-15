#!bin/bash

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}

mkdir -p $CI_PROJECT_DIR/deps
cd $CI_PROJECT_DIR/deps
wget https://archives.boost.io/release/1.81.0/source/boost_1_81_0.tar.bz2
tar -jxvf boost_1_81_0.tar.bz2
cd boost_1_81_0 && ./bootstrap.sh && ./b2 install
