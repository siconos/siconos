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
git clone https://github.com/HDFGroup/hdf5.git
mkdir build
cd build
#sh /opt/intel/oneapi/setvars.sh
export CC=mpiicc
export CXX=mpiicpc
cmake ../hdf5 -DCMAKE_C_COMPILER="mpiicc;-cc=icx" -DCMAKE_CXX_COMPILER="mpiicpc;-cxx=icpx" -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/install/hdf5 -DHDF5_ENABLE_PARALLEL=ON
make -j $nbprocs
make install
