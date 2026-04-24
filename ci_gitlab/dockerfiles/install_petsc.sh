#!/bin/bash

# --- Script used to configure, build and install petsc ---
#
# Warning: this script is usually called by CI jobs inside a Docker container
# with default values set by CI (template or .gitlab-ci.yml)
# 
# Those values might be different from the default ones when this script is executed manualy, on the command line

set -e

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi

: ${PETSC_INSTALL_DIR:?"Please set environment variable PETSC_INSTALL_DIR with the place where petsc must be installed."}

WORK_DIR="${WORK_DIR:=$HOME/build}"
mkdir -p ${WORK_DIR}
cd ${WORK_DIR}
curl -L https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-lite-3.21.4.tar.gz -o petsc-3.21.4.tar.gz
tar -xzf petsc-3.21.4.tar.gz
cd petsc-3.21.4

# Configure (run from source dir, output directed to install prefix)
# All the flags come from the guix definition of petsc.
./configure --prefix="$PETSC_INSTALL_DIR" --with-64-bit-indices --with-mpi=0 \
            --with-openmp=1 --with-blas-lapack=1 --with-superlu=0 --with-debugging=0 \
            --COPTFLAGS="-g -O3" --CXXOPTFLAGS="-g -O3" --FOPTFLAGS="-g -O3"

# Build and install
make -j "$nbprocs" > /dev/null
echo "----> Installing PETSc..."
make install > /dev/null

# Clean up source and build artefacts
cd ${WORK_DIR}
rm petsc-3.21.4.tar.gz
rm -rf petsc-3.21.4
