#!bin/bash

# --- Standard install of siconos (all components, without oce) ---

# Note : this script takes osname (from docker image in gitlab-ci script) as arg.
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}

# Install python packages
python3 -m pip install -U -r $CI_PROJECT_DIR/ci_gitlab/requirements.txt

# Build Siconos
mkdir $CI_PROJECT_DIR/build
cd $CI_PROJECT_DIR/build

# Build Siconos, run tests but do not submit to cdash
#ctest -S ${CI_PROJECT_DIR}/ci_gitlab/ctest_driver_install_siconos.cmake -Dmodel=$CTEST_MODEL -DSICONOS_INSTALL_DIR=${CI_PROJECT_DIR}/install-siconos -DOSNAME=$1 -DNO_SUBMIT=TRUE -V 
cmake $CI_PROJECT_DIR
make -j 4

# Install
make install
