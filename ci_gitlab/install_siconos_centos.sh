#!bin/bash

# --- Build, test (opt) and install Siconos software ---
#
# Usage :
# > export CI_PROJECT_DIR=<path-to-siconos-repository>
# > export ctest_build_model=Experimental or Continuous or Nightly
# > export IMAGE_NAME="some name for cdash build site"
# 
# > sh install_siconos.sh user_option_filename
#
# will build Siconos using the configuration defined in user_option_filename.
#
# - user_option_filename is optional. If not set, siconos build will use cmake/siconos_default.cmake file.
#   Use absolute path or path relative to $CI_PROJECT_DIR/build
# - export for CI_PROJECT_DIR, IMAGE_NAME and ctest_build_model is not needed when this script is called by gitlab-ci.
#

: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}
: ${ctest_build_model:?"Please set Dashboard client mode (environment variable ctest_build_model). Choose among Experimental, Continuous or Nightly."}
: ${IMAGE_NAME:?"Please set environment variable IMAGE_NAME. It will be used to name cdash build site."}
: ${cdash_submit:?"Please set environment variable cdash_submit to TRUE or FALSE. If true, ctests results will be submitted to cdash server."}

# Create build dir
mkdir -p $CI_PROJECT_DIR/build
cd $CI_PROJECT_DIR/build
#tmp fix
# --- Run ctest for Siconos ---
# configure, build, test and submit to cdash.
# 
# Input variables are :
# - model (from gitlab-ci file), Dashboard client mode can be Continuous, Nightly, Experimental, check https://cmake.org/cmake/help/latest/manual/ctest.1.html#ctest-start-step
# - SICONOS_INSTALL_DIR : where Siconos will be installed
# - USER_FILE : user options file.
# - OSNAME : set to IMAGE_NAME
cmake $CI_PROJECT_DIR -DUSER_OPTIONS_FILE=$CI_PROJECT_DIR/ci_gitlab/siconos_default.cmake -DBOOST_LIBRARYDIR=/usr/lib64/boost169 -DBOOST_INCLUDEDIR=/usr/include/boost169 -DCLAPACK_LIBRARY=/usr/lib64/libopenblas.so -DCMAKE_CXX_STANDARD=11
make -j 4
# Install
#make install
