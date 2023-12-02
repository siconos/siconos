#!bin/bash

# --- Script used in continuous integration to configure, build and test siconos software ---
#
# Usage :
# > export CI_PROJECT_DIR=<path-to-siconos-repository>
# > export ctest_build_model=Experimental or Continuous or Nightly
# > export cdash_submit=1 or 0
# > export allow_parallel_build=1 or 0. Set to 1 to allow -jN, 0 to restrict to -j1.
# > export BUILD_MODE=Configure, Build or Test
# > sh ctest_siconos.sh
#
# - ctest_mode : choose among 'Configure', 'Build', 'Test' or 'all'
# - user_option_filename is optional. If not set, siconos build will use <siconos repository>/cmake/default_options.cmake file.
#
# - Will execute :
#    - ctest_configure (cmake) if ctest_mode=Configure
#    - ctest_build (make) if ctest_mode=Build
#    - ctest_test (test) if ctest_mode=Test
#    - ctest... all steps if ctest_mode=all
#  Results will be submitted to cdash if cdash_submit=1.
# 
# 
#  Use absolute path or path relative to $CI_PROJECT_DIR/build
#
# The default installation path for siconos is /home/install-siconos.
# Use -DSICONOS_INSTALL_DIR=<something else> as ctest option to change this location.

: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}
: ${ctest_build_model:?"Please set Dashboard client mode. Choose among Experimental, Continuous or Nightly."}
: ${cdash_submit:?"Please set environment variable cdash_submit to TRUE or FALSE. If true, ctests results will be submitted to cdash server."}
: ${allow_parallel_build:?"Please set environment variable allow_parallel_build to TRUE or FALSE. If true, ctests will use paralle build option (-jN)".}
: ${BUILD_MODE:?"Please choose build mode among configure, build or test."}

# set default config file
CONF_FILE="${CONF_FILE:=$CI_PROJECT_DIR/config_samples/siconos_ci_default.cmake}"


echo "${BUILD_MODE} and ${CONF_FILE}"
if [ $BUILD_MODE = "configure" ] || [ $BUILD_MODE = "all" ]
then
    rm -rf $CI_PROJECT_DIR/build
    mkdir -p $CI_PROJECT_DIR/build
fi
   
# --- Run ctest for Siconos ---
cd $CI_PROJECT_DIR/build

ctest -S ${CI_PROJECT_DIR}/ci_gitlab/ctest_driver_install_siconos.cmake -Dmodel=$ctest_build_model -DALLOW_PARALLEL_BUILD=$allow_parallel_build -DCDASH_SUBMIT=$cdash_submit -V -DCTEST_MODE=$BUILD_MODE -DUSER_FILE=$CONF_FILE --output-junit test_results.xml
