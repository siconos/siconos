#!bin/bash

# --- Script used to configure, build and test siconos software ---
#
# Warning: this script is call by CI jobs with default values set by CI (template or .gitlab-ci.yml)
#
# Those values might be different from the default ones when this script is executed manualy, on the command line
# 
# 
# Usage
#
# Required: 
# > export CI_PROJECT_DIR=<path-to-siconos-repository>
# > export BUILD_MODE=configure, build or test
#
# Optional:
# > export CONF_FILE=some_conf.cmake, default = CI_PROJECT_DIR/config_samples/default.cmake
# > export CTEST_BUILD_MODE=Experimental or Continuous or Nightly
# > export CDASH_SUBMIT=1 or 0
# > export PARALLEL_BUILD=1 or 0. Set to 1 to allow -jN, 0 to restrict to -j1.
#
# Run:
# > sh ctest_siconos.sh
#
#
# - Will execute :
#    - ctest_configure (cmake) if ctest_mode=configure
#    - ctest_build (make) if ctest_mode=build
#    - ctest_test (test) if ctest_mode=test
#    - ctest... all steps if ctest_mode=all
#  Results will be submitted to cdash if CDASH_SUBMIT=1.
# 
# 

: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos' repository (absolute) path."}
: ${BUILD_MODE:?"Please choose build mode among configure, build or test."}

set -e

# set default config file
CONF_FILE="${CONF_FILE:=$CI_PROJECT_DIR/config_samples/siconos_ci_default.cmake}"
# Default build dir, if not set
BUILD_DIR="${BUILD_DIR:=$HOME/build}" 
# Default ctest mode
CTEST_BUILD_MODEL="${CTEST_BUILD_MODEL:=Experimental}"
# Set to 1 to allow -jN, 0 to restrict to -j1.
PARALLEL_BUILD="${PARALLEL_BUILD=:=1}"
# Default: submit to cdash
CDASH_SUBMIT="${CDASH_SUBMIT=:=1}"

# Read conf file from previous step, if any
# The name of the conf. file is required to set CTEST_BUILD_NAME and ensure proper cdash submissions
if [ $BUILD_MODE != "configure" ] && [ $BUILD_MODE != "all" ] && test -f "$BUILD_DIR/options.env"; then
    export CONF_FILE="$(cat $BUILD_DIR/options.env)"
fi

ctest -S ${CI_PROJECT_DIR}/ci_gitlab/ctest_driver_install_siconos.cmake \
     -Dmodel=$CTEST_BUILD_MODEL -DALLOW_PARALLEL_BUILD=$PARALLEL_BUILD -DCDASH_SUBMIT=$CDASH_SUBMIT \
     -DCTEST_MODE=$BUILD_MODE -DUSER_OPTIONS_FILE=$CONF_FILE --output-junit test_results.xml \
     -DCTEST_BINARY_DIRECTORY=$BUILD_DIR -DCTEST_SOURCE_DIRECTORY=$CI_PROJECT_DIR \
    --output-log $BUILD_DIR/siconos-ctest-$BUILD_MODE.log -VV


echo "\n\n============= CTEST Conf ==============\n"
echo "- Ctest mode: ${BUILD_MODE}"
echo "- Options file: ${CONF_FILE}"
echo "- Log file: siconos-ctest-$BUILD_MODE.log"
echo "\n\n=======================================\n"


# Save conf name in a file that can be used in next CI step. This is useful to ensure the same site/build name for CDash between jobs.
echo $CONF_FILE > $BUILD_DIR/options.env # keep the name of options file for next stages
