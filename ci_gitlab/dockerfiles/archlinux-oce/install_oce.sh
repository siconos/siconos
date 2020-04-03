#!bin/bash
#
# 
# Usage :
#
# export CI_PROJECT_DIR=<PATH TO SICONOS TUTORIALS REPOSITORY> (or any 'reference path')
# sh install_oce clone_oce
# OR
# sh install_oce
#
# Remark : when call from gitlab ci script, the export CI_PROJECT_DIR step is useless.
#
#
# Result :
# - clone oce (from github) in $CI_PROJECT_DIR/build/
# - configure and build oce in $CI_PROJECT_DIR/build/oce-last
# - install oce in CI_PROJECT_DIR/install/oce
# - clone pythonocc (from github) in $CI_PROJECT_DIR/build/
# - configure and build pythonocc in $CI_PROJECT_DIR/build/pythonocc
# - install pythonocc package in $CI_PROJECT_DIR/install/site-packages
# - set PYTHONPATH properly to allow pythonocc usage (--> test python -c 'import OCC')
#
# When called without 'clone_oce' arg, the first three steps are ignored (install of oce).
# In that case pythonocc is built and installed assuming a former installation of oce (e.g. with package manager)
#

# Get number of procs
if  [ -x "$(command -v nproc)" ]; then
   export nbprocs=`nproc --all`  # linux
elif  [ -x "$(command -v sysctl)" ]; then
   export nbprocs=`sysctl -n hw.ncpu` # macos
else
   export nbprocs=2
fi
# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with 'siconos-tutorials' repository (absolute) path."}
# Create directory to install libs.
mkdir -p $CI_PROJECT_DIR/install/

# Build dir
mkdir -p $CI_PROJECT_DIR/build

cd $CI_PROJECT_DIR

echo "Clone last version of oce ..."
git clone https://github.com/tpaviot/oce.git > /dev/null

# --- OCE ---
mkdir $CI_PROJECT_DIR/build/oce-last
cd $CI_PROJECT_DIR/build/oce-last
# Warning : install in 'user' path, that will be transfered between jobs (artifacts)
cmake $CI_PROJECT_DIR/oce -DOCE_INSTALL_PREFIX=$CI_PROJECT_DIR/install/oce  -Wno-deprecated -Wno-dev -DCMAKE_BUILD_TYPE=Release
make -j $nbprocs > /dev/null
echo "----> install oce ..."
make install > /dev/null
# Save path to OCEConfig.cmake, required to configure pythonocc
export OCE_INSTALL=`grep OCEConfig.cmake install_manifest.txt| sed 's/OCEConfig.cmake//g'`


# -- Python occ --
# Clone last pythonocc version.
# We assume it is complient with the installed oce version.
# Maybe we should clone specific tags for oce and pythonocc? 
cd $CI_PROJECT_DIR
git clone https://github.com/tpaviot/pythonocc-core.git  > /dev/null

mkdir $CI_PROJECT_DIR/build/pythonocc
cd $CI_PROJECT_DIR/build/pythonocc

# Requires (in calling script):
# installpath=`python3 -c "import site;print(site.USER_SITE)"`# Unfortunately, this cannot work, artifacts must be
# in CI_PROJECT_DIR ...
export pyocc_installpath=$CI_PROJECT_DIR/install/site-packages
# Mind the OCC at the end of the install path!
cmake $CI_PROJECT_DIR/pythonocc-core -DCMAKE_BUILD_TYPE=Release -Wno-deprecated -DOCE_DIR=$OCE_INSTALL -DPYTHONOCC_INSTALL_DIRECTORY=$pyocc_installpath/OCC
echo "----> install pythonocc ..."

make install -j $nbprocs > /dev/null

cd $CI_PROJECT_DIR
# test ...
export PYTHONPATH=$pyocc_installpath
python3 -c 'import OCC; print(OCC.__file__)'
