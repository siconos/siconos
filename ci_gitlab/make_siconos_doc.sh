# --- Config and build of Siconos documentation ---
# 
# Action :
# - fetch and install all packages and tools required to build siconos doc
# - configure siconos project with documentation on
# - build siconos and build doc.
#
# The script is supposed to be called by ci (job in .gitlab-ci.yml)
# and executed on a runner.
#
# It assumes that CI_PROJECT_DIR is set to siconos repository path (which is the case when called from gitlab-ci).
#
# Build directory will be $CI_PROJECT_DIR/build.
# Documentation output path will be $CI_PROJECT_DIR/build/docs/build/html
# Siconos configuration is described in siconos/CI/siconos_docs.cmake file.
# 
# This last path is supposed to be automatically fetched and published to https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos
# thanks to job pages in gitlab-ci.yml.

# Check if CI_PROJECT_DIR is set AND not empty
: ${CI_PROJECT_DIR:?"Please set environment variable CI_PROJECT_DIR with the path to 'siconos' repository (absolute) path."}

DOC_CI_CONFIG_DIR=$CI_PROJECT_DIR/CI/docs
mkdir $CI_PROJECT_DIR/build
cd $CI_PROJECT_DIR/build
#export LANG=C.UTF-8 # Required, else doxy2swig fails!
cmake $CI_PROJECT_DIR -DUSER_OPTIONS_FILE=$CI_PROJECT_DIR/CI/siconos_docs.cmake
make doc -j 4
