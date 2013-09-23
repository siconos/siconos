#!/bin/csh -x
source ~/.tcshrc
setenv SiconosMechanisms_BUILD /scratch/Vincent/Build/Mechanics/Debug/src/mechanisms/
setenv SiconosMechanisms_DIR  /Users/acary/siconos/Mechanics/src/mechanisms/
setenv PYTHONPATH $PYTHONPATH":$SiconosMechanisms_BUILD/frontEnd/MBTB:$SiconosMechanisms_BUILD/frontEnd/CADMBTB:"

setenv CSF_GraphicShr /Library/OpenCASCADE/6.3.0/lib/libTKOpenGl.dylib
setenv CSF_GraphicShr /usr/local/lib/libTKernel.dylib

env
