#!/bin/sh -x
export SiconosMechanisms_BUILD=/scratch/siconos/Build/Mechanics/Debug/src/mechanisms/
export SiconosMechanisms_DIR=/home/siconos/siconos/Mechanics/src/mechanisms/
export PYTHONPATH="$SiconosMechanisms_BUILD/frontEnd/MBTB:$SiconosMechanisms_BUILD/frontEnd/CADMBTB:"

export CSF_GraphicShr=/usr/lib/libTKOpenGl-6.3.0.so


env
