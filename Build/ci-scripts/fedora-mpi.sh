#!/bin/bash

MODULEPATH=/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles
export MODULEPATH
eval `/usr/bin/modulecmd bash load mpi/openmpi-x86_64`

cmake "$@"
