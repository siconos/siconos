#!/bin/bash
source /applis/site/nix.sh

frequency=$1
output_freq=$2
#OAR --project siconos
#OAR --name Fretless_Bass
#OAR -l /nodes=1/core=1,walltime=240:00:00
### #OAR -t timesharing=perignon,*

siconos_dir=$HOME/Softs/myfork/siconos/examples/Mechanics/Music
rundir=/nfs_scratch/$USER/Music/Fretless/F_${frequency}_id_${OAR_JOB_ID}

mkdir -p $rundir
cd $rundir
# 
siconos $siconos_dir/run_fretless_bis.py $frequency $output_freq

