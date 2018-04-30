#!/bin/bash
source /applis/site/nix.sh

frequency=$1
output_freq=$2
restit=$3
#OAR --project siconos
#OAR --name Fretless_Bass_09
#OAR -l /nodes=1/core=1,walltime=1240:00:00
### #OAR -t timesharing=perignon,*

siconos_dir=$HOME/Softs/myfork/siconos/examples/Mechanics/Music
rundir=/nfs_scratch/$USER/Music/F_${frequency}_id_${OAR_JOB_ID}
echo $frequency ${OAR_JOB_ID} ${HOSTNAME} $output_freq $restit ${OAR_JOB_NAME}>> jobs_params
mkdir -p $rundir
cd $rundir
# 
siconos $siconos_dir/run_fretless.py $frequency $output_freq $restit

