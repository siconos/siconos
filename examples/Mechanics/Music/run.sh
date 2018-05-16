#!/bin/bash
source /applis/site/nix.sh

frequency=$1
output_freq=$2
restit=$3
final_time=$4
case=$5
#OAR --project siconos
#####OAR --name Guitar_Bass
#OAR -l /nodes=1/core=1,walltime=1240:00:00
### #OAR -t timesharing=perignon,*

siconos_dir=$HOME/Softs/myfork/siconos/examples/Mechanics/Music
rundir=/nfs_scratch/$USER/Music/F_${frequency}_id_${OAR_JOB_ID}
echo $frequency ${OAR_JOB_ID} ${HOSTNAME} $output_freq $restit ${OAR_JOB_NAME}>> jobs_params

mkdir -p $rundir
cd $rundir
# 
siconos $siconos_dir/run.py $frequency $output_freq $restit $final_time $case

