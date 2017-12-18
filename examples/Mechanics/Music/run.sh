#!/bin/bash -fx
source /applis/site/nix.sh

# specific variables to the test bench set
frequency=$1

#OAR --project siconos
#OAR --name Guitar_Bass
#OAR -p network_address='luke42'
#OAR -l /nodes=1,walltime=240:00:00
#OAR -t timesharing=acary,*

siconos_dir=$HOME/siconos/


rundir=/nfs_scratch/acary/Music/F=${frequency}_${OAR_JOB_ID}
mkdir -p $rundir
cd $rundir
. $rundir

# 

#rsync -av $siconos_dir/examples/Mechanics/Music/ .
cp $siconos_dir/examples/Mechanics/Music/run.py .
cp $siconos_dir/examples/Mechanics/Music/guitar.py .
cp $siconos_dir/examples/Mechanics/Music/numpywrappers.py .
cp -r $siconos_dir/examples/Mechanics/Music/donnees_siconos .
python run.py $frequency

