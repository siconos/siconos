#!/bin/sh -fx

# specific variables to the test bench set
frequency=$1

#OAR --project siconos
#OAR --name Guitar_51200
#OAR -p cpumarch != 'haswell'
#OAR -l /nodes=1,walltime=48:00:00

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

