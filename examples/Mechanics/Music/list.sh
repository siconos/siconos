#!/bin/bash
oarsub -S "./run.sh 64225280 4096"
oarsub -S "./run.sh 32112640 2048"
oarsub -S "./run.sh 16056320 1024"
oarsub -S "./run.sh 8028160 512"
oarsub -S "./run.sh 4014080 256"
oarsub -S "./run.sh 2007040 128"
oarsub -S "./run.sh 1003520 64"
oarsub -S "./run.sh 501760 32"
oarsub -S "./run.sh 250880 16"
oarsub -S "./run.sh 125440 8"
oarsub -S "./run.sh 62720 4"
oarsub -S "./run.sh 31360 2"
oarsub -S "./run.sh 15680 1"


