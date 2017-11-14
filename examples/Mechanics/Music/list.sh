#!/bin/bash

oarsub -S "./run.sh 15680 1"
oarsub -S "./run.sh 31360 1"
oarsub -S "./run.sh 62720 1"
oarsub -S "./run.sh 125440 1"
oarsub -S "./run.sh 250880 1"
oarsub -S "./run.sh 501760 1"
oarsub -S "./run.sh 1003520 1"
oarsub -S "./run.sh 2007040 1"
oarsub -S "./run.sh 4014080 1"
oarsub -S "./run.sh 8028160 1"
oarsub -S "./run.sh 16056320 1"
oarsub -S "./run.sh 32112640 1"
oarsub -S "./run.sh 64225280 1"



