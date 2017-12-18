#!/bin/sh -fx

oarsub -S "./run.sh 1960"
oarsub -S "./run.sh 3920"
oarsub -S "./run.sh 7840"
oarsub -S "./run.sh 15680"
oarsub -S "./run.sh 31360"
oarsub -S "./run.sh 62720"
oarsub -S "./run.sh 125440"
oarsub -S "./run.sh 250880"
oarsub -S "./run.sh 501760"
oarsub -S "./run.sh 1003520"
oarsub -S "./run.sh 2007040"
oarsub -S "./run.sh 4014080"
oarsub -S "./run.sh 8028160"
oarsub -S "./run.sh 16056320"
oarsub -S "./run.sh 32112640"
oarsub -S "./run.sh 64225280"



