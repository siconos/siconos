#!/bin/bash
oarsub -S "./run_fretless.sh 64225280 4096" -p "network_address='luke42'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh  32112640 2048" -p "network_address='luke43'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 16056320 1024" -p "network_address='luke38'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 8028160 512" -p "network_address='luke39'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 4014080 256" -p "network_address='luke43'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 2007040 128" -p "network_address='luke43'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 1003520 64" -p "network_address='luke43'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 501760 32" -p "network_address='luke38'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 250880 16" -p "network_address='luke42'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 125440 8" -p "network_address='luke42'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 62720 4" -p "network_address='luke42'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 31360 2" -p "network_address='luke42'" -t timesharing=perignon,*
oarsub -S "./run_fretless.sh 15680 1" -p "network_address='luke42'" -t timesharing=perignon,*



