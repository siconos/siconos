#!/bin/bash
oarsub -S "./run2.sh 64225280 4096" -p "network_address='luke42'"
oarsub -S "./run2.sh 32112640 2048" -p "network_address='luke42'"
oarsub -S "./run2.sh 16056320 1024" -p "network_address='luke42'"
oarsub -S "./run2.sh 8028160 512" -p "network_address='luke42'"
oarsub -S "./run2.sh 4014080 256" -p "network_address='luke42'"
oarsub -S "./run2.sh 2007040 128" -p "network_address='luke42'"
oarsub -S "./run2.sh 1003520 64" -p "network_address='luke42'"
oarsub -S "./run2.sh 501760 32" -p "network_address='luke42'"
oarsub -S "./run2.sh 250880 16" -p "network_address='luke42'"
oarsub -S "./run2.sh 125440 8" -p "network_address='luke42'"
oarsub -S "./run2.sh 62720 4" -p "network_address='luke42'"
oarsub -S "./run2.sh 31360 2" -p "network_address='luke42'"
oarsub -S "./run2.sh 15680 1" -p "network_address='luke42'"


# oarsub -S "./run2.sh 64225280 4096" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 32112640 2048" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 16056320 1024" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 8028160 512" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 4014080 256" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 2007040 128" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 1003520 64" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 501760 32" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 250880 16" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 125440 8" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 62720 4" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 31360 2" -p "network_address='luke42'" -t timesharing=perignon,*
# oarsub -S "./run2.sh 15680 1" -p "network_address='luke42'" -t timesharing=perignon,*


