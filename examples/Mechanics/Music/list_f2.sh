#!/bin/bash
oarsub -S "./run_fretless_bis.sh 64225280 4096" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 32112640 2048" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 16056320 1024" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 8028160 512" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 4014080 256" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 2007040 128" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 1003520 64" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 501760 32" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 250880 16" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 125440 8" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 62720 4" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 31360 2" -p "network_address='luke44'"
oarsub -S "./run_fretless_bis.sh 15680 1" -p "network_address='luke44'"


# oarsub -S "./run_fretless_bis.sh 64225280 4096" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 32112640 2048" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 16056320 1024" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 8028160 512" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 4014080 256" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 2007040 128" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 1003520 64" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 501760 32" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 250880 16" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 125440 8" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 62720 4" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 31360 2" -p "network_address='luke44'" -t timesharing=perignon,*
# oarsub -S "./run_fretless_bis.sh 15680 1" -p "network_address='luke44'" -t timesharing=perignon,*


