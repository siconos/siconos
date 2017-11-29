28/11/2017
----------

# Bass guitar, e=0.9
sh ./list.sh
oarsub -S "./run.sh 64225280 4096" -p "network_address='luke42'"
OAR_JOB_ID=3875402
oarsub -S "./run.sh  32112640 2048" -p "network_address='luke42'"
OAR_JOB_ID=3875403
oarsub -S "./run.sh 16056320 1024" -p "network_address='luke42'"
OAR_JOB_ID=3875404
oarsub -S "./run.sh 8028160 512" -p "network_address='luke42'" 
OAR_JOB_ID=3875405
oarsub -S "./run.sh 4014080 256" -p "network_address='luke42'"
OAR_JOB_ID=3875406
oarsub -S "./run.sh 2007040 128" -p "network_address='luke42'"
OAR_JOB_ID=3875407
oarsub -S "./run.sh 1003520 64" -p "network_address='luke42'"
OAR_JOB_ID=3875408
oarsub -S "./run.sh 501760 32" -p "network_address='luke42'"
OAR_JOB_ID=3875409
oarsub -S "./run.sh 250880 16" -p "network_address='luke42'"
OAR_JOB_ID=3875410
oarsub -S "./run.sh 125440 8" -p "network_address='luke42'"
OAR_JOB_ID=3875411
oarsub -S "./run.sh 62720 4" -p "network_address='luke42'"
OAR_JOB_ID=3875412
oarsub -S "./run.sh 31360 2" -p "network_address='luke42'"
OAR_JOB_ID=3875413
oarsub -S "./run.sh 15680 1" -p "network_address='luke42'"
OAR_JOB_ID=3875414



# coef =0
oarsub -S "./run2.sh 64225280 4096" -p "network_address='luke44'"
OAR_JOB_ID=3875417
oarsub -S "./run2.sh 32112640 2048" -p "network_address='luke44'"
OAR_JOB_ID=3875418
oarsub -S "./run2.sh 16056320 1024" -p "network_address='luke44'"
OAR_JOB_ID=3875419
oarsub -S "./run2.sh 8028160 512" -p "network_address='luke44'"
OAR_JOB_ID=3875420
oarsub -S "./run2.sh 4014080 256" -p "network_address='luke44'"
OAR_JOB_ID=3875421
oarsub -S "./run2.sh 2007040 128" -p "network_address='luke44'"
OAR_JOB_ID=3875422
oarsub -S "./run2.sh 1003520 64" -p "network_address='luke44'"
OAR_JOB_ID=3875423
oarsub -S "./run2.sh 501760 32" -p "network_address='luke44'"
OAR_JOB_ID=3875424
oarsub -S "./run2.sh 250880 16" -p "network_address='luke44'"
OAR_JOB_ID=3875425
oarsub -S "./run2.sh 125440 8" -p "network_address='luke44'"
OAR_JOB_ID=3875426
oarsub -S "./run2.sh 62720 4" -p "network_address='luke44'"
OAR_JOB_ID=3875427
oarsub -S "./run2.sh 31360 2" -p "network_address='luke44'"
OAR_JOB_ID=3875428
oarsub -S "./run2.sh 15680 1" -p "network_address='luke44'"
OAR_JOB_ID=3875429

# Bass coeff : 1.
sh ./list3.sh
oarsub -S "./run3.sh 64225280 4096" -p "network_address='luke43'"
OAR_JOB_ID=3875484
oarsub -S "./run3.sh 32112640 2048" -p "network_address='luke43'"
OAR_JOB_ID=3875485
oarsub -S "./run3.sh 16056320 1024" -p "network_address='luke43'"
OAR_JOB_ID=3875486
oarsub -S "./run3.sh 8028160 512" -p "network_address='luke43'"
OAR_JOB_ID=3875487
oarsub -S "./run3.sh 4014080 256" -p "network_address='luke43'"
OAR_JOB_ID=3875488
oarsub -S "./run3.sh 2007040 128" -p "network_address='luke43'"
OAR_JOB_ID=3875489
oarsub -S "./run3.sh 1003520 64" -p "network_address='luke43'"
OAR_JOB_ID=3875490
oarsub -S "./run3.sh 501760 32" -p "network_address='luke43'"
OAR_JOB_ID=3875491
oarsub -S "./run3.sh 250880 16" -p "network_address='luke43'"
OAR_JOB_ID=3875492
oarsub -S "./run3.sh 125440 8" -p "network_address='luke43'"
OAR_JOB_ID=3875493
oarsub -S "./run3.sh 62720 4" -p "network_address='luke43'"
OAR_JOB_ID=3875494
oarsub -S "./run3.sh 31360 2" -p "network_address='luke43'"
OAR_JOB_ID=3875495
oarsub -S "./run3.sh 15680 1" -p "network_address='luke43'"
OAR_JOB_ID=3875496

# fretless, coeff = 0.9
oarsub -S "./run_fretless.sh 64225280 4096" -p "network_address='luke42'"
OAR_JOB_ID=3876050
oarsub -S "./run_fretless.sh 32112640 2048" -p "network_address='luke42'"
OAR_JOB_ID=3876051
oarsub -S "./run_fretless.sh 16056320 1024" -p "network_address='luke42'"
OAR_JOB_ID=3876052
oarsub -S "./run_fretless.sh 8028160 512" -p "network_address='luke42'"
OAR_JOB_ID=3876053
oarsub -S "./run_fretless.sh 4014080 256" -p "network_address='luke42'"
OAR_JOB_ID=3876054
oarsub -S "./run_fretless.sh 2007040 128" -p "network_address='luke42'"
OAR_JOB_ID=3876055
oarsub -S "./run_fretless.sh 1003520 64" -p "network_address='luke42'"
OAR_JOB_ID=3876056
oarsub -S "./run_fretless.sh 501760 32" -p "network_address='luke42'"
OAR_JOB_ID=3876057
oarsub -S "./run_fretless.sh 250880 16" -p "network_address='luke42'"
OAR_JOB_ID=3876058
oarsub -S "./run_fretless.sh 125440 8" -p "network_address='luke42'"
OAR_JOB_ID=3876059
oarsub -S "./run_fretless.sh 62720 4" -p "network_address='luke42'"
OAR_JOB_ID=3876060
oarsub -S "./run_fretless.sh 31360 2" -p "network_address='luke42'"
OAR_JOB_ID=3876061
oarsub -S "./run_fretless.sh 15680 1" -p "network_address='luke42'"
OAR_JOB_ID=3876062

	

# fretless, coeff = 0.
oarsub -S "./run_fretless_bis.sh 64225280 4096" -p "network_address='luke44'"
OAR_JOB_ID=3876063
oarsub -S "./run_fretless_bis.sh 32112640 2048" -p "network_address='luke44'"
OAR_JOB_ID=3876064
oarsub -S "./run_fretless_bis.sh 16056320 1024" -p "network_address='luke44'"
OAR_JOB_ID=3876065
oarsub -S "./run_fretless_bis.sh 8028160 512" -p "network_address='luke44'"
OAR_JOB_ID=3876066
oarsub -S "./run_fretless_bis.sh 4014080 256" -p "network_address='luke44'"
OAR_JOB_ID=3876067
oarsub -S "./run_fretless_bis.sh 2007040 128" -p "network_address='luke44'"
OAR_JOB_ID=3876068
oarsub -S "./run_fretless_bis.sh 1003520 64" -p "network_address='luke44'"
OAR_JOB_ID=3876069
oarsub -S "./run_fretless_bis.sh 501760 32" -p "network_address='luke44'"
OAR_JOB_ID=3876070
oarsub -S "./run_fretless_bis.sh 250880 16" -p "network_address='luke44'"
OAR_JOB_ID=3876071
oarsub -S "./run_fretless_bis.sh 125440 8" -p "network_address='luke44'"
OAR_JOB_ID=3876072
oarsub -S "./run_fretless_bis.sh 62720 4" -p "network_address='luke44'"
OAR_JOB_ID=3876073
oarsub -S "./run_fretless_bis.sh 31360 2" -p "network_address='luke44'"
OAR_JOB_ID=3876074
oarsub -S "./run_fretless_bis.sh 15680 1" -p "network_address='luke44'"
OAR_JOB_ID=3876075

