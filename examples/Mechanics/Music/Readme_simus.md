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

oarsub -S "./run_fretless.sh 32112640 2048" -p "network_address='luke42'"

oarsub -S "./run_fretless.sh 16056320 1024" -p "network_address='luke42'"

oarsub -S "./run_fretless.sh 8028160 512" -p "network_address='luke42'"

oarsub -S "./run_fretless.sh 4014080 256" -p "network_address='luke42'"

oarsub -S "./run_fretless.sh 2007040 128" -p "network_address='luke42'"

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

oarsub -S "./run_fretless_bis.sh 32112640 2048" -p "network_address='luke44'"

oarsub -S "./run_fretless_bis.sh 16056320 1024" -p "network_address='luke44'"

oarsub -S "./run_fretless_bis.sh 8028160 512" -p "network_address='luke44'"

oarsub -S "./run_fretless_bis.sh 4014080 256" -p "network_address='luke44'"

oarsub -S "./run_fretless_bis.sh 2007040 128" -p "network_address='luke44'"

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


oarsub -S "./run_fretless.sh 64225280 4096" -p "network_address='luke42'"
OAR_JOB_ID=3877878
oarsub -S "./run_fretless.sh 32112640 2048" -p "network_address='luke42'"
OAR_JOB_ID=3877879
oarsub -S "./run_fretless.sh 16056320 1024" -p "network_address='luke42'"
OAR_JOB_ID=3877880
oarsub -S "./run_fretless.sh 8028160 512" -p "network_address='luke42'"
OAR_JOB_ID=3877881
oarsub -S "./run_fretless.sh 4014080 256" -p "network_address='luke42'"
OAR_JOB_ID=3877882
oarsub -S "./run_fretless.sh 2007040 128" -p "network_address='luke42'"
OAR_JOB_ID=3877883
oarsub -S "./run_fretless.sh 1003520 64" -p "network_address='luke42'"
OAR_JOB_ID=3877884
oarsub -S "./run_fretless.sh 501760 32" -p "network_address='luke42'"
OAR_JOB_ID=3877885
oarsub -S "./run_fretless.sh 250880 16" -p "network_address='luke42'"
OAR_JOB_ID=3877886
oarsub -S "./run_fretless.sh 125440 8" -p "network_address='luke42'"
OAR_JOB_ID=3877887
oarsub -S "./run_fretless.sh 62720 4" -p "network_address='luke42'"
OAR_JOB_ID=3877888
oarsub -S "./run_fretless.sh 31360 2" -p "network_address='luke42'"
OAR_JOB_ID=3877889
oarsub -S "./run_fretless.sh 15680 1" -p "network_address='luke42'"
OAR_JOB_ID=3877890




oarsub -S "./run_fretless_bis.sh 64225280 4096" -p "network_address='luke43'"
OAR_JOB_ID=3877891
oarsub -S "./run_fretless_bis.sh 32112640 2048" -p "network_address='luke43'"
OAR_JOB_ID=3877892
oarsub -S "./run_fretless_bis.sh 16056320 1024" -p "network_address='luke43'"
OAR_JOB_ID=3877893
oarsub -S "./run_fretless_bis.sh 8028160 512" -p "network_address='luke43'"
OAR_JOB_ID=3877894
oarsub -S "./run_fretless_bis.sh 4014080 256" -p "network_address='luke43'"
OAR_JOB_ID=3877895
oarsub -S "./run_fretless_bis.sh 2007040 128" -p "network_address='luke43'"
OAR_JOB_ID=3877896
oarsub -S "./run_fretless_bis.sh 1003520 64" -p "network_address='luke43'"
OAR_JOB_ID=3877897
oarsub -S "./run_fretless_bis.sh 501760 32" -p "network_address='luke43'"
OAR_JOB_ID=3877898
oarsub -S "./run_fretless_bis.sh 250880 16" -p "network_address='luke43'"
OAR_JOB_ID=3877899
oarsub -S "./run_fretless_bis.sh 125430 8" -p "network_address='luke43'"
OAR_JOB_ID=3877900
oarsub -S "./run_fretless_bis.sh 62720 4" -p "network_address='luke43'"
OAR_JOB_ID=3877901
oarsub -S "./run_fretless_bis.sh 31360 2" -p "network_address='luke43'"
OAR_JOB_ID=3877902
oarsub -S "./run_fretless_bis.sh 15680 1" -p "network_address='luke43'"
OAR_JOB_ID=3877903

oarsub -S "./run.sh 63825280 4096" -p "network_address='luke38'"
OAR_JOB_ID=3877919
oarsub -S "./run.sh  32112640 2048" -p "network_address='luke38'"
OAR_JOB_ID=3877920
oarsub -S "./run.sh 16056320 1024" -p "network_address='luke38'"
OAR_JOB_ID=3877921
oarsub -S "./run.sh 8028160 512" -p "network_address='luke38'" 
OAR_JOB_ID=3877922
oarsub -S "./run.sh 4014080 256" -p "network_address='luke38'"
OAR_JOB_ID=3877923
oarsub -S "./run.sh 2007040 128" -p "network_address='luke38'"
OAR_JOB_ID=3877924
oarsub -S "./run.sh 1003520 64" -p "network_address='luke38'"
OAR_JOB_ID=3877925
oarsub -S "./run.sh 501760 32" -p "network_address='luke38'"
OAR_JOB_ID=3877926
oarsub -S "./run.sh 250880 16" -p "network_address='luke38'"
OAR_JOB_ID=3877927
oarsub -S "./run.sh 125440 8" -p "network_address='luke38'"
OAR_JOB_ID=3877928
oarsub -S "./run.sh 62720 4" -p "network_address='luke38'"
OAR_JOB_ID=3877929
oarsub -S "./run.sh 31360 2" -p "network_address='luke38'"
OAR_JOB_ID=3877930
oarsub -S "./run.sh 15680 1" -p "network_address='luke38'"
OAR_JOB_ID=3877931


oarsub -S "./run2.sh 64225280 4096" -p "network_address='luke44'"
OAR_JOB_ID=3877938
oarsub -S "./run2.sh 32112640 2048" -p "network_address='luke44'"
OAR_JOB_ID=3877939
oarsub -S "./run2.sh 16056320 1024" -p "network_address='luke44'"
OAR_JOB_ID=3877940
oarsub -S "./run2.sh 8028160 512" -p "network_address='luke44'"
OAR_JOB_ID=3877941
oarsub -S "./run2.sh 4014080 256" -p "network_address='luke44'"
OAR_JOB_ID=3877942
oarsub -S "./run2.sh 2007040 128" -p "network_address='luke44'"
OAR_JOB_ID=3877943
oarsub -S "./run2.sh 1003520 64" -p "network_address='luke44'"
OAR_JOB_ID=3877944
oarsub -S "./run2.sh 501760 32" -p "network_address='luke44'"
OAR_JOB_ID=3877945
oarsub -S "./run2.sh 250880 16" -p "network_address='luke44'"
OAR_JOB_ID=3877946
oarsub -S "./run2.sh 125440 8" -p "network_address='luke44'"
OAR_JOB_ID=3877947
oarsub -S "./run2.sh 62720 4" -p "network_address='luke44'"
OAR_JOB_ID=3877948
oarsub -S "./run2.sh 31360 2" -p "network_address='luke44'"
OAR_JOB_ID=3877949
oarsub -S "./run2.sh 15680 1" -p "network_address='luke44'"
OAR_JOB_ID=3877950

