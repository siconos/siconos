/!\ Draft /!\


Vibrating string with impact based on 'normal modes' spatial discretisation.


sh ./all_bass.sh
sh ./all_fretless.sh
oarsub -S "./run.sh 4014080 1 0.0 0.4 bass_guitar/pb2" -p "network_address='luke42'" --name=bass_full_2018_0.0
oarsub -S "./run.sh 4014080 1 0.9 0.4 bass_guitar/pb2" -p "network_address='luke42'" --name=bass_full_2018_0.9
oarsub -S "./run.sh 4014080 1  0.0 0.4 fretless_bass_guitar/bsf" -p "network_address='luke42'" --name=fretless_full_2018_0.0
oarsub -S "./run.sh 4014080 1  0.9 0.4 fretless_bass_guitar/bsf" -p "network_address='luke42'" --name=fretless_full_2018_0.9
oarsub -S "./run.sh 4014080 1 1.0 0.4 fretless_bass_guitar/bsf" -p "network_address='luke42'" --name=fretless_full_2018_1.0
oarsub -S "./run.sh 4014080 1 1.0 0.4 bass_guitar/pb2" -p "network_address='luke42'" --name=bass_full_2018_1.0
sh ./all_one_contact.sh



Names of interest:

fretless_full_2018 : 3 simus, 4MHz, 0., 0.9, 1.0, output every time step (!! 31GB files)
bass_full_2018 : 3 simus, 4MHz, 0., 0.9, 1.0, output every time step (!! 11GB files)
bass_2018 : 45 files
fretless_2018 : 39 files
one_contact_2018 : 30 files, 0., 0.5, 0.9, 1.
