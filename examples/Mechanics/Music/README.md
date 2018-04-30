# Simulation of the dynamics of a guitar string with siconos


3 differents cases are available :

- one string/one contact

- 'real' bass guitar : one string, 20 frets (and so 20 contacts)

- fretless guitar : one string, no frets, many contacts between string and neck.



## Usage

* with default setup:

siconos run.py 

* with input args

siconos run.py frequency output_frequency restit final_time matlab_input

- frequency : sample freq
- output_frequency : freq. on output files writing
- restit : restitution coeff
- final_time : duration
- matlab_input : path/prefix where path is the location of matlab files describing the model. e.g. : bass_guitar/pb2


## Usage on Luke cluster
On Luke, nix and oar config are required.

To submit a job, use run.sh file:

Example :

oarsub -S "./run.sh 15680 1 0.0 0.4 fretless_bass_guitar/bsf" -p "network_address='luke42'" --name=fretless_2018_0.0

parameters following run.sh are the same as those given above, when using siconos run.py.


To prepare a campaign (i.e. a set of jobs for differents frequencies and restitution coeff) use 'make_runners.py' file
to create a driver.sh. See make_runners.py for details.

Remark : each time run.sh is called, a line is appended to 'jobs_params' file, with all the values of the job parameters.
This jobs_params file will be used for postprocessing.





## Files

guitar.py : description of the dynamical system (string), the frets (interactions)
