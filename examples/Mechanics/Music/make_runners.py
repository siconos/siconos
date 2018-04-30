"""Tools to generate 'sh' files to run simulations on Luke, using oar
with the proper parameters.
"""

import numpy as np
import os
import glob
import pickle


# Create oar file ...
def create_jobs(outputfile, runner, restit_node, job_name, freqs,
                final_time, matlab_input):
    """
    Usage :

    1) create sh file, e.g.

    create_jobs('list.sh', './run.sh', {0.: 38, 0.9:42, 1.:44}, 'bass_guitar', [1000, 2000], 1., "one_contact/pb1")
    
    2) launch oar jobs

    sh ./list.sh
    

    - list.sh : output file name
    - run.sh : oar driver
    - restit_node : dict, keys==restit coeff values, values == Luke node number.
      e.g. {0.9 : 42} means that all jobs for restitution coeff equal to 0.9 will be launched on Luke42.
    - job_name : oar job name (ie arg for oarsub --name ...)
    - freqs : list of frequencies to be considered.
    - final_time : duration
    - matlab_input : path/prefix where matlab description of parameters is to be found.

    """
    outputs = 2 ** np.arange(len(freqs))
    command = 'oarsub -S \"' + runner
    endcommand = ' -p \"network_address=\'luke'
    file = open(outputfile, 'w')
    file.write('#!/bin/bash\n')
    for e in restit_node.keys():
        file.write('# e = ' + str(e) + '\n')
        node = restit_node[e]
        for i in range(len(freqs)):
            oarline = command + ' ' + str(freqs[i]) + ' ' + str(outputs[i])
            oarline += ' ' + str(e) + ' ' + str(final_time) + ' ' + matlab_input
            oarline += '\"'
            oarline += endcommand
            oarline += str(node) + '\'\" --name=' + job_name + str(e) + '\n'  
            print(oarline)
            file.write(oarline)
        
    file.close()





### Build jobs driver for single contact case ###

final_time = 4.
freqs = [5000, 10000, 100000, 1000000, 10000000, 100000000]
matlab_input = "one_contact/pb1"

#create_jobs('all_one_contact.sh', './run_one_contact.sh', {0. : 38, 0.5: 38, 1.: 39}, 'one_contact_2018_', freqs, final_time, matlab_input)

indices = np.arange(3, 18)   
freqs_bass = 2 ** indices * 1960
final_time = 0.4
matlab_input = "one_contact/guitare_obst0"
create_jobs('all_single.sh', './run.sh', {0. : 42, 0.9: 43, 1.: 44}, 'bass_one_c_', freqs_bass, final_time, matlab_input)

### Build jobs driver for bass case ###

indices = np.arange(3, 18)   
freqs_bass = 2 ** indices * 1960
final_time = 0.4
matlab_input = "bass_guitar/pb2"
create_jobs('all_bass.sh', './run.sh', {0. : 42, 0.9: 43, 1.: 44}, 'bass_2018_', freqs_bass, final_time, matlab_input)


### Build jobs driver for fretless case ###

indices = np.arange(3, 16)   
freqs_fretless = 2 ** indices * 1960
matlab_input = "fretless_bass_guitar/bsf"

create_jobs('all_fretless.sh', './run.sh', {0. : 42, 0.9: 44, 1.: 40}, 'fretless_2018_', freqs_fretless, final_time, matlab_input)


