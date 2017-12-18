import numpy as np
import os
import glob
import pickle


# Create oar file ...
def create_jobs(outputfile, runner, restit_node, job_name, freqs,
                final_time, matlab_input):
    """
    Usage :

    create_jobs('list.sh', './run.sh', {0.: 38, 0.9:42, 1.:44}, 'bass_guitar')
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





# single contact

final_time = 4.
freqs = [5000, 10000, 100000, 1000000, 10000000, 100000000]
matlab_input = "one_contact/pb1"
create_jobs('all_one_contact.sh', './run_one_contact.sh', {0. : 42, 0.5: 42, 1.: 42}, 'one_contact', freqs, final_time, matlab_input)

# bass and fretless campaigns

indices = np.arange(3, 18)   
freqs_bass = 2 ** indices * 1960
final_time = 0.4
matlab_input = "bass_guitar/pb2"
create_jobs('all_bass.sh', './run.sh', {0. : 42, 0.9: 42, 1.: 42}, 'bass', freqs_bass, final_time, matlab_input)

indices = np.arange(3, 16)   
freqs_fretless = 2 ** indices * 1960
matlab_input = "fretless_bass_guitar/bsf"
create_jobs('all_fretless.sh', './run.sh', {0. : 42, 0.9: 42, 1.: 42}, 'fretless', freqs_fretless, final_time, matlab_input)


