"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.


Bass guitar with frets.
Restitution coeff = 0.9

"""
import sys
import time
import os
from model_tools import save_simu_to_hdf5, build_frets_from_file
from guitar import StringDS, Guitar, Fret

# ---- Read (command line inputs) parameters ---
# - sample freq
# - output (time) frequency
# - restitution coeff
# - final time
# - path to matlab inputs
# sample freq and  time-discretisation

### Default parameters of the model (string and contacts) and of the simulation ###

## Case 1 : one contact
one_contact = {'fs': 5000, 'output_freq' : 1, 'restit': 1., 'final_time': 3.,
               'matlab_input': 'one_contact/pb1', 'nb_modes': 999, 'length': 1.,
               'output_name': 'single_e', 'filt_frets': True}

## Case 2 : bass guitar, with frets
bass_guitar = {'fs': 15680, 'output_freq' : 1, 'restit': 1., 'final_time': 1.,
               'matlab_input': 'bass_guitar/pb2', 'nb_modes': 862, 'length': .863,
               'output_name': 'bass_e', 'filt_frets': True}

## Case 3 : fretless bass guitar
one_contact = {'fs': 15680, 'output_freq' : 1, 'restit': 1.,
               'final_time': 1., 'matlab_input':'fretless_bass_guitar/bsf',
               'nb_modes': 862, 'length': .863,
               'output_name': 'fretless_e', 'filt_frets': False}


### Select case ###

ocp = one_contact

""" Read some values from command line

 If set, command line arguments must be:

 python run.py fs output_freq restitution_coeff final_time matlab_input

"""

if len(sys.argv) > 1:
    fs = float(sys.argv[1])
    output_freq = int(sys.argv[2])
    restit = float(sys.argv[3])
    final_time = float(sys.argv[4])
    matlab_input = sys.argv[5]

else:
    # No inputs from command --> use default values
    fs = ocp['fs']
    output_freq = ocp['output_freq']
    restit = ocp['restit']
    final_time = ocp['final_time']
    matlab_input = ocp['matlab_input']

   
number_of_modes = ocp['nb_modes']
initial_time = 0.
max_coords = None#(1., 0.5)
filt_frets = ocp['filt_frets']
output_name = ocp['output_name']
output_name += str(restit)
visu = True # plot frets/neck 
# -- Geometry and material --
# Indeed, since parameters are read from matlab input, only length matter.
G_string = {
    'length': ocp['length']
    # diameter = equivalent diameter (A5)
    #'diameter': 1.14e-3,
    #'density': 6.69e-3,
    #'B': 3.5e-5,
    #'tension': 191.6,
    }


# -- The dynamical system --
string = StringDS(number_of_modes, geometry_and_material=G_string,
                  max_coords=max_coords,
                  matlab_input=matlab_input)


current_path = os.path.dirname(os.path.realpath(__file__))

# -- The interactions --
frets_file = os.path.join(current_path, matlab_input) + '_h.mat'
interactions = build_frets_from_file(string, restit, frets_file, filt_frets, visu)

frets = list(interactions.keys())
nb_frets = len(frets)

# -- The nsds --
guitar_model = Guitar(interactions, [initial_time, final_time],
               fs, output_freq,
               interactions_output=2) # 2 to save y and lambda

if not guitar_model.save_interactions:
    print("Warning! No interactions output!")
# Save initial state
guitar_model.time[0] = initial_time
guitar_model.save_ds_state_modal(0, string)

if guitar_model.save_interactions:
    buff = guitar_model.data_interactions
    for i in range(nb_frets):
        buff[frets[i]][0][0] = frets[i].y(0) 
        buff[frets[i]][1][0] = frets[i].lambda_(1) 

print('Ready to start simulation for frequency {0}.'.format(fs))
print('Save output every {0} time steps.'.format(output_freq))
msg = 'Read data from files :\n'
msg += '- neck profile:' + frets_file
msg += '\n- eigenfrequencies: ' + matlab_input + '_frequs.mat\n'
msg += '- damping: ' + matlab_input + '_amortissements.mat\n'
print(msg)

# Get simulation 
simu = guitar_model.simulation
# -- Model setup with dynamics, interaction and simulation --

k = 1
print("Start simulation ...")
start_time = time.clock()
pos = 1

while simu.hasNextEvent():
    if k % 100000 == 0:
        print('step = ', k, '---- time = ',
              simu.nextTime(),
              '------------------------')
    simu.computeOneStep()

    # -- save data every output_freq time step --
    if k % guitar_model.output_freq == 0:
        # current time
        guitar_model.time[pos] = simu.nextTime()
        # modal positions
        guitar_model.save_ds_state_modal(pos, string)
        
        # interactions
        if guitar_model.save_interactions:
            buff = guitar_model.data_interactions
            for i in range(nb_frets):
                buff[frets[i]][0][pos] = frets[i].y(0) 
                buff[frets[i]][1][pos] = frets[i].lambda_(1) 
                #buff[frets[i]][2][pos = frets[i].y(1) 

        pos += 1
    k += 1
    simu.nextStep()
print('End of simulation process. Duration: ', time.clock() - start_time)
print("nb steps", k)
# --- Output dir for results ---
result_dir = os.getcwd()# + '/temp'
if not os.path.exists(result_dir):
    os.mkdir(result_dir)
filename = output_name + '_' + str(number_of_modes) + '_' + str(int(fs)) + '.h5'
filename = os.path.join(result_dir, filename)
start = time.clock()
save_simu_to_hdf5(guitar_model, string,
                   matlab_data=matlab_input,
                   filename=filename, filt_frets=filt_frets,
                   restit=restit)
print('write (hdf5) file ' + filename + ': {0}.'.format(time.clock() - start))
