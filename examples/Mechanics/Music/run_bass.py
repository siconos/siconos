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

from guitar import StringDS, Guitar

# ---- Read (command line inputs) parameters ---
# - sample freq
# - output (time) frequency
# - restitution coeff
# - final time
# - path to matlab inputs
# sample freq and  time-discretisation

current_path = os.path.dirname(os.path.realpath(__file__))
# if freq is set as input arg ...
if len(sys.argv) > 1:
    fs = float(sys.argv[1])
    output_freq = int(sys.argv[2])
    restit = float(sys.argv[3])
    final_time = float(sys.argv[4])
    matlab_input = sys.argv[5]

else:
    fs = 15680.
    output_freq = 1
    restit = 1.
    final_time = 1.
    matlab_input = 'bass_guitar/pb2'

number_of_modes = 862
    
assert fs > 14000.
# lower frequencies require quadruple prec for exp computation.


initial_time = 0.

max_coords = (3.6e-3, 0.64)

#max_coords = (7.8e-3, 0.64)

if matlab_input.find('fretless') >= 0:
    filt_frets = False # To select only 'real' frets
    output_name = 'fretless_e'
else:
    filt_frets = True
    output_name = 'bass_e'
output_name += str(restit)

visu = False # plot frets/neck 

# -- Geometry and material --
# Indeed, since parameters are read from matlab input, only length matter.
G_string = {
    'length': 0.863,
    # diameter = equivalent diameter (A5)
    #'diameter': 1.14e-3,
    #'density': 6.69e-3,
    #'B': 3.5e-5,
    #'tension': 191.6,
    }


# -- The dynamical system --
string = StringDS(number_of_modes, geometry_and_material=G_string,
                  max_coords=max_coords,
                  matlab_input=matlab_input,
                  )

# -- The interactions --
frets_file = os.path.join(current_path, matlab_input) + '_h.mat'
interactions = build_frets_from_file(string, restit, frets_file, filt_frets, visu)

frets = list(interactions.keys())
nb_frets = len(frets)

# -- The nsds --
guitar_model = Guitar(interactions, [initial_time, final_time],
               fs, output_freq,
                      interactions_output=2) # 2 to save y and lambda

# Save initial state
guitar_model.time[0] = initial_time
guitar_model.save_ds_state_modal(0, string)

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
simu = guitar_model.simulation()
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
