"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
from guitar import StringDS, Fret, Guitar
#import siconos.kernel as sk
import time

import numpy as np
import scipy.io
import sys,os

visu=False
if visu:
    import matplotlib.pyplot as plt    

# ======== Description of the string(s) ==========
# -- Geometry and material --
G_string = {
    'length': 0.863,
    'diameter': 0.43e-3,
    'density': 6.69e-3,
    'B': 3.5e-5,
    'tension': 191.6,
}

# A dictionnary with parameters required to compute quality factor
damping_parameters = {
    'nu_air': 1.8e-5,
    'rho_air': 1.2,
    'delta_ve': 0.01,
    '1/qte': 6e-6}

# -- Spatial discretisation (modal proj) and initial conditions --
number_of_modes = 864
# position (index) of the max value of initial state
# --> middle of the string
#imiddle = int((number_of_modes + 2) / 2)

# -- The dynamical system(s) --

# Warning: 'real dofs' numbers start from 0 to number_of_modes + 1
# but DS size is number_of_modes, boundary points are ignored.
guitar_string = StringDS(number_of_modes, geometry_and_material=G_string,
                         damping_parameters=damping_parameters,
                         max_coords=(7.8e-3, .64))

# -- The interaction(s) between strings and frets --
# One Interaction is needed for each contact point.

frets_file = './donnees_siconos/pb2_h.mat'
frets_positions = scipy.io.loadmat(frets_file)['h'][:600, 0]
nb_frets = frets_positions.size
dx = guitar_string.space_step
frets_indices = np.arange(1, nb_frets + 1)
#array(np.round(frets_positions / dx), np.int32)
#frets_y = np.linspace(-0.5e-4, -1e-4, 20)

frets = []
interactions = {}
for i in range(nb_frets):
    frets.append(Fret(guitar_string,
                      contact_positions=(frets_indices[i],
                                         frets_positions[i]),
                      restitution_coeff=0.))
    interactions[frets[-1]] = guitar_string
# contact at a point close to left boundary.
# guitar_fret_left = Fret(guitar_string, contact_positions=(3, -1.e-4),
#                         restitution_coeff=0.9)

# -- Model and simulaton --


# sample freq and time-discretisation

# if freq is set as input arg ...
if len(sys.argv) > 1:
    fe = float(sys.argv[1])
else:
    fe = 1960
initial_time = 0.
final_time = 0.3
final_time = 3.00

guitar_model = Guitar(interactions,
                      #{guitar_fret_middle: guitar_string,
                      # guitar_fret_left: guitar_string
                      # },
                      [initial_time, final_time],
                      fe)


simu = guitar_model.simulation()
# sk.TimeStepping(time_discr, osi, osnspb)

# # -- Model setup with dynamics, interaction and simulation --

# -- Save inital state --
# Note about savings:
#  For each iteration k, we need :
#  - the current time value, saved in guitar_model.time[k]
#  - data for ds (positions, velocities ...):
#     use save_ds_state(k, ds) for each required ds
#  - data for interactions (impulsion at impact, distance ...)
#     use save_interaction_state(k, interaction)
guitar_model.time[0] = initial_time
guitar_model.save_ds_state(0, guitar_string)
for i in range(len(frets)):
    guitar_model.save_interaction_state(0, frets[i])

k = 1
print("Start simulation ...")
start_time = time.clock()
while simu.hasNextEvent():
    if k % 100 == 0:
        print('step = ', k, '---- time = ',
              simu.nextTime(),
              '------------------------')
    simu.computeOneStep()
    guitar_model.time[k] = simu.nextTime()
    guitar_model.save_ds_state(k, guitar_string)
    for i in range(len(frets)):
        guitar_model.save_interaction_state(k, frets[i])
    k += 1
    simu.nextStep()
print('End of simulation process. Duration: ', time.clock() - start_time)


# -- Save results for ds in numpy file +--
result_dir = 'results'
if not os.path.exists(result_dir):
        os.mkdir(result_dir)

output = guitar_model.data_ds[guitar_string]
filename = os.path.join(result_dir,'data_ds_'+str(number_of_modes)+'_'+str(fe))
np.save(filename, output)

# -- Save results for interaction in numpy file --
for i in range(len(frets)):
    output = guitar_model.data_interactions[frets[i]]
    filename = os.path.join(result_dir,'data_interactions_'+str(i)+'_'+str(number_of_modes)+'_'+str(fe))
    np.save(filename, output)
# to plot results, call:
#guitar_model.plot_ds_state(some_ds, indices, fig_number)
# --> plot some_ds attributes (position/time ...)
#guitar_model.plot_interaction(some_interaction, fig_number)
# --> plot data relative to some_interaction
# guitar_model.plot_mode(some_ds, filename)
# --> create animation for some_ds mode


#guitar_model.plot_ds_state(guitar_string)
