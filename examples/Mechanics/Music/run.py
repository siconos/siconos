"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
from guitar import StringDS, Fret, Guitar
import siconos.kernel as sk
import time
import matplotlib.pyplot as plt
import numpy as np


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

# contact in the middle of the string (defined by imiddle)


nb_frets = 4
#frets_pos = np.linspace(0.05, 0.6, nb_frets)
frets_pos = np.asarray([0.2, 0.3, 0.5, 0.6])
dx = guitar_string.space_step
frets_ind = np.asarray(np.round(frets_pos / dx), np.int32)
frets_y = np.asarray([0., -1.e-5, -1.3e-5, -2.e-5])
nb_frets = 20
frets_pos = np.linspace(0.05, 0.6, 20)
dx = guitar_string.space_step
frets_ind = np.asarray(np.round(frets_pos / dx), np.int32)
frets_y = np.linspace(-0.5e-4, -1e-4, 20)
frets = []
interactions = {}
for i in range(nb_frets):
    frets.append(Fret(guitar_string,
                      contact_positions=(frets_ind[i], frets_y[i]),
                      restitution_coeff=0.))
    interactions[frets[-1]] = guitar_string
# contact at a point close to left boundary.
# guitar_fret_left = Fret(guitar_string, contact_positions=(3, -1.e-4),
#                         restitution_coeff=0.9)

# -- Model and simulaton --


# sample freq and time-discretisation
fe = 51200
initial_time = 0.
final_time = 0.2

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

i0 = frets[0]
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
    count = 0
    for i in range(len(frets)):
        guitar_model.save_interaction_state(k, frets[i])
    #   guitar_model.save_interaction_state(k, guitar_fret_middle)
        if len(np.where(frets[i].lambda_(1)[0] > 1e-5)[0])>0:
            count += 1
    print("time %s, %s", str(guitar_model.time[k]), str(count))
    k += 1
    simu.nextStep()
print('End of simulation process. Duration: ', time.clock() - start_time)

# -- Save results for ds in numpy file --
#output = guitar_model.data_ds[guitar_string]
#np.save('data1001', output)


# to plot results, call:
#guitar_model.plot_ds_state(some_ds, indices, fig_number)
# --> plot some_ds attributes (position/time ...)
#guitar_model.plot_interaction(some_interaction, fig_number)
# --> plot data relative to some_interaction
# guitar_model.plot_mode(some_ds, filename)
# --> create animation for some_ds mode


#guitar_model.plot_ds_state(guitar_string)
