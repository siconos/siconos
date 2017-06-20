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
    'length': 1. , #1.002,
    'diameter': 0.43e-3,
    'density': 1., #1.17e-3,
    'B': 0., #1.97e-5,
    'tension': 1.#180.5,
}

# A dictionnary with parameters required to compute quality factor
damping_parameters = {
    'nu_air': 1.8e-5,
    'rho_air': 1.2,
    'delta_ve': 4.5e-3,
    '1/qte': 2.03e-4}

# -- Spatial discretisation (modal proj) and initial conditions --
number_of_modes = 3
# position (index) of the max value of initial state
# --> middle of the string
imiddle = int((number_of_modes + 2) / 2)

# -- The dynamical system(s) --

# Warning: 'real dofs' numbers start from 0 to number_of_modes + 1
# but DS size is number_of_modes, boundary points are ignored.
guitar_string = StringDS(number_of_modes, geometry_and_material=G_string,
                         damping_parameters=damping_parameters,
                         umax=1.44e-3, imax=imiddle)

# -- The interaction(s) between strings and frets --
# One Interaction is needed for each contact point.

# contact in the middle of the string (defined by imiddle)

guitar_fret_middle = Fret(guitar_string, contact_positions=(imiddle, 0.),
                          restitution_coeff=1.)
# contact at a point close to left boundary.
guitar_fret_left = Fret(guitar_string, contact_positions=(3, 1.e-4),
                        restitution_coeff=0.9)

# -- Model and simulaton --


# sample freq and time-discretisation
fe = 10 * 51200
initial_time = 0.
final_time = 0.05

guitar_model = Guitar({guitar_fret_middle: guitar_string,
                       #guitar_fret_left: guitar_string
                       },
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
guitar_model.save_interaction_state(0, guitar_fret_middle)

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
    guitar_model.save_interaction_state(k, guitar_fret_middle)
    k += 1
    simu.nextStep()
print('End of simulation process. Duration: ', time.clock() - start_time)

# -- Save results for ds in numpy file --
#output = guitar_model.data_ds[guitar_string]
#np.save('data1001', output)


# to plot results, call:
# guitar_model.plot_ds_state(some_ds, indices, fig_number)
# --> plot some_ds attributes (position/time ...)
# guitar_model.plot_interaction(some_interaction, fig_number)
# --> plot data relative to some_interaction
# guitar_model.plot_mode(some_ds, filename)
# --> create animation for some_ds mode



