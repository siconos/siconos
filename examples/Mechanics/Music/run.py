"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
from guitar import StringDS, Fret
import siconos.kernel as sk
import matplotlib.pyplot as plt
import time
import numpy as np
import cmath
import math


# ---- Description of the string ---
# -- Geometry and material --
G_string = {
    'length': 1.002,
    'diameter': 0.43e-3,
    'density': 1.17e-3,
    'B': 1.97e-5,
    'tension': 180.5,
}

# A dictionnary with parameters required to compute quality factor

damping_parameters = {
    'nu_air': 1.8e-5,
    'rho_air': 1.2,
    'delta_ve': 4.5e-3,
    '1/qte': 2.03e-4}

# damping_parameters = {
#     'nu_air': 0.,
#     'rho_air': 0.,
#     'delta_ve': 0.,
#     '1/qte': 0.}

# -- Spatial discretisation (modal proj) and initial conditions --
number_of_modes = 10
ndof = number_of_modes #+ 2
imax = int((number_of_modes + 2) / 2)
# -- The dynamical system(s) --
guitar_string = StringDS(ndof, geometry_and_material=G_string,
                         damping_parameters=damping_parameters,
                         umax=1.44e-3, imax=imax)

# -- The interaction --
# contact in the middle of the string
guitar_fret = Fret(guitar_string, position=[imax - 1, 0.],
                   restitution_coeff=0.9)

# -- Model and simulaton --


# # (1) OneStepIntegrators
osi = sk.MoreauJeanBilbaoOSI()
# (2) Time discretisation --
# sample freq and time-discretisation
fe = 10 * 51200
initial_time = 0.
time_step = 1. / fe
final_time = 0.05
time_discr = sk.TimeDiscretisation(initial_time, time_step)
nb_time_steps = (int)((final_time - initial_time) / time_step)
# (3) one step non smooth problem
osnspb = sk.LCP()
# (4) Simulation setup with (1) (2) (3)
simu = sk.TimeStepping(time_discr, osi, osnspb)

# -- Model setup with dynamics, interaction and simulation --

guitar_model = sk.Model(initial_time, final_time)
# insert ds into model ...
guitar_model.nonSmoothDynamicalSystem().insertDynamicalSystem(guitar_string)
# link the interaction(s) and the dynamical system(s)
guitar_model.nonSmoothDynamicalSystem().link(guitar_fret, guitar_string)
# simulation initialisation
guitar_model.setSimulation(simu)
guitar_model.initialize()


k = 1
print("Start simulation ...")
while simu.hasNextEvent():
    if k % 100 == 0:
        print('step = ', k, '---- time = ',
              simu.nextTime(),
              '------------------------')
    simu.computeOneStep()
    #model.save_state(k)
    k += 1
    simu.nextStep()
print('End of simulation process.')

