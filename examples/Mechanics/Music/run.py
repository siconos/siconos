from string_ds import StringDS
import math
import numpy as np
import siconos.kernel as sk


# initial conditions
number_of_modes = 20
ndof = number_of_modes + 2
q0 = np.zeros((ndof, ), dtype=np.float64)
q0[math.floor(ndof / 2)] = 1.8e-3
v0 = np.zeros_like(q0)
length = 1.002
diameter = 0.43e-3
density = 1.17e-3
B = 1.78e-5
T = 180.5
# A dictionnary with parameters required to compute quality factor
damping_parameters = {
    'nu_air': 1.8e-5,
    'rho_air': 1.2,
    'delta_ve': 4.5e-3,
    '1/qte': 2.03e-4
    }

# The dynamical system (geometry, material ...)
guitar = StringDS(ndof, B, density, diameter, length,
                  q0, v0, T, damping_parameters)


print guitar.K()
print guitar.C()
print guitar.q()


# Interaction
# string-floor
H = np.zeros((ndof, ), dtype=np.float64)
b = [0.]
ic = ndof / 2
H[ic] = 1
e = 0.9
nslaw = sk.NewtonImpactNSL(e)
relation = sk.LagrangianLinearTIR(H, b)
inter = sk.Interaction(nslaw, relation)

# Model
t0 = 0.
tend = 3.
guitar_model = sk.Model(t0, tend)

# add the dynamical system to the non smooth dynamical system
guitar_model.nonSmoothDynamicalSystem().insertDynamicalSystem(guitar)

# link the interaction and the dynamical system
guitar_model.nonSmoothDynamicalSystem().link(inter, guitar)

# Simulation
# (1) OneStepIntegrators
theta = 0.5
OSI = sk.MoreauJeanOSI(theta)

# (2) Time discretisation --
Fs = 1960
time_step = 1. / Fs
t = sk.TimeDiscretisation(t0, time_step)

# (3) one step non smooth problem
osnspb = sk.LCP()

# (4) Simulation setup with (1) (2) (3)
simu = sk.TimeStepping(t, OSI, osnspb)


# end of model definition

#
# computation
#

# simulation initialization
guitar_model.setSimulation(simu)
guitar_model.initialize()

while simu.hasNextEvent():
    simu.computeOneStep()
    #    k += 1
    simu.nextStep()

print guitar.q()
# print mass
# print mass.shape
# #guitar.computeMass()
