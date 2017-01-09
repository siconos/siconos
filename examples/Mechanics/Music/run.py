"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
from string_ds import StringDS
from fret import Fret, Guitar

# ---- Description of the string ---
# -- Geometry and material --
length = 1.002
diameter = 0.43e-3
density = 1.17e-3
# B is a function of Young Modulus (E),
# moment of inertia ...
# B = pi^2EI / (TL^2)
B = 1.78e-5
# tension
T = 180.5
# A dictionnary with parameters required to compute quality factor
damping_parameters = {
    'nu_air': 1.8e-5,
    'rho_air': 1.2,
    'delta_ve': 4.5e-3,
    '1/qte': 2.03e-4}


damping_parameters = {
    'nu_air': 0.,
    'rho_air': 0.,
    'delta_ve': 0.,
    '1/qte': 0.}

# -- Spatial discretisation (modal proj) and initial conditions --
number_of_modes = 1
ndof = number_of_modes + 2
imax = int(ndof / 2)
# -- The dynamical system(s) --
# ie guitar strings
guitar_string = StringDS(ndof, B, density, diameter, length,
                         umax=1.8e-3, imax=imax,
                         tension=T, damping_parameters=damping_parameters)

# -- The obstacles (Interactions) --
# ie guitar fret(s)
fret = Fret(guitar_string, position=[imax, 0.])

# -- The model to gather frets and strings and simulate the dynamics --
t0 = 0.
tend = 0.01
guitar_model = Guitar(guitar_string, fret, [t0, tend], fs=8000.)

# -- Run the simulation --

k = 1
print("Start simulation ...")
while guitar_model.simu.hasNextEvent():
    guitar_model.simu.computeOneStep()
    guitar_model.save_state(k)
    k += 1
    guitar_model.simu.nextStep()

print 'End of simulation process.'

# -- plot state vs time --
guitar_model.plot_state()

# -- create string animation --
guitar_model.plot_modes('string.mp4')
