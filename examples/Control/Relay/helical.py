"""Helical example,
from B. Caillaud,
http://www.irisa.fr/prive/Benoit.Caillaud/cours-hybride-2014/Cours_modelisation_des_systemes_hybrides/Cours/Cours.html
"""

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import siconos.kernel as sk
import siconos.numerics as sn
from math import ceil

# == User-defined parameters ==
ndof = 3  # number of degrees of freedom of your system
t0 = 0.0
T = 400  # Total simulation times
h = 0.05  # Time step
rho = 10.0
x00 = rho  # Initial position
x01 = rho
x02 = 0.0
alpha = 0.05  # angle of the square helical
beta = 0.01  # thread of the helical
gamma = 0.  # thread variation

# -- Dynamical system --
# dx / dt = A.x + b + r
A = np.zeros((ndof, ndof), dtype=np.float64)
x0 = np.zeros(3, dtype=np.float64)
x0.flat[...] = [x00, x01, x02]
b = np.zeros_like(x0)
b[2] = beta
particle = sk.FirstOrderLinearDS(x0, A, b)

# -- Interaction --
# y = C.x + D.lambda
# r = B.lambda
ninter = 2
B = np.zeros((ndof, ninter), dtype=np.float64)
B[0, 0] = -alpha * 0.5
B[0, 1] = 1 + alpha * 0.5
B[1, 0] = -B[0, 1]
B[1, 1] = B[0, 0]
B[2, 0] = gamma * 0.5
B[2, 1] = gamma * 0.5

C = np.zeros((ninter, ndof), dtype=np.float64)
C[0, 0] = C[1, 1] = -1.

particle_relation = sk.FirstOrderLinearR(C, B)

nslaw = sk.RelayNSL(ninter)

particle_interaction = sk.Interaction(nslaw, particle_relation)

# -- The Model --
filippov = sk.Model(t0, T)
nsds = filippov.nonSmoothDynamicalSystem()
nsds.insertDynamicalSystem(particle)
nsds.link(particle_interaction, particle)

# -- Simulation --
td = sk.TimeDiscretisation(t0, h)
simu = sk.TimeStepping(td)
# osi
theta = 0.5
myIntegrator = sk.EulerMoreauOSI(theta)
simu.insertIntegrator(myIntegrator)

# osns
osnspb = sk.Relay(sn.SICONOS_RELAY_LEMKE)
simu.insertNonSmoothProblem(osnspb)

filippov.setSimulation(simu)
filippov.initialize()

# -- Get the values to be plotted --
output_size = 1 + ndof + 2 * ninter
nb_time_steps = int((T - t0) / h) + 1
data_plot = np.empty((nb_time_steps, output_size))
data_plot[0, 0] = filippov.t0()
data_plot[0, 1:4] = particle.x()
data_plot[0, 4:] = 0.

# time loop
k = 1
while simu.hasNextEvent():
    simu.computeOneStep()
    data_plot[k, 0] = simu.nextTime()
    data_plot[k, 1:4] = particle.x()
    data_plot[k, 4:6] = particle_interaction.lambda_(0)[0:2]
    data_plot[k, 6:8] = particle_interaction.y(0)[0:2]
    k += 1
    simu.nextStep()


# save to disk
np.savetxt('helical.dat', data_plot)


def plot_results():
    """Plot 3d curve z = f(x, y)
    with ds_state = [x, y, z]
    """
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x = data_plot[:, 1]
    y = data_plot[:, 2]
    z = data_plot[:, 3]
    ax.plot(x, y, z, label='z = f(x, y)')
    ax.legend()
    plt.savefig('helical.png')


# --- Uncomment lines below to plot interesting stuff ---
plot_results()
# subplot(311)
# title('position')
# plot(data_plot[:,0], data_plot[:,1])
# grid()
# subplot(312)
# title('velocity')
# plot(data_plot[:,0], data_plot[:,2])
# grid()
# subplot(313)
# plot(data_plot[:,0], data_plot[:,3])
# title('lambda')
# grid()
# show()

# plot(data_plot[:,1], data_plot[:,2])
# grid()
# show()

