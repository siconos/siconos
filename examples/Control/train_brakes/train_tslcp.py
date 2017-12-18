"""Train brake example
from B. Caillaud,  <benoit.caillaud@inria.fr>, Inria, Rennes, September 2013

Time-Stepping with LCP.
Equations are ::

    for all i=1...n
    j_{i-1} - k_i - l_i - j_i = 0
    w_i = u_i + v_i
    w_{i-1} - R j_{i-1} - w_i = 0
    h_i + S k_i = v_i
    u'_i = 1/C (k_i + l_i)
    v'_i = 1/D l_i
    0 <= k_i complementarity -h_i >= 0

    w_0 is a constant
    j_n is a constant

"""

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import siconos.kernel as sk
from math import ceil

# == User-defined parameters ==
number_of_cars = 3
t0 = 0.0
T = 10  # Total simulation times
tau = 1e-3  # Time step
tscale = 1.

# Resistors (Ohm)
R = 100.
S = 10.

# Capacitors (Farad)
C = 1000.e-6
D = 100.e-6

# Initial state
u0 = 7.5  # C capacitors voltates (V)
v0 = - u0  # D capacitors voltages (V)
rho = 0.6667  # factor
epsilon = rho * u0 / number_of_cars  # increment of capacitor voltages (V)

# Constant perturbation
w0 = 10.0  # Voltage at the head of the brake line
jn = 0.0  # current at the tail of the brake line


# -- Dynamical system --
# *** Linear time-invariant dynamical system (LTIDS) ***
# q' = A.q + b + r with :
#
# q = [ u1 ... u_n v_1 ... v_n ]
# q(t0) = q0 = [ u0 ... u0 v0 ... v0 ]
ndof = 2 * number_of_cars
q0 = np.zeros(ndof, dtype=np.float64)
q0[:number_of_cars] = (1. - rho) * u0 + epsilon * np.arange(number_of_cars)
q0[number_of_cars:] = (1. - rho) * v0 -\
    epsilon * np.arange(number_of_cars, ndof)

A = np.zeros((ndof, ndof), dtype=np.float64)
val = -2. / (R * C)
np.fill_diagonal(A[:number_of_cars, :number_of_cars], val)
np.fill_diagonal(A[:number_of_cars, number_of_cars:ndof], val)
val2 = -2. / (R * D)
np.fill_diagonal(A[number_of_cars:ndof, :number_of_cars], val2)
np.fill_diagonal(A[number_of_cars:ndof, number_of_cars:ndof], val2)
A[number_of_cars - 1, number_of_cars - 1] *= 0.5
A[number_of_cars - 1, ndof - 1] *= 0.5
A[ndof - 1, number_of_cars - 1] *= 0.5
A[ndof - 1, ndof - 1] *= 0.5
# extra-diag values
for i in xrange(1, number_of_cars):
    A[i, i - 1] = A[i - 1, i] = 1. / (R * C)
    A[i, i + number_of_cars - 1] = A[i - 1, i + number_of_cars] = 1. / (R * C)
    A[i + number_of_cars, i - 1] = A[i + number_of_cars - 1, i] = 1. / (R * D)
    A[i + number_of_cars, i + number_of_cars - 1] = 1. / (R * D)
    A[i + number_of_cars - 1, i + number_of_cars] = 1. / (R * D)

b = np.zeros(ndof, np.float64)
b[0] = w0 / (R * C)
b[number_of_cars - 1] = - jn / C
b[number_of_cars] = w0 / (R * D)
b[ndof - 1] = -jn / D

RC = sk.FirstOrderLinearTIDS(q0, A, b)

# -- Interaction --
# *** Linear time invariant relation (LTIR) ***
# y = N.q + M.x and r = B.x
# y = [ -h_1 ... -h_n ]
# x = [ k_1 ... k_n ]
# + complementarity(x,y)

ninter = number_of_cars
B = np.zeros((ndof, ninter), dtype=np.float64)
np.fill_diagonal(B[number_of_cars:, :number_of_cars], -1. / D)

M = np.zeros((ninter, ninter), dtype=np.float64)
np.fill_diagonal(M[:number_of_cars, :number_of_cars], S)

N = np.zeros((ninter, ndof), dtype=np.float64)
np.fill_diagonal(N[:number_of_cars, number_of_cars:], -1.)

relation = sk.FirstOrderLinearTIR(N, B)
relation.setDPtr(M)

nslaw = sk.ComplementarityConditionNSL(ninter)

interaction = sk.Interaction(nslaw, relation)

# -- The Model --
circuit = sk.Model(t0, T, 'train')
nsds = circuit.nonSmoothDynamicalSystem()
nsds.insertDynamicalSystem(RC)
nsds.link(interaction, RC)

# -- Simulation --
td = sk.TimeDiscretisation(t0, tau)
simu = sk.TimeStepping(td)
# osi
theta = 0.50000000000001
osi = sk.EulerMoreauOSI(theta)
simu.insertIntegrator(osi)

# osns
osnspb = sk.LCP()
simu.insertNonSmoothProblem(osnspb)

circuit.setSimulation(simu)
circuit.initialize()

# -- Get the values to be plotted --
output_size = 2 * number_of_cars + 1
nb_time_steps = int((T - t0) / tau) + 1
data_plot = np.empty((nb_time_steps, output_size))
data_plot[0, 0] = circuit.t0() / tscale
data_plot[0, 1:] = RC.x()

# time loop
k = 1
while simu.hasNextEvent():
    # update input
    current = 2 * k - nb_time_steps
    if current >= 0 and current <= 1:
        RC.b()[0] = 0.0
        RC.b()[number_of_cars] = 0

    simu.computeOneStep()
    data_plot[k, 0] = simu.nextTime() / tscale
    data_plot[k, 1:] = RC.x()
    k += 1
    simu.nextStep()


# save to disk
np.savetxt('train_tslcp.dat', data_plot)


def plot_results():
    """Plot 3d curve z = f(x, y)
    with ds_state = [x, y, z]
    """
    plt.subplot(211)
    time = data_plot[:, 0]
    q = data_plot[:, 1:number_of_cars + 1]
    v = data_plot[:, number_of_cars + 1:]
    plt.plot(time, q, label='u')
    plt.ylabel('u')
    plt.subplot(212)
    plt.plot(time, v, label='v')
    plt.ylabel('v')
    plt.xlabel('time')
    plt.savefig('train_brakes.png')


# --- Uncomment lines below to plot interesting stuff ---
plot_results()
