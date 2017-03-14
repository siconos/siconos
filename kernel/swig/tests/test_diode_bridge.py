"""Diode Bridge simulation. See examples manual.
Test purpose --> compare with reference results.
"""
import os
from siconos.tests_setup import working_dir
from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
    ComplementarityConditionNSL, Interaction, Model, EulerMoreauOSI, \
    TimeDiscretisation, LCP, TimeStepping
from numpy import empty
from siconos.kernel import SimpleMatrix, getMatrix
from numpy.linalg import norm


# Common data
t0 = 0.0
total_time = 5.0e-3
time_step = 1.0e-6
inductance = 1e-2
capacitance = 1e-6
resistance = 1e3
initial_voltage = 10.0
model_title = "DiodeBridge"
init_state = [initial_voltage, 0]

A = [[0, -1.0 / capacitance], [1.0 / inductance, 0]]

C = [[0., 0.], [0, 0.], [-1., 0.], [1., 0.]]

D = [[1. / resistance, 1. / resistance, -1., 0.],
     [1. / resistance, 1. / resistance, 0., -1.],
     [1., 0., 0., 0.],
     [0., 1., 0., 0.]]

B = [[0., 0., -1. / capacitance, 1. / capacitance],
     [0., 0., 0., 0.]]


def test_diode_bridge():
    """Build diode bridge model"""
    # dynamical system
    bridge_ds = FirstOrderLinearDS(init_state, A)
    # interaction
    diode_bridge_relation = FirstOrderLinearTIR(C, B)
    diode_bridge_relation.setDPtr(D)

    nslaw = ComplementarityConditionNSL(4)
    bridge_interaction = Interaction(nslaw, diode_bridge_relation)

    # Model
    diode_bridge = Model(t0, total_time, model_title)

    #  add the dynamical system in the non smooth dynamical system
    diode_bridge.nonSmoothDynamicalSystem().insertDynamicalSystem(bridge_ds)

    #   link the interaction and the dynamical system
    diode_bridge.nonSmoothDynamicalSystem().link(bridge_interaction, bridge_ds)

    # Simulation

    # (1) OneStepIntegrators
    theta = 0.5
    integrator = EulerMoreauOSI(theta)
    # (2) Time discretisation
    time_discretisation = TimeDiscretisation(t0, time_step)

    # (3) Non smooth problem
    non_smooth_problem = LCP()

    # (4) Simulation setup with (1) (2) (3)
    bridge_simulation = TimeStepping(time_discretisation,
                                     integrator, non_smooth_problem)

    # simulation initialization
    diode_bridge.setSimulation(bridge_simulation)
    diode_bridge.initialize()
    k = 0
    h = bridge_simulation.timeStep()
    # Number of time steps
    N = int((total_time - t0) / h)

    # Get the values to be plotted
    # ->saved in a matrix dataPlot
    data_plot = empty([N, 8])

    x = bridge_ds.x()
    print("Initial state : ", x)
    y = bridge_interaction.y(0)
    print("First y : ", y)
    lambda_ = bridge_interaction.lambda_(0)

    # For the initial time step:
    # time
    data_plot[k, 0] = t0

    #  inductor voltage
    data_plot[k, 1] = x[0]

    # inductor current
    data_plot[k, 2] = x[1]

    # diode R1 current
    data_plot[k, 3] = y[0]

    # diode R1 voltage
    data_plot[k, 4] = - lambda_[0]

    # diode F2 voltage
    data_plot[k, 5] = - lambda_[1]

    # diode F1 current
    data_plot[k, 6] = lambda_[2]

    # resistor current
    data_plot[k, 7] = y[0] + lambda_[2]

    k += 1
    while k < N:
        bridge_simulation.computeOneStep()
        #non_smooth_problem.display()
        data_plot[k, 0] = bridge_simulation.nextTime()
        #  inductor voltage
        data_plot[k, 1] = x[0]
        # inductor current
        data_plot[k, 2] = x[1]
        # diode R1 current
        data_plot[k, 3] = y[0]
        # diode R1 voltage
        data_plot[k, 4] = - lambda_[0]
        # diode F2 voltage
        data_plot[k, 5] = - lambda_[1]
        # diode F1 current
        data_plot[k, 6] = lambda_[2]
        # resistor current
        data_plot[k, 7] = y[0] + lambda_[2]
        k += 1
        bridge_simulation.nextStep()

    #
    # comparison with the reference file
    #
    ref = getMatrix(SimpleMatrix(os.path.join(working_dir,
                                              "data/diode_bridge.ref")))
    assert norm(data_plot - ref) < 1e-12
    return ref, data_plot
