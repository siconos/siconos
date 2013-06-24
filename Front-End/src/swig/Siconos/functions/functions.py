""" Useful set of functions using some siconos functions """

import Siconos.Kernel as SK
import numpy as np


def compute_dt_matrices(A, B, h, TV=False):
    # variable declaration
    t0 = 0.0  # start time
    T = 1  # end time
    n, m = B.shape

     # Matrix declaration
    x0 = np.random.random(n)
    Csurface = np.random.random((m, n))

    # Declaration of the Dynamical System
    if TV:
        processDS = SK.FirstOrderLinearDS(x0, A)
    else:
        processDS = SK.FirstOrderLinearTIDS(x0, A)
    # Model
    process = SK.Model(t0, T)
    process.nonSmoothDynamicalSystem().insertDynamicalSystem(processDS)
    # time discretisation
    processTD = SK.TimeDiscretisation(t0, h)
    # Creation of the Simulation
    processSimulation = SK.TimeStepping(processTD, 0)
    processSimulation.setName("plant simulation")
    # Declaration of the integrator
    processIntegrator = SK.ZeroOrderHold(processDS)
    processSimulation.insertIntegrator(processIntegrator)

    rel = SK.FirstOrderLinearTIR(Csurface, B)
    nslaw = SK.RelayNSL(m)
    inter = SK.Interaction(m, nslaw, rel)

    #process.nonSmoothDynamicalSystem().insertInteraction(inter, True)
    process.nonSmoothDynamicalSystem().link(inter, processDS)
    process.nonSmoothDynamicalSystem().setControlProperty(inter, True)
    # Initialization
    process.initialize(processSimulation)

    # Main loop
    processSimulation.computeOneStep()
    Ad = SK.getMatrix(processIntegrator.Ad(processDS))
    Bd = SK.getMatrix(processIntegrator.Bd(processDS))

    return (Ad, Bd)
