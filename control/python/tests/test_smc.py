# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.
#
# Copyright 2026 INRIA.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# this test is taken almost ve@rbatim from RelayBiSimulation_OT2_noCplugin.py

import siconos.modeling as sm
import numpy as np
from math import ceil, sin
import siconos.integrators
import siconos.simulation
import scipy.linalg as la


def test_smc1():
    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(sm.FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50 * time)
            # XXX fix this !
            u = [t, -t]
            self.setbPtr(u)

    # variable declaration
    ndof = 2  # Number of degrees of freedom of your system
    t0 = 0.0  # start time
    T = 1  # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2  # time step for control
    Xinit = 1.0  # initial position
    N = int(ceil((T - t0) / h + 10))  # number of time steps
    outputSize = 4  # number of variable to store at each time step

    # Matrix declaration
    A = np.zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
    Brel = np.array([[0], [1]])
    sensorC = np.eye(ndof)
    sensorD = np.zeros((ndof, ndof))
    Csurface = [[0, 1.0]]

    # Simple check
    if h > hControl:
        print("hControl must be bigger than h")
        exit(1)

    # Declaration of the Dynamical System
    processDS = MyFOLDS(x0, A)
    # XXX b is not automatically created ...
    #    processDS.setb([0, 0])
    # Model
    process = sm.NonSmoothDynamicalSystem(t0, T)
    process.insertDynamicalSystem(processDS)
    # time discretization
    processTD = siconos.simulation.TimeDiscretisation(t0, h)
    tSensor = siconos.simulation.TimeDiscretisation(t0, hControl)
    tActuator = siconos.simulation.TimeDiscretisation(t0, hControl)
    # Creation of the Simulation
    processSimulation = siconos.simulation.TimeStepping(process, processTD, 0)
    processSimulation.setName("plant simulation")

    # Declaration of the integrator
    processIntegrator = siconos.integrators.ZeroOrderHoldOSI()
    processSimulation.associate(processIntegrator, processDS)
    # Actuator, Sensor & ControlManager
    control = siconos.control.ControlManager(processSimulation)
    sens = siconos.control.LinearSensor(processDS, sensorC, sensorD)

    control.addSensorPtr(sens, tSensor)
    act = siconos.control.LinearSMCOT2(sens)
    act.setCsurface(Csurface)
    act.setB(Brel)
    control.addActuatorPtr(act, tActuator)

    # Initialization.
    control.initialize(process)
    # This is not working right now
    # eventsManager = s.eventsManager()

    # Matrix for data storage
    dataPlot = np.empty((3 * (N + 1), outputSize))
    dataPlot[0, 0] = t0
    dataPlot[0, 1] = processDS.x()[0]
    dataPlot[0, 2] = processDS.x()[1]
    dataPlot[0, 3] = act.u()[0]

    # Main loop
    k = 1
    while processSimulation.hasNextEvent():
        if (
            processSimulation.eventsManager().nextEvent().getType()
            == siconos.simulation.TD_EVENT
        ):
            processSimulation.computeOneStep()
        dataPlot[k, 0] = processSimulation.nextTime()
        dataPlot[k, 1] = processDS.x()[0]
        dataPlot[k, 2] = processDS.x()[1]
        dataPlot[k, 3] = act.u()[0]
        k += 1
        processSimulation.nextStep()
    #    print processSimulation.nextTime()
    # Resize matrix
    ### dataPlot.resize(k, outputSize)


# Same test, but with the simplified interface
def test_smc2(datafile):  # uses datafile pytest fixture

    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(sm.FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50 * time)
            u = [t, -t]
            self.setbPtr(u)

    # variable declaration
    ndof = 2  # Number of degrees of freedom of your system
    t0 = 0.0  # start time
    T = 1  # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2  # time step for control
    Xinit = 1.0  # initial position

    # Matrix declaration
    A = np.zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
    sensorC = np.eye(ndof)
    sensorD = np.zeros((ndof, ndof))
    Brel = np.array([[0], [1]])
    Csurface = [[0, 1.0]]

    # Simple check
    if h > hControl:
        print("hControl must be bigger than h")
        exit(1)

    # Declaration of the Dynamical System
    processDS = MyFOLDS(x0, A)
    # XXX b is not automatically created ...
    processDS.setbPtr([0, 0])
    sim = siconos.control.ControlZOHSimulation(t0, T, h)
    sim.addDynamicalSystem(processDS)
    # time discretisation
    # Actuator, Sensor & ControlManager
    sens = siconos.control.LinearSensor(processDS, sensorC, sensorD)
    sim.addSensor(sens, hControl)
    act = siconos.control.LinearSMCOT2(sens)
    act.setCsurface(Csurface)
    act.setB(Brel)
    sim.addActuator(act, hControl)

    # Run the simulation
    sim.initialize()
    sim.run()
    # get the results
    tmpData = sim.data()
    dataPlot = tmpData
    # compare with the reference
    ref = np.loadtxt(datafile("smc_2.ref.gz"))
    np.savetxt("smc_2.dat", dataPlot)
    print("%e" % la.norm(dataPlot - ref))
    if la.norm(dataPlot - ref) > 5e-12:
        print(dataPlot - ref)
        print("ERROR: The result is rather different from the reference file.")
    assert la.norm(dataPlot - ref) < 5e-12


if __name__ == "__main__":
    print("test_smc1")
    test_smc1()
    test_smc2()
