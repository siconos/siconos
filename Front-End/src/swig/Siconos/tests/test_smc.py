#!/usr/bin/env python


# this test is taken almost verbatim from RelayBiSimulation_OT2_noCplugin.py
def test_smc1():
    from Siconos.Kernel import FirstOrderLinearDS, Model, TimeDiscretisation, \
        TimeStepping, ZeroOrderHold, ControlManager, LinearSensor, LinearSMCOT2, \
        TD_EVENT
    from numpy import eye, empty, zeros
    import numpy as np
    from math import ceil, sin

    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50*time)
            # XXX fix this !
            u = [t, -t]
            self.setb(u)

    # variable declaration
    ndof = 2   # Number of degrees of freedom of your system
    t0 = 0.0   # start time
    T = 1    # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2  # time step for control
    Xinit = 1.0  # initial position
    N = ceil((T-t0)/h + 10)  # number of time steps
    outputSize = 4  # number of variable to store at each time step

    # Matrix declaration
    A = zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
    Brel = np.array([[0], [1]])
    sensorC = eye(ndof)
    sensorD = zeros((ndof, ndof))
    Csurface = [[0, 1.0]]

    # Simple check
    if h > hControl:
        print "hControl must be bigger than h"
        exit(1)

    # Declaration of the Dynamical System
    processDS = MyFOLDS(x0, A)
    # XXX b is not automatically created ...
#    processDS.setb([0, 0])
    # Model
    process = Model(t0, T)
    process.nonSmoothDynamicalSystem().insertDynamicalSystem(processDS)
    # time discretization
    processTD = TimeDiscretisation(t0, h)
    tSensor = TimeDiscretisation(t0, hControl)
    tActuator = TimeDiscretisation(t0, hControl)
    # Creation of the Simulation
    processSimulation = TimeStepping(processTD, 0)
    processSimulation.setName("plant simulation")
    # Declaration of the integrator
    processIntegrator = ZeroOrderHold(processDS)
    processSimulation.insertIntegrator(processIntegrator)
    # Actuator, Sensor & ControlManager
    control = ControlManager(processSimulation)
    sens = LinearSensor(processDS, sensorC, sensorD)

    control.addSensorPtr(sens, tSensor)
    act = LinearSMCOT2(sens)
    act.setCsurfacePtr(Csurface)
    act.setBPtr(Brel)
    control.addActuatorPtr(act, tActuator)

    # Initialization.
    process.initialize(processSimulation)
    control.initialize(process)
    # This is not working right now
    #eventsManager = s.eventsManager()

    # Matrix for data storage
    dataPlot = empty((3*(N+1), outputSize))
    dataPlot[0, 0] = t0
    dataPlot[0, 1] = processDS.x()[0]
    dataPlot[0, 2] = processDS.x()[1]
    dataPlot[0, 3] = act.u()[0]

    # Main loop
    k = 1
    while processSimulation.hasNextEvent():
        if processSimulation.eventsManager().nextEvent().getType() == TD_EVENT:
            processSimulation.computeOneStep()
        dataPlot[k, 0] = processSimulation.nextTime()
        dataPlot[k, 1] = processDS.x()[0]
        dataPlot[k, 2] = processDS.x()[1]
        dataPlot[k, 3] = act.u()[0]
        k += 1
        processSimulation.nextStep()
    #    print processSimulation.nextTime()
    # Resize matrix
    dataPlot.resize(k, outputSize)


#Same test, but with the simplified interface
def test_smc2():
    from Siconos.Kernel import FirstOrderLinearDS, TimeDiscretisation, \
        ControlFirstOrderLinearS, LinearSensor, ioMatrix_write, \
        LinearSMCOT2, getMatrix, SimpleMatrix
    from numpy import eye, zeros
    import numpy as np
    from math import sin
    from numpy.linalg import norm

    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50*time)
            u = [t, -t]
            self.setb(u)

    # variable declaration
    ndof = 2   # Number of degrees of freedom of your system
    t0 = 0.0   # start time
    T = 1    # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2  # time step for control
    Xinit = 1.0  # initial position

    # Matrix declaration
    A = zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
    sensorC = eye(ndof)
    sensorD = zeros((ndof, ndof))
    Brel = np.array([[0], [1]])
    Csurface = [[0, 1.0]]

    # Simple check
    if h > hControl:
        print "hControl must be bigger than h"
        exit(1)

    # Declaration of the Dynamical System
    processDS = MyFOLDS(x0, A)
    # XXX b is not automatically created ...
    processDS.setb([0, 0])
    controlProcess = ControlFirstOrderLinearS(t0, T, h, x0, A)
    controlProcess.setProcessDS(processDS)
    controlProcess.initialize()
    # time discretisation
    tSensor = TimeDiscretisation(t0, hControl)
    tActuator = TimeDiscretisation(t0, hControl)
    # Actuator, Sensor & ControlManager
    control = controlProcess.CM()
    sens = LinearSensor(processDS, sensorC, sensorD)
    control.addSensorPtr(sens, tSensor)
    act = LinearSMCOT2(sens)
    act.setCsurfacePtr(Csurface)
    act.setBPtr(Brel)
    control.addActuatorPtr(act, tActuator)

    # Initialization
    control.initialize(controlProcess.model())

    # Run the simulation
    controlProcess.run()
    # get the results
    tmpData = controlProcess.data()
    dataPlot = tmpData
    # compare with the reference
    ref = getMatrix(SimpleMatrix("smc_2.ref"))
    np.savetxt("smc_2.dat", dataPlot)
    print("%e" % norm(dataPlot - ref))
    if (norm(dataPlot - ref) > 5e-12):
        print(dataPlot - ref)
        print("ERROR: The result is rather different from the reference file.")
#    assert norm(dataPlot - ref) < 5e-12
