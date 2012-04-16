#!/usr/bin/env python

# this test is taken almost verbatim from RelayBiSimulation_OT2_noCplugin.py
def test_smc1():
    from Siconos.Kernel import FirstOrderLinearDS, Model, TimeDiscretisation, \
            TimeStepping, Moreau, ControlManager, LinearSensor, LinearSMCOT2
    from numpy import eye, empty, zeros
    from math import ceil, sin

    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50*time)
            tmpz = self.z()
            # XXX fix this !
            if len(tmpz) != 2:
                print("DEBUG z has length ", len(tmpz))
                return
            # XXX we need to find a smarter way to do things here
            # we need to convert from vector (sage) to arrayish
            u = [t, -t] + tmpz
            self.setb(u)

    # variable declaration
    ndof = 2   # Number of degrees of freedom of your system
    t0 = 0.0   # start time
    T = 1    # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2 # time step for control
    Xinit = 1.0 # initial position
    theta = 0.5
    N = ceil((T-t0)/h + 10) # number of time steps
    outputSize = 5 # number of variable to store at each time step

    # Matrix declaration
    A = zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
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
    processDS.setb([0, 0])
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
    processIntegrator = Moreau(processDS, theta)
    processSimulation.insertIntegrator(processIntegrator)
    # Actuator, Sensor & ControlManager
    control = ControlManager(process)
    sens = LinearSensor(tSensor, processDS, sensorC, sensorD)

    control.addSensorPtr(sens)
    act = LinearSMCOT2(tActuator, processDS)
    act.setCsurfacePtr(Csurface)
    act.addSensorPtr(sens)
    control.addActuatorPtr(act)

    # Initialization.
    process.initialize(processSimulation)
    control.initialize()
    # This is not working right now
    #eventsManager = s.eventsManager()

    # Matrix for data storage
    dataPlot = empty((3*(N+1), outputSize))
    dataPlot[0, 0] = t0
    dataPlot[0, 1] = processDS.x()[0]
    dataPlot[0, 2] = processDS.x()[1]
    dataPlot[0, 3] = processDS.z()[0]
    dataPlot[0, 4] = processDS.z()[1]

    # Main loop
    k = 1
    while(processSimulation.nextTime() < T):
        processSimulation.computeOneStep()
        dataPlot[k, 0] = processSimulation.nextTime()
        dataPlot[k, 1] = processDS.x()[0]
        dataPlot[k, 2] = processDS.x()[1]
        dataPlot[k, 3] = processDS.z()[0]
        dataPlot[k, 4] = processDS.z()[1]
        k += 1
        processSimulation.nextStep()
    #    print processSimulation.nextTime()
    # Resize matrix
    dataPlot.resize(k, outputSize)

#Same test, but with the simplified interface
def test_smc2():
    from Siconos.Kernel import FirstOrderLinearDS, TimeDiscretisation, \
            ControlFirstOrderLinearDS, LinearSensor, \
            LinearSMCOT2, getMatrix, SimpleMatrix
    from numpy import eye, zeros
    from math import sin
    from numpy.linalg import norm

    # Derive our own version of FirstOrderLinearDS
    class MyFOLDS(FirstOrderLinearDS):
        def computeb(self, time):
            t = sin(50*time)
            tmpz = self.z()
            # XXX fix this !
            if len(tmpz) != 2:
                print("DEBUG z has length ", len(tmpz))
                return
            u = [t, -t] + tmpz
            self.setb(u)

    # variable declaration
    ndof = 2   # Number of degrees of freedom of your system
    t0 = 0.0   # start time
    T = 1    # end time
    h = 1.0e-4  # time step for simulation
    hControl = 1.0e-2 # time step for control
    Xinit = 1.0 # initial position

    # Matrix declaration
    A = zeros((ndof, ndof))
    x0 = [Xinit, -Xinit]
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
    processDS.setb([0, 0])
    controlProcess = ControlFirstOrderLinearDS(t0, T, h, x0, A)
    controlProcess.setProcessDS(processDS)
    controlProcess.initialize()
    # time discretisation
    tSensor = TimeDiscretisation(t0, hControl)
    tActuator = TimeDiscretisation(t0, hControl)
    # Actuator, Sensor & ControlManager
    control = controlProcess.CM()
    sens = LinearSensor(tSensor, processDS, sensorC, sensorD)
    control.addSensorPtr(sens)
    act = LinearSMCOT2(tActuator, processDS)
    act.setCsurfacePtr(Csurface)
    act.addSensorPtr(sens)
    control.addActuatorPtr(act)

    # Initialization
    control.initialize()

    # Run the simulation
    controlProcess.run()
    # get the results
    tmpData = controlProcess.data()
    dataPlot = tmpData
    # compare with the reference
    ref = getMatrix(SimpleMatrix("smc_2.ref"))
    print("%e" % norm(dataPlot - ref))
    if (norm(dataPlot - ref) > 1e-12):
        print(dataPlot - ref)
        print("ERROR: The result is rather different from the reference file.")
