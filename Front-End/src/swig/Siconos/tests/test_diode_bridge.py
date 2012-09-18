

# need a ref file.
def test_diodebridge1():
    t0 = 0.0
    T = 5.0e-3       # Total simulation time
    h_step = 1.0e-6  # Time step
    Lvalue = 1e-2  # inductance
    Cvalue = 1e-6   # capacitance
    Rvalue = 1e3    # resistance 
    Vinit = 10.0    # initial voltage
    Modeltitle = "DiodeBridge"

    withPlot=True
    if (withPlot) :
        from matplotlib.pyplot import subplot, title, plot, grid, show

    from Siconos.Kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
                               ComplementarityConditionNSL, Interaction,\
                               Model, Moreau, TimeDiscretisation, LCP,  \
                               TimeStepping

    #
    # dynamical system
    #
    init_state = [Vinit,0] 

    A = [[0,          -1.0/Cvalue],
         [1.0/Lvalue, 0          ]]

    LSDiodeBridge=FirstOrderLinearDS(init_state, A)

    #
    # Interactions
    #

    C = [[0.,   0.],
         [0,    0.],
         [-1.,  0.],
         [1.,   0.]]

    D = [[1./Rvalue, 1./Rvalue, -1.,  0.],
         [1./Rvalue, 1./Rvalue,  0., -1.],
         [1.,        0.,         0.,  0.],
         [0.,        1.,         0.,  0.]]

    B = [[0.,        0., -1./Cvalue, 1./Cvalue],
         [0.,        0.,  0.,        0.       ]]

    LTIRDiodeBridge=FirstOrderLinearTIR(C,B)
    LTIRDiodeBridge.setDPtr(D)

    nslaw=ComplementarityConditionNSL(4)
    InterDiodeBridge=Interaction(4, nslaw,LTIRDiodeBridge,1)
    InterDiodeBridge.insert(LSDiodeBridge)

    #
    # Model
    #
    DiodeBridge=Model(t0,T,Modeltitle)

    #   add the dynamical system in the non smooth dynamical system
    DiodeBridge.nonSmoothDynamicalSystem().insertDynamicalSystem(LSDiodeBridge)

    #   link the interaction and the dynamical system
    DiodeBridge.nonSmoothDynamicalSystem().link(InterDiodeBridge,LSDiodeBridge)


    #
    # Simulation
    #

    # (1) OneStepIntegrators
    theta = 0.5
    aOSI = Moreau(LSDiodeBridge,theta)

    # (2) Time discretisation
    aTiDisc = TimeDiscretisation(t0,h_step)

    # (3) Non smooth problem
    aLCP = LCP()

    # (4) Simulation setup with (1) (2) (3)
    aTS = TimeStepping(aTiDisc,aOSI,aLCP)

    # end of model definition

    #
    # computation
    #

    # simulation initialization
    DiodeBridge.initialize(aTS)

    k = 0
    h = aTS.timeStep();
    print "Timestep : ",h
    # Number of time steps
    N = (T-t0)/h
    print "Number of steps : ",N

    # Get the values to be plotted 
    # ->saved in a matrix dataPlot

    from numpy import empty
    dataPlot = empty([N+1,8])

    x = LSDiodeBridge.x()
    print "Initial state : ",x
    y = InterDiodeBridge.y(0)
    print "First y : ",y
    lambda_ = InterDiodeBridge.lambda_(0)

    while (k < N):    
        aTS.computeOneStep()
        #aLCP.display()
        dataPlot[k, 0] = aTS.nextTime()
        #  inductor voltage
        dataPlot[k, 1] = x[0]
        # inductor current
        dataPlot[k, 2] = x[1]
        # diode R1 current
        dataPlot[k, 3] = y[0]
        # diode R1 voltage
        dataPlot[k, 4] = - lambda_[0]
        # diode F2 voltage 
        dataPlot[k, 5] = - lambda_[1]
        # diode F1 current
        dataPlot[k, 6] = lambda_[2]
        # resistor current
        dataPlot[k, 7] = y[0] + lambda_[2]
        k += 1
        aTS.nextStep()


def test_diodebridge2():
    t0 = 0.0
    T = 5.0e-3       # Total simulation time
    h_step = 1.0e-6  # Time step
    Lvalue = 1e-2  # inductance
    Cvalue = 1e-6   # capacitance
    Rvalue = 1e3    # resistance 
    Vinit = 10.0    # initial voltage
    Modeltitle = "DiodeBridge"

    withPlot=True
    if (withPlot) :
        from matplotlib.pyplot import subplot, title, plot, grid, show

    from Siconos.Kernel import FirstOrderLinearDS, FirstOrderLinearR, \
                               ComplementarityConditionNSL, Interaction,\
                               Model, Moreau, TimeDiscretisation, LCP,  \
                               TimeStepping

    #
    # dynamical system
    #
    init_state = [Vinit,0] 

    A = [[0,          -1.0/Cvalue],
         [1.0/Lvalue, 0          ]]

    LSDiodeBridge=FirstOrderLinearDS(init_state, A)

    #
    # Interactions
    #

    C = [[0.,   0.],
         [0,    0.],
         [-1.,  0.],
         [1.,   0.]]

    D = [[1./Rvalue, 1./Rvalue, -1.,  0.],
         [1./Rvalue, 1./Rvalue,  0., -1.],
         [1.,        0.,         0.,  0.],
         [0.,        1.,         0.,  0.]]

    B = [[0.,        0., -1./Cvalue, 1./Cvalue],
         [0.,        0.,  0.,        0.       ]]


    class VoltageSource(FirstOrderLinearR):

        def __init__(self, *args):
            super(VoltageSource, self).__init__(*args)

        def computeE(self, time, inter):
            print time
            print inter.data(Interaction.z)


    LTIRDiodeBridge=VoltageSource(C,B)
    LTIRDiodeBridge.setDPtr(D)

    nslaw=ComplementarityConditionNSL(4)
    InterDiodeBridge=Interaction(4, nslaw,LTIRDiodeBridge,1)
    InterDiodeBridge.insert(LSDiodeBridge)

    #
    # Model
    #
    DiodeBridge=Model(t0,T,Modeltitle)

    #   add the dynamical system in the non smooth dynamical system
    DiodeBridge.nonSmoothDynamicalSystem().insertDynamicalSystem(LSDiodeBridge)

    #   link the interaction and the dynamical system
    DiodeBridge.nonSmoothDynamicalSystem().link(InterDiodeBridge,LSDiodeBridge)


    #
    # Simulation
    #

    # (1) OneStepIntegrators
    theta = 0.5
    aOSI = Moreau(LSDiodeBridge,theta)

    # (2) Time discretisation
    aTiDisc = TimeDiscretisation(t0,h_step)

    # (3) Non smooth problem
    aLCP = LCP()

    # (4) Simulation setup with (1) (2) (3)
    aTS = TimeStepping(aTiDisc,aOSI,aLCP)

    # end of model definition

    #
    # computation
    #

    # simulation initialization
    DiodeBridge.initialize(aTS)

    k = 0
    h = aTS.timeStep();
    print "Timestep : ",h
    # Number of time steps
    N = (T-t0)/h
    print "Number of steps : ",N

    # Get the values to be plotted 
    # ->saved in a matrix dataPlot

    from numpy import empty
    dataPlot = empty([N+1,8])

    x = LSDiodeBridge.x()
    print "Initial state : ",x
    y = InterDiodeBridge.y(0)
    print "First y : ",y
    lambda_ = InterDiodeBridge.lambda_(0)

    while (k < N):    
        aTS.computeOneStep()
        #aLCP.display()
        dataPlot[k, 0] = aTS.nextTime()
        #  inductor voltage
        dataPlot[k, 1] = x[0]
        # inductor current
        dataPlot[k, 2] = x[1]
        # diode R1 current
        dataPlot[k, 3] = y[0]
        # diode R1 voltage
        dataPlot[k, 4] = - lambda_[0]
        # diode F2 voltage 
        dataPlot[k, 5] = - lambda_[1]
        # diode F1 current
        dataPlot[k, 6] = lambda_[2]
        # resistor current
        dataPlot[k, 7] = y[0] + lambda_[2]
        k += 1
        aTS.nextStep()


     #
    # comparison with the reference file
    #
    from Siconos.Kernel import SimpleMatrix, getMatrix
    from numpy.linalg import norm

    ref = getMatrix(SimpleMatrix("diode_bridge.ref"))

    assert (norm(dataPlot - ref) < 1e-12)

