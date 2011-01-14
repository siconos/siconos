#!/usr/bin/env python

from numpy import *
import Siconos.Kernel as K

def test_bouncing_ball():
    t0 = 0      # start time
    T = 5      # end time
    h = 0.005   # time step
    r = 0.1     # ball radius
    g = 9.81    # gravity
    m = 1       # ball mass
    e = 0.9     # restitution coeficient
    theta = 0.5 # theta scheme



    #
    # dynamical system
    #
    x = array([1,0,0]) # initial position
    v = array([0,0,0]) # initial velocity
    mass = eye(3)      # mass matrix
    mass[2,2]=3./5 * r * r
                    
    # the dynamical system
    ball = K.LagrangianLinearTIDS(x,v,mass)

    # set external forces 
    weight = array([-m * g, 0, 0])
    ball.setFExtPtr(weight)

    #
    # Interactions
    #

    # ball-floor
    H = array([[1,0,0]])

    nslaw = K.NewtonImpactNSL(e)
    relation = K.LagrangianLinearTIR(H)
    inter = K.Interaction(1, nslaw, relation)

    #
    # Model
    #
    bouncingBall = K.Model(t0,T)

    # add the dynamical system to the non smooth dynamical system
    bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

    # link the interaction and the dynamical system
    bouncingBall.nonSmoothDynamicalSystem().link(inter,ball);


    #
    # Simulation
    #

    # (1) OneStepIntegrators
    OSI = K.Moreau(theta)
    OSI.insertDynamicalSystem(ball)

    # (2) Time discretisation --
    t = K.TimeDiscretisation(t0,h)
    
    # (3) one step non smooth problem
    osnspb = K.LCP()
    
    # (4) Simulation setup with (1) (2) (3)
    s = K.TimeStepping(t)
    s.insertIntegrator(OSI)
    s.insertNonSmoothProblem(osnspb)

    # end of model definition
    
    #
    # computation
    #

    # simulation initialization
    bouncingBall.initialize(s)


    # the number of time steps
    N = (T-t0)/h

    # Get the values to be plotted 
    # ->saved in a matrix dataPlot

    dataPlot = empty((N+1,5))

    lambda_ = inter.lambda_(1)

    dataPlot[0, 0] = t0
    dataPlot[0, 1] = ball.q()[0]
    dataPlot[0, 2] = ball.velocity()[0]
    dataPlot[0, 3] = ball.p(2)[0]
    dataPlot[0, 4] = inter.lambda_(1)

    k = 1

    # time loop
    while(s.nextTime() < T):
        s.computeOneStep()

        dataPlot[k, 0] = s.nextTime()
        dataPlot[k, 1] = ball.q()[0]
        dataPlot[k, 2] = ball.velocity()[0]
        dataPlot[k, 3] = ball.p(2)[0]
        dataPlot[k, 4] = inter.lambda_(1)[0]

        k += 1
        s.nextStep()
        print s.nextTime()


