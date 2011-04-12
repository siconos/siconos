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


    # a ball with its own computeFExt
    class Ball(K.LagrangianLinearTIDS):
        def computeFExt(self, t):
            print "computing FExt at t=",t
            weight = array([-m * g, 0, 0])
            self.setFExtPtr(weight)

    ball_d = Ball(array([1, 0, 0]),array([0, 0, 0]),mass)

    ball_d.setFExtPtr(array([0, 0, 0]))

    #
    # Interactions
    #

    # ball-floor
    H = array([[1,0,0]])

    nslaw = K.NewtonImpactNSL(e)
    nslaw_d = K.NewtonImpactNSL(e)
    
    relation = K.LagrangianLinearTIR(H)
    relation_d = K.LagrangianLinearTIR(H)

    inter = K.Interaction(1, nslaw, relation)
    inter_d = K.Interaction(1, nslaw_d, relation_d)

    #
    # Model
    #
    bouncingBall = K.Model(t0,T)

    bouncingBall_d = K.Model(t0,T)

    # add the dynamical system to the non smooth dynamical system
    bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)
    bouncingBall_d.nonSmoothDynamicalSystem().insertDynamicalSystem(ball_d)

    # link the interaction and the dynamical system
    bouncingBall.nonSmoothDynamicalSystem().link(inter,ball);
    bouncingBall_d.nonSmoothDynamicalSystem().link(inter_d,ball_d);

    


    #
    # Simulation
    #

    # (1) OneStepIntegrators
    OSI = K.Moreau(theta)
    OSI.insertDynamicalSystem(ball)

    OSI_d = K.Moreau(theta)
    OSI_d.insertDynamicalSystem(ball_d)




    # (2) Time discretisation --
    t = K.TimeDiscretisation(t0,h)

    t_d = K.TimeDiscretisation(t0,h)
    
    # (3) one step non smooth problem
    osnspb = K.LCP()

    osnspb_d = K.LCP()
    
    # (4) Simulation setup with (1) (2) (3)
    s = K.TimeStepping(t)
    s.insertIntegrator(OSI)
    s.insertNonSmoothProblem(osnspb)

    s_d = K.TimeStepping(t_d)
    s_d.insertIntegrator(OSI_d)
    s_d.insertNonSmoothProblem(osnspb_d)

    # end of model definition
    
    #
    # computation
    #

    # simulation initialization
    bouncingBall.initialize(s)

    bouncingBall_d.initialize(s_d)


    # the number of time steps
    N = (T-t0)/h

    # Get the values to be plotted 
    # ->saved in a matrix dataPlot

    s_d.computeOneStep()

    dataPlot = empty((N+1,5))
    dataPlot_d = empty((N+1,5))

    lambda_ = inter.lambda_(1)

    dataPlot[0, 0] = t0
    dataPlot[0, 1] = ball.q()[0]
    dataPlot[0, 2] = ball.velocity()[0]
    dataPlot[0, 3] = ball.p(2)[0]
    dataPlot[0, 4] = inter.lambda_(1)

    dataPlot_d[0, 0] = t0
    dataPlot_d[0, 1] = ball_d.q()[0]
    dataPlot_d[0, 2] = ball_d.velocity()[0]
    dataPlot_d[0, 3] = ball_d.p(2)[0]
    dataPlot_d[0, 4] = inter_d.lambda_(1)

    k = 1

    # time loop
    while(s.nextTime() < T):
        s.computeOneStep()
        s_d.computeOneStep()

        dataPlot[k, 0] = s.nextTime()
        dataPlot[k, 1] = ball.q()[0]
        dataPlot[k, 2] = ball.velocity()[0]
        dataPlot[k, 3] = ball.p(2)[0]
        dataPlot[k, 4] = inter.lambda_(1)[0]

        dataPlot_d[k, 0] = s_d.nextTime()
        dataPlot_d[k, 1] = ball_d.q()[0]
        dataPlot_d[k, 2] = ball_d.velocity()[0]
        dataPlot_d[k, 3] = ball_d.p(2)[0]
        dataPlot_d[k, 4] = inter_d.lambda_(1)[0]

        assert dataPlot[k,1] == dataPlot_d[k,1]

        k += 1
        s.nextStep()
        s_d.nextStep()
        print s.nextTime()


def test_C():
    class Rel(K.LagrangianScleronomousR):
        pass
    
    r = Rel()
    j = array([[1,2,3],[4,5,6]])
    r.setJachqPtr(j)

    import numpy as np
    assert np.max(r.C() - array([[1,2,3],[4,5,6]])) == 0.
    assert np.max(r.C() - array([[0,2,3],[4,5,6]])) == 1.
