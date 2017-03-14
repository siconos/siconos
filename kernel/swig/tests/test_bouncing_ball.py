#!/usr/bin/env python

from siconos.tests_setup import working_dir
import siconos.kernel as sk
import numpy as np
import os


def test_bouncing_ball1():
    """Run a complete simulation (Bouncing ball example)
    LagrangianLinearTIDS,  no plugins.
    """

    t0 = 0.      # start time
    tend = 10.   # end time
    h = 0.005    # time step
    r = 0.1      # ball radius
    g = 9.81     # gravity
    m = 1        # ball mass
    e = 0.9      # restitution coeficient
    theta = 0.5  # theta scheme

    #
    # dynamical system
    #
    x = np.zeros(3, dtype=np.float64)
    x[0] = 1.
    v = np.zeros_like(x)
    # mass matrix
    mass = np.eye(3, dtype=np.float64)
    mass[2, 2] = 3. / 5 * r * r

    # the dynamical system
    ball = sk.LagrangianLinearTIDS(x, v, mass)

    # set external forces
    weight = np.zeros_like(x)
    weight[0] = -m * g
    ball.setFExtPtr(weight)

    #
    # Interactions
    #

    # ball-floor
    H = np.zeros((1, 3), dtype=np.float64)
    H[0, 0] = 1.

    nslaw = sk.NewtonImpactNSL(e)
    relation = sk.LagrangianLinearTIR(H)
    inter = sk.Interaction(nslaw, relation)

    #
    # Model
    #
    bouncing_ball = sk.Model(t0, tend)

    # add the dynamical system to the non smooth dynamical system
    bouncing_ball.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

    # link the interaction and the dynamical system
    bouncing_ball.nonSmoothDynamicalSystem().link(inter, ball)

    #
    # Simulation
    #

    # (1) OneStepIntegrators
    OSI = sk.MoreauJeanOSI(theta)

    # (2) Time discretisation --
    t = sk.TimeDiscretisation(t0, h)

    # (3) one step non smooth problem
    osnspb = sk.LCP()

    # (4) Simulation setup with (1) (2) (3)
    s = sk.TimeStepping(t)
    s.insertIntegrator(OSI)
    s.insertNonSmoothProblem(osnspb)

    # end of model definition

    #
    # computation
    #

    # simulation initialization
    bouncing_ball.setSimulation(s)
    bouncing_ball.initialize()

    #
    # save and load data from xml and .dat
    #
    try:
        from siconos.io import save
        save(bouncing_ball, "bouncingBall.xml")
        save(bouncing_ball, "bouncingBall.bin")

    except:
        print("Warning : could not import save from siconos.io")

    # the number of time steps
    nb_time_steps = int((tend - t0) / h + 1)

    # Get the values to be plotted
    # ->saved in a matrix dataPlot

    data = np.empty((nb_time_steps, 5))

    #
    # numpy pointers on dense Siconos vectors
    #
    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)
    lambda_ = inter.lambda_(1)

    #
    # initial data
    #
    data[0, 0] = t0
    data[0, 1] = q[0]
    data[0, 2] = v[0]
    data[0, 3] = p[0]
    data[0, 4] = lambda_[0]

    k = 1

    # time loop
    while(s.hasNextEvent()):
        s.computeOneStep()

        data[k, 0] = s.nextTime()
        data[k, 1] = q[0]
        data[k, 2] = v[0]
        data[k, 3] = p[0]
        data[k, 4] = lambda_[0]

        k += 1
        #print(s.nextTime())
        s.nextStep()

    #
    # comparison with the reference file
    #

    ref = sk.getMatrix(sk.SimpleMatrix(
        os.path.join(working_dir, "data/result.ref")))
    assert (np.linalg.norm(data - ref) < 1e-12)


def xtest_bouncing_ball_from_xml():

    assert False  # just have to load from xml...


def xtest_bouncing_ball_from_binary():

    assert False  # just have to load from .dat...


def test_bouncing_ball2():
    """Run a complete simulation (Bouncing ball example)
    LagrangianLinearTIDS,  plugged Fext.
    """

    t0 = 0       # start time
    T = 5        # end time
    h = 0.005    # time step
    r = 0.1      # ball radius
    g = 9.81     # gravity
    m = 1        # ball mass
    e = 0.9      # restitution coeficient
    theta = 0.5  # theta scheme

    #
    # dynamical system
    #
    x = np.zeros(3, dtype=np.float64)
    x[0] = 1.
    v = np.zeros_like(x)
    # mass matrix
    mass = np.eye(3, dtype=np.float64)
    mass[2, 2] = 3. / 5 * r * r

    # the dynamical system
    ball = sk.LagrangianLinearTIDS(x, v, mass)
    weight = np.zeros(ball.dimension())
    weight[0] = -m * g
    ball.setFExtPtr(weight)

    # a ball with its own computeFExt
    class Ball(sk.LagrangianLinearTIDS):

        def computeFExt(self, t):
            """External forces operator computation
            """
            print("computing FExt at t=", t)
            #self._fExt[0] = -m * g
            weight = np.zeros(self.dimension())
            weight[0] = -m * g
            self.setFExtPtr(weight)

    ball_d = Ball(x.copy(), v.copy(), mass)
    ball_d.computeFExt(t0)
    # Interactions
    #

    # ball-floor
    H = np.zeros((1, 3), dtype=np.float64)
    H[0, 0] = 1.

    nslaw = sk.NewtonImpactNSL(e)
    nslaw_d = sk.NewtonImpactNSL(e)

    relation = sk.LagrangianLinearTIR(H)
    relation_d = sk.LagrangianLinearTIR(H)

    inter = sk.Interaction(nslaw, relation)
    inter_d = sk.Interaction(nslaw_d, relation_d)

    #
    # Model
    #
    bouncing_ball = sk.Model(t0, T)

    bouncing_ball_d = sk.Model(t0, T)

    # add the dynamical system to the non smooth dynamical system
    bouncing_ball.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)
    bouncing_ball_d.nonSmoothDynamicalSystem().insertDynamicalSystem(ball_d)

    # link the interaction and the dynamical system
    bouncing_ball.nonSmoothDynamicalSystem().link(inter, ball)
    bouncing_ball_d.nonSmoothDynamicalSystem().link(inter_d, ball_d)

    #
    # Simulation
    #

    # (1) OneStepIntegrators
    OSI = sk.MoreauJeanOSI(theta)

    OSI_d = sk.MoreauJeanOSI(theta)

    # (2) Time discretisation --
    t = sk.TimeDiscretisation(t0, h)
    t_d = sk.TimeDiscretisation(t0, h)

    # (3) one step non smooth problem
    osnspb = sk.LCP()

    osnspb_d = sk.LCP()

    # (4) Simulation setup with (1) (2) (3)
    s = sk.TimeStepping(t)
    s.insertIntegrator(OSI)
    s.insertNonSmoothProblem(osnspb)

    s_d = sk.TimeStepping(t_d)
    s_d.insertIntegrator(OSI_d)
    s_d.insertNonSmoothProblem(osnspb_d)

    # end of model definition

    #
    # computation
    #

    # simulation initialization
    bouncing_ball.setSimulation(s)
    bouncing_ball_d.setSimulation(s_d)
    bouncing_ball.initialize()
    bouncing_ball_d.initialize()

    # the number of time steps
    nb_time_steps = int((T - t0) / h + 1)

    # Get the values to be plotted
    # ->saved in a matrix data

    s_d.computeOneStep()

    data = np.empty((nb_time_steps + 1, 5))
    data_d = np.empty((nb_time_steps + 1, 5))

    data[0, 0] = t0
    data[0, 1] = ball.q()[0]
    data[0, 2] = ball.velocity()[0]
    data[0, 3] = ball.p(1)[0]
    data[0, 4] = inter.lambda_(1)

    data_d[0, 0] = t0
    data_d[0, 1] = ball_d.q()[0]
    data_d[0, 2] = ball_d.velocity()[0]
    data_d[0, 3] = ball_d.p(1)[0]
    data_d[0, 4] = inter_d.lambda_(1)

    k = 1

    # time loop
    while(s.hasNextEvent()):
        s.computeOneStep()
        s_d.computeOneStep()

        data[k, 0] = s.nextTime()
        data[k, 1] = ball.q()[0]
        data[k, 2] = ball.velocity()[0]
        data[k, 3] = ball.p(1)[0]
        data[k, 4] = inter.lambda_(1)[0]

        data_d[k, 0] = s_d.nextTime()
        data_d[k, 1] = ball_d.q()[0]
        data_d[k, 2] = ball_d.velocity()[0]
        data_d[k, 3] = ball_d.p(1)[0]
        data_d[k, 4] = inter_d.lambda_(1)[0]
    
        assert np.allclose(data[k, 1], data_d[k, 1])

        #print(s.nextTime())
        k += 1
        s.nextStep()
        s_d.nextStep()
