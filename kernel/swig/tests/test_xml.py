"""Test the xml input
"""

try:
    import pytest

    xfail = pytest.mark.xfail
except ImportError:
    import py.test

    xfail = py.test.mark.xfail

from siconos.fromXml import buildModelXML
import siconos.kernel as SK
import numpy as np
import siconos


def test_xml1(datafile):
    """the BouncingBall"""

    bouncingBall, s = buildModelXML(datafile("BBallTS.xml"))

    dsN = SK.dynamicalSystems(bouncingBall.topology().dSG(0))[0].number()
    ball = bouncingBall.dynamicalSystem(dsN)

    N = 2000  # Number of time steps
    # saved in a matrix dataPlot

    outputSize = 4
    dataPlot = np.zeros((N + 1, outputSize))

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)

    dataPlot[0, 0] = bouncingBall.t0()
    dataPlot[0, 1] = q[0]
    dataPlot[0, 2] = v[0]
    dataPlot[0, 3] = p[0]

    print("====> Start computation ...")
    # --- Time loop  ---
    k = 1
    while s.hasNextEvent():
        s.computeOneStep()
        # --- Get values to be plotted ---
        dataPlot[k, 0] = s.nextTime()
        dataPlot[k, 1] = q[0]
        dataPlot[k, 2] = v[0]
        dataPlot[k, 3] = p[0]
        s.nextStep()
        k += 1

    print("End of computation - Number of iterations done: {:}".format(k))
    print("====> Output file writing ...")
    dataPlot.resize(k, outputSize)
    # np.savetxt("BBallTS.dat", dataPlot)

    # Comparison with a reference file
    dataPlotRef = SK.getMatrix(
        SK.SimpleMatrix(datafile("BBallTSXML.ref"))
    )
    if np.linalg.norm(dataPlot - dataPlotRef, ord=np.inf) > 1e-12:
        print(dataPlot - dataPlotRef)
        print("ERROR: The result is rather different from the reference file.")


def test_xml2(datafile):
    """BallInBowl"""
    # --- buildModelXML loading from xml file ---
    bouncingBall, s = buildModelXML(datafile("BallInBowl.xml"))

    # --- Get the simulation ---
    k = 0
    T = bouncingBall.finalT()
    t0 = bouncingBall.t0()
    h = s.timeStep()
    N = int((T - t0) / h)

    # --- Get the values to be plotted ---
    # . saved in a matrix dataPlot
    dataPlot = np.zeros((N + 1, 6))

    print("Prepare data for plotting ... ")
    # For the initial time step:
    # time
    dataPlot[k, 0] = bouncingBall.t0()
    # state q for the first dynamical system (ball)
    dsN = SK.dynamicalSystems(bouncingBall.topology().dSG(0))[0].number()
    ball = bouncingBall.dynamicalSystem(dsN)

    q = ball.q()
    v = ball.velocity()
    p = ball.p(1)

    dataPlot[k, 1] = q[0]
    dataPlot[k, 2] = v[0]
    dataPlot[k, 3] = q[1]
    dataPlot[k, 4] = v[1]
    dataPlot[k, 5] = p[0]

    # --- Compute elapsed time ---
    print("Computation ... ")
    # --- Time loop  ---
    while s.hasNextEvent():
        # solve ...
        s.computeOneStep()

        # --- Get values to be plotted ---
        # time
        dataPlot[k, 0] = s.nextTime()
        # Ball: state q
        dataPlot[k, 1] = q[0]
        # Ball: velocity
        dataPlot[k, 2] = v[0]
        # Ground: q
        dataPlot[k, 3] = q[1]
        # Ground: velocity
        dataPlot[k, 4] = v[1]
        # Reaction
        dataPlot[k, 5] = p[0]
        #  dataPlot[k, 6] = osi.computeResidu()
        # transfer of state i+1 into state i and time incrementation
        s.nextStep()
        k += 1

    # Number of time iterations
    print("Number of iterations done: ")

    # dataPlot (ascii) output
    # ioMatrix::write(dataPlot,"noDim")
    # np.savetxt("BallInBowl.dat", dataPlot)


def test_xml3(datafile):
    """DryFriction"""
    # --- buildModelXML loading from xml file ---
    oscillator, s = buildModelXML(datafile("DryFriction.xml"))

    # --- Get the simulation ---
    k = 0
    T = oscillator.finalT()
    t0 = oscillator.t0()
    h = s.timeStep()
    N = int((T - t0) / h)
    # --- Get the values to be plotted ---
    # . saved in a matrix dataPlot
    dataPlot = np.zeros((N + 1, 5))

    print("Prepare data for plotting ... ")
    # For the initial time step:
    # time
    dataPlot[k, 0] = t0
    # state q for the first dynamical system (ball)
    dsN = SK.dynamicalSystems(oscillator.topology().dSG(0))[0].number()
    oscillo = oscillator.dynamicalSystem(dsN)
    inter = SK.interactions(oscillator.topology().indexSet(0))[0]

    dataPlot[k, 1] = oscillo.q()[0]
    # velocity for the oscillo
    dataPlot[k, 2] = oscillo.velocity()[0]
    dataPlot[k, 3] = inter.lambda_(1)[0]
    dataPlot[k, 4] = inter.lambda_(1)[1]

    # --- Compute elapsed time ---
    print("Computation ... ")
    # --- Time loop  ---
    while k < N:
        # get current time step
        k += 1
        #  print( " Pas " << k
        # solve ...
        s.computeOneStep()
        # --- Get values to be plotted ---
        # time
        dataPlot[k, 0] = s.nextTime()
        # Oscillo: state q
        dataPlot[k, 1] = oscillo.q()[0]
        # Oscillo: velocity
        dataPlot[k, 2] = oscillo.velocity()[0]

        dataPlot[k, 3] = inter.lambda_(1)[0]
        dataPlot[k, 4] = inter.lambda_(1)[1]
        # transfer of state i+1 into state i and time incrementation
        s.nextStep()

    # Number of time iterations
    print("Number of iterations done: {:}".format(k))

    # dataPlot (ascii) output
    # np.savetxt("DryFriction.dat", dataPlot)


@xfail
def test_xml4(datafile):
    """CamFollower"""
    # --- buildModelXML loading from xml file ---
    CamFollower, s = buildModelXML(datafile("CamFollower_TIDS.xml"))

    # --- Get and initialize the simulation ---
    k = 0
    T = CamFollower.finalT()
    t0 = CamFollower.t0()
    h = s.timeStep()
    N = int((T - t0) / h)

    # --- Get the values to be plotted ---
    # . saved in a matrix dataPlot
    dataPlot = np.zeros((N + 1, 8))

    print("Prepare data for plotting ... ")
    # For the initial time step:
    # time
    dataPlot[k, 0] = t0

    # state q for the Follower
    dsN = CamFollower.topology().dSG(0).dynamicalSystems()[0].number()
    Follower = CamFollower.dynamicalSystem(dsN)
    inter = CamFollower.topology().dSG(0).interactions()[0]
    # Position of the Follower
    dataPlot[k, 1] = Follower.q()[0]
    # Velocity for the Follower
    dataPlot[k, 2] = Follower.velocity()[0]
    # Reaction
    dataPlot[k, 3] = inter.lambda_(1)[0]
    # External Forcing
    dataPlot[k, 4] = Follower.fExt()[0]

    # State of the Cam
    rpm = 358

    CamEqForce = CamState(t0, rpm, CamPosition, CamVelocity, CamAcceleration)
    # Position of the Cam
    dataPlot[k, 5] = CamPosition
    # Velocity of the Cam
    dataPlot[k, 6] = CamVelocity
    # Acceleration of the Cam
    dataPlot[k, 7] = CamPosition + Follower.q()[0]

    print("Computation ... ")
    # --- Time loop  ---

    while k < N:
        # get current time step
        k += 1

        s.computeOneStep()
        # --- Get values to be plotted ---
        dataPlot[k, 0] = s.nextTime()
        #  dataPlot[k, 1] = Follower.q()[0]
        #  dataPlot[k, 2] = ball.velocity()[0]
        dataPlot[k, 1] = Follower.q()[0]
        dataPlot[k, 2] = Follower.velocity()[0]
        dataPlot[k, 3] = inter.lambda_(1)[0]
        dataPlot[k, 4] = Follower.fExt()[0]

        CamEqForce = CamState(
            s.nextTime(), rpm, CamPosition, CamVelocity, CamAcceleration
        )

        dataPlot[k, 5] = CamPosition
        dataPlot[k, 6] = CamVelocity
        dataPlot[k, 7] = CamPosition + Follower.q()[0]
        # transfer of state i+1 into state i and time incrementation

        s.nextStep()

    # Number of time iterations
    print("Number of iterations done: {:}".format(k))

    # dataPlot (ascii) output
    # np.savetxt("CamFollower.dat", dataPlot)


if siconos.WITH_FORTRAN:

    def test_xml5(datafile):
        """Bouncing Ball ED"""
        # --- buildModelXML loading from xml file ---
        bouncingBall, s = buildModelXML(datafile("BBallED.xml"))

        # --- Get and initialize the simulation ---
        dsN = SK.dynamicalSystems(bouncingBall.topology().dSG(0))[0].number()
        ball = bouncingBall.dynamicalSystem(dsN)

        # --- Get the values to be plotted ---
        # . saved in a matrix dataPlot

        N = 12368  # Number of saved points: depends on the number of events ...
        outputSize = 5
        dataPlot = np.zeros((N + 1, outputSize))

        q = ball.q()
        v = ball.velocity()

        dataPlot[0, 0] = bouncingBall.t0()
        dataPlot[0, 1] = q[0]
        dataPlot[0, 2] = v[0]
        dataPlot[0, 3] = 0.0
        dataPlot[0, 4] = 0.0

        print("====> Start computation ... ")
        # --- Time loop  ---
        eventsManager = s.eventsManager()
        numberOfEvent = 0
        k = 0
        nonSmooth = False

        while s.hasNextEvent():
            k += 1

            s.advanceToEvent()
            p = ball.p(1)
            f = ball.p(2)

            if eventsManager.nextEvent().getType() == 2:
                nonSmooth = True

            s.processEvents()
            # If the treated event is non smooth, the pre-impact state has been solved
            # in memory vectors during process.
            if nonSmooth:
                dataPlot[k, 0] = s.startingTime()
                dataPlot[k, 1] = ball.qMemory().getSiconosVector(1)[0]
                dataPlot[k, 2] = ball.velocityMemory().getSiconosVector(1)[0]
                k += 1
                nonSmooth = False
            dataPlot[k, 0] = s.startingTime()
            dataPlot[k, 1] = q[0]
            dataPlot[k, 2] = v[0]
            dataPlot[k, 3] = p[0]
            dataPlot[k, 4] = f[0]
            numberOfEvent += 1

        # --- Output files ---
        dataPlot.resize(k, outputSize)
        # np.savetxt("BBallED.dat", dataPlot)
        # Comparison with a reference file
        dataPlotRef = SK.getMatrix(
            SK.SimpleMatrix(datafile("BouncingBallEDXml.ref")))

        if np.linalg.norm(dataPlot - dataPlotRef, ord=np.inf) > 1e-11:
            print("Warning. The results is rather different from the reference file.")
            print(np.linalg.norm(dataPlot - dataPlotRef, ord=np.inf))
            exit(1)


def test_xml6(datafile):
    """BeadPlan"""
    # --- buildModelXML loading from xml file ---
    oscillator, s = buildModelXML(datafile("BeadPlan.xml"))

    # --- Get and initialize the simulation ---

    k = 0
    t0 = oscillator.t0()
    T = oscillator.finalT()
    h = s.timeStep()
    N = int((T - t0) / h)  # Number of time steps

    # --- Get the values to be plotted ---
    # . saved in a matrix dataPlot
    dataPlot = np.zeros((N, 3))

    print("Prepare data for plotting ... ")
    # For the initi)al time step:
    # XXX fix this crap
    dsN = SK.dynamicalSystems(oscillator.topology().dSG(0))[0].number()
    oscillo = oscillator.dynamicalSystem(dsN)
    q = oscillo.q()
    v = oscillo.velocity()
    p = oscillo.p(1)

    dataPlot[k, 0] = t0
    dataPlot[0, 1] = q[0]
    # velocity for the oscillo
    dataPlot[0, 2] = v[0]

    print("Computation ... ")
    # --- Time loop  ---
    # while (s.hasNextEvent()
    # {
    #   // solve ...
    #   s.computeOneStep()
    #   // --- Get values to be plotted ---
    #   dataPlot[k, 0] = s.nextTime()
    #   dataPlot[k, 1] = oscillo.q()[0]
    #   dataPlot[k, 2] = oscillo.velocity()[0]
    #   // transfer of state i+1 into state i and time incrementation
    #   s.nextStep()
    #   k += 1
    #   }

    # Number of time iterations
    print("Number of iterations done: {:}".format(k))

    # dataPlot (ascii) output
    # np.savetxt("BeadPlan.dat", dataPlot)


if __name__ == "__main__":
    print("test_xml1")
    test_xml1()
