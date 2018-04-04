Simulation
----------

Usage : example of the simulation of the bouncing ball
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the definition of a model for a simple bouncing ball
(:ref:`bouncingball-model`) we are going to run a simulation with a 
Moreau-Jean time stepping scheme (see http://siconos.gforge.inria.fr/UsersGuide/docSimuMoreauTS.html)

We first import the needed classes:

.. testcode::

   from siconos.kernel import \
       TimeStepping, MoreauJeanOSI, TimeDiscretisation, LCP


Prior to the construction of a time stepping scheme, we define a time
discretization by the construction of a TimeDiscretisation object. We
provide the starting time :math:`0` and a fixed time step
:math:`0.05`:

.. testcode::

   td = TimeDiscretisation(0, 0.05)


A time stepping simulation object is built with the time discretization:

.. testcode::

   simulation = TimeStepping(td)


We have to provide to the time stepping scheme an one-step
integrator. We choose a Moreau-Jean integrator for which we have to
give the value of the :math:`\theta` parameter:

.. testcode::

   osi = MoreauJeanOSI(0.5)

Then, we attach the previously defined dynamical system to the integrator, which also adds the integrator to the simulation and prepares required work vectors:

.. testcode::

   simulation.prepareIntegratorForDS(osi, ball, bouncingBall, 0)

In the Moreau-Jean time-stepping scheme, the unilateral constraints, after
being reformulated at the velocity level, lead at each timestep to nonsmooth
optimization problems. In our case, it is a linear complementarity problem
(LCP):

.. testcode::

   lcp = LCP()

The default solver for an LCP is Lemke. As the one step integrator
object, it needs to be attached to the simulation:

.. testcode::

   simulation.insertNonSmoothProblem(lcp)


At this stage, the simulation object does not know the interactions we have
defined for the dynamical systems. An initialization phase remains to be done:

.. testcode::

   bouncingBall.setSimulation(simulation)
   bouncingBall.initialize()

The simulation is now ready for execution.

The simulation object provides methods in order to do the computation at each timestep:

  - `simulation.hasNextEvent()` to check if some computation remains to be done.
  - `simulation.computeOneStep()` to perform the computation at the current timestep.
  - `simulation.nextStep()` to increment the current timestep.


The following loop takes care of running the simulation:

.. testcode::

   while simulation.hasNextEvent():
       simulation.computeOneStep()

       # the current time is simulation.nextTime()
       # the current ball position is ball.q()
       # the current ball velocity is ball.velocity()
       # the current reaction force is ball.p(1)

       simulation.nextStep()


Example of the simulation of a diodes bridge:
+++++++++++++++++++++++++++++++++++++++++++++

The diodes bridge model is given in the Modelisation section :
:ref:`diodes-bridge-model`

The same classes as for the simulation of the bouncing ball are
imported:

.. testcode::

   from siconos.kernel import \
       TimeStepping, MoreauJeanOSI, TimeDiscretisation, LCP

And the time stepping simulation object built with a time discretization
object and the one-step integrator object are the same:

.. testcode::

   td = TimeDiscretisation(0, 0.05)
   simulation = TimeStepping(td)
   osi = MoreauJeanOSI(0.5)

The diodes bridge dynamical system built in the modelisation section is
attached to the integrator which is attached to the simulation:

.. testcode::

   simulation.prepareIntegratorForDS(osi, diodesBridge, DiodesBridgeModel, 0)

At each time step, the Moreau-Jean time stepping scheme needs a linear
complementarity problem to be solved. The `LCP` object is attached to
the simulation.

.. testcode::

   lcp = LCP()
   simulation.insertNonSmoothProblem(lcp)

The simulation is initialized with the modelisation part:

.. testcode::

   DiodesBridgeModel.setSimulation(simulation)
   DiodesBridgeModel.initialize()


It is now ready for execution.

We may set some pointers on interesting values like this:

.. testcode::

   x = diodesBridge.x()
   y = InterDiodeBridge.y(0)
   lambda_ = InterDiodeBridge.lambda_(0)

Then, the execution loop is similar to the one of the bouncing ball:

.. testcode::

   while simulation.hasNextEvent():
       simulation.computeOneStep()

       # inductor voltage is in x[0]
       # inductor current is in x[1]
       # diode R1 current is in y[0]
       # diode R1 voltage is in lambda_[0]
       # diode F2 voltage is in lambda_[1]
       # diode F1 current is in lambda_[2]
       # resistor current is y[0] + lambda_[2]

       simulation.nextStep()


Simulation API
^^^^^^^^^^^^^^

.. automodule:: siconos.kernel
  :members: :eval:`under_directory(['src/simulationTools'])`
