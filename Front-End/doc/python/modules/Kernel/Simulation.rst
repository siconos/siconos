Simulation
----------

Usage : example of the simulation of the bouncing ball
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the definition of a model for a simple bouncing ball
(:ref:`bouncingball-model`) we are going to run a simulation with a 
Moreau-Jean time stepping scheme (see http://siconos.gforge.inria.fr/UsersGuide/docSimuMoreauTS.html)

We first import the needed classes:

.. testcode::

   from Siconos.Kernel import \
       TimeStepping, MoreauJeanOSI, TimeDiscretisation, LCP


Prior to the construction of a time stepping scheme, we define a time
discretisation by the construction of a TimeDiscretisation object. We
provide the starting time :math:`0` and a fixed time step
:math:`0.05`:

.. testcode::

   td = TimeDiscretisation(0, 0.05)


A time stepping simulation object is built with the time discretisation:

.. testcode::

   simulation = TimeStepping(td)


We have to provide to the time stepping scheme an one-step
integrator. We choose a Moreau-Jean integrator for which we have to
give the value of the :math:`\theta` parameter:

.. testcode::

   OSI = MoreauJeanOSI(0.5)

Then, we attach the previously defined dynamical system to the integrator:

.. testcode::

   bouncingBall.nonSmoothDynamicalSystem().setOSI(ball, OSI)

And we attach this one step integrator to the simulation:

.. testcode::

   simulation.insertIntegrator(OSI)
   
In the Moreau-Jean time-stepping scheme, the unilateral constraints, after
being reformulated at the velocity level, lead at each timestep to nonsmooth
optimization problems. In our case, it is a linear complementarity problem
(LCP):

.. testcode::

   lcp = LCP()

The default solver for LCP is Lemke. As the one step integrator
object, it needs to be attached to the simulation:

.. testcode::

   simulation.insertNonSmoothProblem(lcp)


At this stage, the simulation object does not know the interactions we have
defined for the dynamical systems. An initialization phase remains to be done:

.. testcode::

   bouncingBall.initialize(simulation)

The simulation is now ready for execution.

The simulation object provides methods in order to do the computation at each timestep:

  - `simulation.hasNextEvent()` to check if some computation remains to be done.
  - `simulation.computeOneStep()` to perform the computation a the current timestep.
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



Simulation API
^^^^^^^^^^^^^^

.. automodule:: Siconos.Kernel
  :members: :eval:`under_directory(['../../../Kernel/src/simulationTools'])`
