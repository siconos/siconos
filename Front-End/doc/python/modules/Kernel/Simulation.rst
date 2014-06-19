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

   OSI.insertDynamicalSystem(ball)

And we attach this one step integrator to the simulation:

.. testcode::

   simulation.insertIntegrator(OSI)
   

In the Moreau-Jean time-stepping scheme, the unilateral constraints,
after being reformulated in terms of velocity, lead at each time step
to nonsmooth optimization problems. In the case of a bouncing ball
without friction (a constraint space with one dimension), the problem
is a linear complementarity problem (LCP):

.. testcode::

   lcp = LCP()

This one step nonsmooth problem object carry the default settings for
the LCP solver, which may be modified. As the one step integrator
object, it needs to be attached to the simulation:

.. testcode::

   simulation.insertNonSmoothProblem(lcp)


At this stage, the simulation object does not know the relations we
have defined for the dynamical systems. An initialization by the model
needs to be done:

.. testcode::

   bouncingBall.initialize(simulation)

The simulation can now be executed by a loop:

.. testcode::

   while(simulation.hasNextEvent()):
       simulation.computeOneStep()

       # the current time is simulation.nextTime()
       # the current ball position is ball.q()
       # the current ball velocity is ball.velocity()
       # the current reaction force is ball.p(1)

       s.nextStep()



Simulation API
^^^^^^^^^^^^^^

.. automodule:: Siconos.Kernel
  :members: :eval:`under_directory(['../../../Kernel/src/simulationTools'])`
