.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _simulation_nsds:

******************************************
Simulation of non-smooth dynamical systems
******************************************

Once the model (i.e. the nonsmooth dynamical system) has been defined as explained in :ref:`modeling_nsds`,
you must define the global simulation strategy used to formulate and solve this system.
Two types of algorithms are currently available in Siconos:

* *Event-driven algorithms*: the dynamics are decomposed in time into modes, time intervals where the dynamics are smooth, and discrete events, times where the dynamics are nonsmooth.
* *Event-capturing algorithms* (a.k.a. time-stepping): the whole system (smooth dynamics, constraints, nonsmooth laws) is discretised in time, which leads to a one-step nonsmooth problem to be solved at each time step.

The details, advantages and drawbacks of both methods are discussed at length in :cite:`Acary.Brogliato2008`.

For both algorithms, the main steps to describe a simulation are:

* define a time discretisation;
* describe how the dynamical systems are integrated, with 'one-step integrators';
* choose a formulation and a solver for the nonsmooth problem, which gives what we call a 'one-step nonsmooth problem', solved with numerics solvers.

The types of integrators, solvers and formulations of course depend strongly on the strategy.
Before getting into details, here are the minimal steps to build a simulation::

  # define a one-step integrator and associate it with a dynamical system
  osi = MoreauJeanOSI(theta)
  osi.insertDynamicalSystem(your_ds)

  # define a one-step nonsmooth problem
  osnspb = LCP()

  # build a time discretisation
  td = TimeDiscretisation(initial_time, time_step)

  # collect all of them into a global simulation
  simu = TimeStepping(td, osi, osnspb)

  # associate this simulation with a previously defined model (ds and interactions)
  # and initialize
  my_model.setSimulation(simu)
  my_model.initialize()


Depending on your problem, you may have to change the integrator (here a :cpp:class:`MoreauJeanOSI`), the nonsmooth problem formulation (:cpp:class:`LCP`) and the
global strategy (:cpp:class:`TimeStepping`). All the possibilities are detailed in the sections below. You can also look at the examples for
templates.

The simulation loop is then::

  while simu.hasNextEvent():

    # integrate, formalize and solve ...
    simu.computeOneStep()

    # do what you need to save data ...
    dataPlot[k, 0] = s.nextTime()
    dataPlot[k, 1] = q[0]

    # advance to next step
    simu.nextStep()

.. toctree::
   :maxdepth: 4

   time_discretisation
   time_integrators
   event_capturing
   event_driven
   osns_problems
