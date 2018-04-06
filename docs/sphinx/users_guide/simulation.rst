.. _simulation_nsds:

******************************************
Simulation of non-smooth dynamical systems
******************************************

Once the model (i.e. the nonsmooth dynamical system) has been properly defined and described as explained in :ref:`modeling_nsds`,
the global simulation strategy to formalize and solve this system must be defined.
At the time, two types of algorithms are available in siconos:

* *Event-driven algorithms*: based on the time-decomposition of the dynamics in modes, time-intervals where the dynamics are smooth, and discrete events, times where the dynamics are nonsmooth.
* *Event-capturing algorithms* (a.k.a time-stepping), where a time-discretisation of the whole system (smooth dynamics, constraints, nonsmooth laws) is written leading to a one-step nonsmooth problem that must be solved at each time step.

Details, advantages and drawbacks of both methods are largely discussed in :cite:`Acary.Brogliato.2008`.

For both algorithms, the main steps to describe a simulation are:

* define a time discretisation
* describe how dynamical systems will be integrated, thanks to 'one-step integrators'
* choose a formalisation and a solver for the nonsmooth problem, which leads to what we call 'one-step nonsmooth problem' based on numerics solvers.

The types of integrators, solvers, formulation obviously strongly depend on the strategy.
To clarify things before getting into details, here are the standard minimal steps to write to build a simulation::

  # define a one-step integrator and associate it to a dynamical system
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


Depending on your problem, you may have to change the integrator (here a :class:`MoreauJeanOSI`), the nonsmooth problem formulation (:class:`LCP`) and the
global strategy (:class:`TimeStepping`). Details on all the possibilities will be given in the sections below. You may also check the examples package to find some
templates.

Then, the simulation loop will be::

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
	      
   event_capturing
   event_driven
   time_discretisation
   time_integrators
   osns_problems
   numerics_solvers
