.. _simulation_nsds:

Nonsmooth dynamical systems simulation
======================================

Once the model (i.e. the nonsmooth problem) has been properly defined and described as explained in :ref:`modeling_nsds`,
the simulation strategy of this system must be defined, which mainly consists in:

* describe a time discretisation
* describe how dynamical systems will be integrated
* choose a formalisation and a solver for the nonsmooth problem

And above all, a global strategy has to be chosen, among :

* event-capturing (time-stepping) schemes
* event-driven schemes

All these steps are described in the pages below.

.. toctree::
   :maxdepth: 2
	      
   simulation_strategies
   time_discretisation
   time_integrators
   osns_problems
   numerics_solvers
