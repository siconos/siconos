.. _simulation_nsds:

Nonsmooth dynamical systems simulation
======================================

.. toctree::
   :maxdepth: 3

	      
A "simulation" is the description of how a nonsmooth dynamical system will be computed in Siconos.
At the time, two strategies are available:

* time-stepping schemes
* event-driven schemes
  
Key points of the simulation are :
* the time discretisation
* the time integration of dynamical systems
* the formulation and resolution of one or more nonsmooth problems
  

Time discretisation
-------------------

The discretisation scheme is characterized by a vector of size nSteps+1, denoted tk, that handles the values at each time-step, nSteps being the number of time steps. At the time, only a constant time-step (denoted h) time discretisation is available.

A TimeDiscretisation must be linked to one and only one Model (required parameter of any constructor). This Model provides the initial time value and possibly the final time value, which is an optional parameter. Thus depending on the constructor you use, you may need to give the time-step, the final time ... \n
Just take care to avoid redundant or conflicting information. See the constructors list in the TimeDiscretisation class documentation for full details. 

Example::

  SP::Model m(new Model(t0)); // only initial time is provided
  double h = ...;
  unsigned int nSteps = ...;
  SP::TimeDiscretisation td(new TimeDiscretisation(h, nSteps,m));
  // tk and final time will be automatically computed. 
  // 
  // Then a simulation is created and associated to m through td. 
  SP::Simulation s(new TimeStepping(td));


Time integration of the dynamics
--------------------------------

Nonsmooth problems formulation and solve
----------------------------------------
