.. _time_discretisation:

Time discretisation
===================

The discretisation scheme is characterized by a vector of size nSteps+1, denoted tk, that handles the values at each time-step, nSteps being the number of time steps. At the time, only a constant time-step (denoted h) time discretisation is available.

A TimeDiscretisation must be linked to one and only one Model (required parameter of any constructor). This Model provides the initial time value and possibly the final time value, which is an optional parameter. Thus depending on the constructor you use, you may need to give the time-step, the final time ... \n
Just take care to avoid redundant or conflicting information. See the constructors list in the TimeDiscretisation class documentation for full details. 

Example::

  double h = ...;
  unsigned int nSteps = ...;
  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(h, nSteps,m);
  // tk and final time will be automatically computed. 
  // 
  // Then a simulation is created and associated to m through td. 
  auto s = std::make_shared<siconos::simulation::TimeStepping>(td);
