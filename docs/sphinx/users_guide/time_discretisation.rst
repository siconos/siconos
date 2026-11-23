.. _time_discretisation:

Time discretisation
===================

The discretisation scheme is characterized by a vector of size nSteps+1, denoted tk, that holds the time value at each time step, nSteps being the number of time steps. At the moment, only a constant time step (denoted h) is available.

A TimeDiscretisation must be linked to one and only one Model (a required parameter of every constructor). This Model provides the initial time and possibly the final time, which is optional. Depending on the constructor you use, you may therefore need to give the time step, the final time, etc.
Avoid giving redundant or conflicting information. See the list of constructors in the TimeDiscretisation class documentation for details.

Example::

  double h = ...;
  unsigned int nSteps = ...;
  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(h, nSteps,m);
  // tk and final time will be automatically computed.
  //
  // Then a simulation is created and associated with m through td.
  auto s = std::make_shared<siconos::simulation::TimeStepping>(td);
