.. _problems_and_solvers:

#####################################################
Nonsmooth problems formulations and available solvers
#####################################################

The core of a Siconos simulation consists in writing a formulation for a given problem
and choose a solver to determine the unknowns.

The different problems and their features are described below, with the available solvers
for each formulation.

Write and solve a problem with Siconos
**************************************

Create a problem
================

Create and describe a solver
============================

Solver parameters are handled by the object SolverOptions. The simplest way to create and use a solver
is to select the corresponding id (check the list of avalaible solver for each problem and the corresponding numbers in the pages below) and initialize the solver with this id.


**Kernel (high-level) interface:**

.. code-block:: c++

   // -- C++ API --
   // use solver id as a parameter for the one-step nonsmooth problem constructor
   SP::OneStepNSProblem problem(new LCP(SICONOS_LCP_LEMKE));
   // get options :
   SP::SolverOptions options = problem->numericsSolverOptions() 


.. code-block:: python

   // -- Python API --
   import siconos.kernel as sk
   lcp = sk.LCP(sk.SICONOS_LCP_LEMKE)
   
**Numerics (low-level) interface:**

.. code-block:: c

   // -- C/C++ API --
   int id = SICONOS_LCP_LEMKE;
   SolverOptions * options =(SolverOptions *)malloc(sizeof(SolverOptions));
   solver_options_create(options, id);

.. code-block:: python

   // -- Python API --
   import siconos.numerics as sn
   options = sn.SolverOptions(sn.SICONOS_LCP_LEMKE)

In any case, the id is the only required input. All the other parameters have default values.

* an id (int) that uniquely identifies the solver,
  
* iparam, an array used to save integer type parameters,

* dparam, an array used to save real (double) type parameters,

* internalSolvers : array of SolverOptions used to describe (optional) internal solvers


dparam indices common to all solvers
------------------------------------

* dparam[SICONOS_DPARAM_TOL] (in): solver tolerance
* dparam[SICONOS_DPARAM_RESIDU] (out): computed error
  
iparam indices common to all solvers
------------------------------------

* iparam[SICONOS_IPARAM_MAX_ITER] (in): maximum number of iterations allowed
* iparam[SICONOS_IPARAM_PREALLOC] (in): keep work space across calls (?), 0 (false) by default.
* iparam[SICONOS_IPARAM_ITER_DONE] (out): number of iterations done

  
Solve a problem
===============

**Numerics (low-level) interface:**

.. code-block:: c

   // -- C/C++ API --

   // Create the solver
   int id = SICONOS_LCP_LEMKE;
   SolverOptions * options =(SolverOptions *)malloc(sizeof(SolverOptions));
   solver_options_create(options, id);

   // Call the driver
   lcp_lexicolemke(problem, z, w, info, options);

   // Clear memory
   solver_options_clear(options);
   free(options);
   
   
.. code-block:: python

   // -- Python API --
   import siconos.numerics as sn
   options = sn.SolverOptions(sn.SICONOS_LCP_LEMKE)
   // ...
   sn.lcp_lexicolemke(problem, z, w, info, options)
   

.. toctree::
   :maxdepth: 3

   vi
   qp
   convex_qp
   lcp
   avi
   friction_contact
   global_friction_contact
   mcp
   mlcp
   ncp
   relay
   soclcp
   
