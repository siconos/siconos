.. _problems_and_solvers:

####################################################
Nonsmooth problem formulations and available solvers
####################################################

The core of a Siconos simulation consists in writing a formulation of a given problem
and choosing a solver to compute the unknowns.

The different problems and their features are described below, with the available solvers
for each formulation.

Write and solve a problem with Siconos
**************************************

Create a problem
================


.. _solver_options:

Create and describe a solver
============================

Solver parameters are handled by the :cpp:class:`SolverOptions` object. It defines which solver is used and its parameters (tolerance...), and may
hold internal work arrays or data.

The simplest way to create and use a solver is to select its id (see the list of available solvers for each problem, and their ids, in the pages below) and initialize the solver with this id.


**Kernel (high-level) interface:**

.. code-block:: c++

   // -- C++ API --
   // use the solver id as a parameter of the one-step nonsmooth problem constructor
   auto problem = std::make_shared<siconos::nonsmooth_formulations::LCP>(SICONOS_LCP_LEMKE);
   // get the options:
   auto options = problem->numericsSolverOptions();


.. code-block:: python

   # -- Python API --
   import siconos.nonsmooth_formulations as snf
   import siconos.numerics as sn
   lcp = snf.LCP(sn.solver_ids.SICONOS_LCP_LEMKE)

**Numerics (low-level) interface:**

.. code-block:: c

   // -- C/C++ API --
   int id = SICONOS_LCP_LEMKE;
   SolverOptions * options = solver_options_create(id);

.. code-block:: python

   # -- Python API --
   import siconos.numerics as sn
   options = sn.SolverOptions(sn.solver_ids.SICONOS_LCP_LEMKE)

In all cases, the id is the only required input. All the other parameters have default values.

To change these default values, set the content of iparam or dparam explicitly,
or use :cpp:func:`solver_options_update_internal` for internal solvers.

e.g.:

.. code-block:: c

   // -- C/C++ API --
   options->dparam[SICONOS_DPARAM_TOL] = 1e-12;

   // Set the first internal solver of options to SICONOS_FRICTION_3D_NSN_AC
   // and change the maximum number of iterations of the internal solver
   int internal_solver_number = 0;
   // Reset the internal solver to SICONOS_FRICTION_3D_NSN_AC with default values for its parameters ...
   solver_options_update_internal(options, internal_solver_number, SICONOS_FRICTION_3D_NSN_AC);
   // and modify the max number of iterations.
   options->internalSolvers[internal_solver_number].iparam[SICONOS_IPARAM_MAX_ITER] = 1000;

.. code-block:: python

   # -- Python API --
   options.dparam[sn.params.SICONOS_DPARAM_TOL] = 1e-12
   options.update_internal(0, sn.solver_ids.SICONOS_FRICTION_3D_NSN_AC);
   options.internalSolvers[0].iparam[sn.params.SICONOS_IPARAM_MAX_ITER] = 1000



A SolverOptions object contains:

* an id (int) that uniquely identifies the solver,

* iparam, an array of integer parameters,

* dparam, an array of real (double) parameters,

* internalSolvers: an array of SolverOptions describing the (optional) internal solvers.


dparam indices common to all solvers
------------------------------------

* dparam[SICONOS_DPARAM_TOL] (in): solver tolerance
* dparam[SICONOS_DPARAM_RESIDU] (out): computed error

iparam indices common to all solvers
------------------------------------

* iparam[SICONOS_IPARAM_MAX_ITER] (in): maximum number of iterations allowed
* iparam[SICONOS_IPARAM_PREALLOC] (in): keep the workspace across calls (?), 0 (false) by default.
* iparam[SICONOS_IPARAM_ITER_DONE] (out): number of iterations done


Solve a problem
===============

**Numerics (low-level) interface:**

.. code-block:: c

   // -- C/C++ API --

   // Create the solver
   int id = SICONOS_LCP_LEMKE;
   SolverOptions * options = solver_options_create(id);

   // Call the driver
   lcp_lexicolemke(problem, z, w, info, options);

   // Clear memory
   solver_options_delete(&options);


.. code-block:: python

   # -- Python API --
   import siconos.numerics as sn
   options = sn.SolverOptions(sn.solver_ids.SICONOS_LCP_LEMKE)
   # ...
   sn.lcp_lexicolemke(problem, z, w, info, options)


.. toctree::
   :maxdepth: 3

   vi
   qp
   convex_qp
   lcp
   avi
   mcp
   mlcp
   ncp
   relay
   soclcp
   friction_contact
   global_friction_contact
   rolling_friction_contact
   generic_mechanical
