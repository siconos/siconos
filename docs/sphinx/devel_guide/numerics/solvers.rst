.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _numerics_solvers:

###############################
Adding a new solver in numerics
###############################

* A solver must be associated with a problem formulation.

* Write the corresponding documentation in:

  docs/sphinx/users_guide/<formulation>.rst: describe the solver and all its parameters.


* Choose an id and complete the corresponding enum:

  * the id must be SICONOS_<formulation>_<solver_name> (all upper-case), e.g. SICONOS_LCP_ENUM.

  * in <formulation>_cst.h, append the solver id to the SICONOS_<formulation>_SOLVER enum.


* Write a driver function for the solver, named <formulation>_<solver_name> in <formulation>_<solver_name>.c

  .. code-block:: c

     void <formulation>_<solver_name>(<NonSmooth problem type>*problem, double * unknown_1, ..., int * info, SolverOptions* options)

* If the solver needs specific values for iparam and dparam, write a <formulation>_<solver_name>_set_default function in the <formulation>_<solver_name>.c file.

  This function must set all the solver-specific values of iparam and dparam,
  e.g.:

  .. code-block::

     void lcp_lexicolemke_set_default(SolverOptions* options)
     {
     options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_LEMKE;
     }


* Use enum values for all indices in iparam and dparam (if required, add new values in <formulation>_cst.h).

* Add a case in the solver_options_create function (SolverOptions.c) for the minimal setup of the solver:

  * choose default values for the tolerance and maxiter and add a line
    calling

    solver_options_initialize(options, solverId, iter_max, tol)

    and a line calling <formulation>_<solver_name>_set_default(options) if (and only if)
    specific values are required.




.. code-block:: c

   // Use the default setup:
   int id = SICONOS_LCP_ENUM
   SolverOptions * options = solver_options_create(id);
   // optional: set parameters
   options->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 1;
   // ...
   // call driver
   lcp_enum(problem, z, w, info, options);
   // ...
   // clear memory
   solver_options_delete(&options);
