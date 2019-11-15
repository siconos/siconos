.. _numerics_solvers:

###############################
Adding a new solver in numerics
###############################

* A solver must be associated to a problem formulation.

* Write the corresponding documentation in :

  docs/sphinx/users_guide/<formulation>.rst: describe solver and all its parameters

  
* Choose an id and complete the corresponding enum:

  * id must be: SICONOS_<formulation>_<solver_name> (all upper-case), e.g. SICONOS_LCP_ENUM.
    
  * in <formulation>_cst.h, append solver id to SICONOS_<formulation>_SOLVER enum.
    

* Write a driver function for the solver, named <formulation>_<solver_name> in <formulation>_<solver_name>.c

  .. code-block:: c

     void <formulation>_<solver_name>(<NonSmooth problem type>*problem, double * unknonw_1, ..., int * info, SolverOptions* options)
  
* If the solver needs some specific values for iparam and dparam, write <formulation>_<solver_name>_set_options function in <formulation>_<solver_name>.c file

  This function must set all solver specific values of iparam and dparam.
  e.g. :

  .. code-block:: 

     void lcp_lexicolemke_set_options(SolverOptions* options)
     {
     options->iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = SICONOS_LCP_PIVOT_LEMKE;
     }

  
* Use enum values for all indices in iparam and dparam (if required, create or append some new value in <formulation>_cst.h variables).
  
* Add a case in solver_options_create function (SolverOptions.c) to deal with the minimal setup for the solver:
  
  * select default values for tolerance and maxiter and add a line to
    call

    solver_options_initialize(options, solverId, iter_max, tol)

    and a line to call <formulation>_<solver_name>_set_options(options) if (and only if)
    some specific values are required.




.. code-block:: c

   // Use default setup :
   int id = SICONOS_LCP_ENUM
   SolverOptions * options = solver_options_create(id);
   // optional : set parameters
   options->iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 1;
   // ...
   // call driver
   lcp_enum(problem, z, w, info, options);
   // ...
   // clear memory
   solver_options_clear(&options);
   
