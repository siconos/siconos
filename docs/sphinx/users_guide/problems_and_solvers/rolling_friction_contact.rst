.. index::
   single: Rolling friction-contact problems

.. contents::

.. _rfc_problem:

Rolling friction-contact problems
*********************************

Problem statement
=================

...

.. _rfc_error:

Implementation in numerics
==========================

Structure to define the problem: :class:`RollingFrictionContactProblem`.

Solvers list : :enum:`FRICTION_SOLVER`  (id contains ROLLING_FRICTION_3D)

The generic drivers for friction-contact problems are:

* :func:`rolling_fc3d_driver`

Error strategy
==============

To set internal solver tolerance (when it makes sense!) use :func:`rolling_fc3d_set_internalsolver_tolerance`.

Check details in :ref:`fc_error`



.. _rfc3d_solvers:

Rolling-Friction 3D available solvers
=====================================

NSGS (:enumerator:`SICONOS_ROLLING_FRICTION_3D_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Non-Smooth Gauss Seidel solver.

**Driver:** :func:`rolling_fc3d_nsgs`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000 : Maximum iteration number
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] : error computation method,
  
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL : Full error computation with velocity computation
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (DEFAULT): Light error computation with incremental values on reaction verification of absolute error at the end
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT : only light error computation (velocity not computed)
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE :  we adapt the frequency of the full erro evaluation.

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0,  error computation frequency

* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] : shuffle the contact indices in the loop
  
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE : no shuffle
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE : shuffle only at the beginning
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP : shuffle in each iteration

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] = 0 : seed for the random generator in shuffling  contacts

* iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] : filter local solution if the local error is greater than 1.0

  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE (default) the filter is not applied
  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE  the filter is applied

* iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] : method uses overrelaxation

  * SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (default) relaxation is not used,
  * SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE  relaxation is used with parameter dparam[8],

  
* dparam[SICONOS_DPARAM_TOL] = 1e-4, user tolerance on the loop
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0
* dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE]  the relaxation parameter omega

Default internal solver : :enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`.

Projection on cone (:enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone`, ...)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


.. csv-table:: Projection on cone solvers
   :header: "Solver id", "Driver"
   :widths: 15, 30

   ":enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone`",":func:`rolling_fc3d_projectionOnCone_solve`"
   ":enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`",":func:`rolling_fc3d_projectionOnConeWithLocalIteration_solve`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]
* iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION]

  * SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE (default for ProjectionOnConeWithLocalIteration)
  * SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE (default for ProjectionOnCone)

* dparam[SICONOS_DPARAM_TOL] =1e-12

out :

iparam[SICONOS_FRICTION_NUMBER_OF_CONTACTS] : number of contacts





