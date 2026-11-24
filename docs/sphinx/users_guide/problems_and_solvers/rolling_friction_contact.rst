.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

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

Structure to define the problem: :cpp:class:`RollingFrictionContactProblem`.

List of solvers: :cpp:enum:`FRICTION_SOLVER` (ids containing ROLLING_FRICTION_3D)

The generic drivers for rolling friction-contact problems are:

* :cpp:func:`rolling_fc3d_driver`

Error strategy
==============

To set the tolerance of the internal solver (when it makes sense), use :cpp:func:`rolling_fc3d_set_internalsolver_tolerance`.

See the details in :ref:`fc_error`.



.. _rfc3d_solvers:

Available 3D rolling friction solvers
=====================================

NSGS (:cpp:enumerator:`SICONOS_ROLLING_FRICTION_3D_NSGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Gauss-Seidel solver.

**Driver:** :cpp:func:`rolling_fc3d_nsgs`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000: maximum number of iterations
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION]: error computation method,

  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL: full error computation, including the velocity
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (DEFAULT): light error computation with incremental values on the reaction, and a check of the absolute error at the end
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT: light error computation only (velocity not computed)
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE: the frequency of the full error evaluation is adapted.

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0: error computation frequency

* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE]: shuffle the contact indices in the loop

  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE: no shuffle
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE: shuffle only at the beginning
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP: shuffle at each iteration

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] = 0: seed of the random generator used to shuffle the contacts

* iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION]: filter the local solution if the local error is greater than 1.0

  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE (default): the filter is not applied
  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE: the filter is applied

* iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION]: use over-relaxation

  * SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (default): relaxation is not used,
  * SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE: relaxation is used with parameter dparam[8],


* dparam[SICONOS_DPARAM_TOL] = 1e-4: user tolerance on the loop
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0
* dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE]: the relaxation parameter omega

Default internal solver: :cpp:enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`.

Projection on cone (:cpp:enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone`, ...)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. list-table:: Projection on cone solvers
   :header-rows: 1
   :widths: 15 30

   * - Solver id
     - Driver
   * - :cpp:enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnCone`
     - :cpp:func:`rolling_fc3d_projectionOnCone_solve`
   * - :cpp:enumerator:`SICONOS_ROLLING_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`
     - :cpp:func:`rolling_fc3d_projectionOnConeWithLocalIteration_solve`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]
* iparam[SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_IPARAM_USE_TRIVIAL_SOLUTION]

  * SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_TRUE (default for ProjectionOnConeWithLocalIteration)
  * SICONOS_FRICTION_3D_NSGS_LOCALSOLVER_USE_TRIVIAL_SOLUTION_FALSE (default for ProjectionOnCone)

* dparam[SICONOS_DPARAM_TOL] =1e-12

Output:

iparam[SICONOS_FRICTION_NUMBER_OF_CONTACTS]: number of contacts
