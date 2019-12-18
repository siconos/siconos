.. index::
   single: Generic mechanical problems (GMP)

.. contents::

.. _gmp_problem:

Generic mechanical problems
***************************

Problem statement
=================

Implementation in numerics
==========================

Structure to define the problem: :class:`GenericMechanicalProblem`.

Solvers list  :enum:`GENERIC_MECHANICAL_SOLVER`

The driver for generic mechanical problem is :func:`gmp_driver`.

GMP available solvers
=====================

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_2D_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

**Driver:** :func:`gmp_driver`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
* iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_ISREDUCED]

  * SICONOS_GENERIC_MECHANICAL_GS_ON_ALLBLOCKS (:func:`gmp_gauss_seidel`)
  * SICONOS_GENERIC_MECHANICAL_SUBS_EQUALITIE  (:func:`gmp_reduced_solve`)
  * SICONOS_GENERIC_MECHANICAL_ASSEMBLE_EQUALITIES (:func:`gmp_reduced_equality_solve`)
  * SICONOS_GENERIC_MECHANICAL_MLCP_LIKE (:func:`gmp_as_mlcp`)

* iparam[SICONOS_GENERIC_MECHANICAL_IPARAM_WITH_LINESEARCH] = 0 (false)
* dparam[SICONOS_DPARAM_TOL] = 1e-4
* dparam[SICONOS_DPARAM_GMP_COEFF_LS] = 1.
  
There are 3 internal solvers :

* internalSolvers[0] = :enumerator:`SICONOS_LCP_LEMKE`
* internalSolvers[1] = :enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`
  with

  * internalSolvers[1]->iparam[SICONOS_IPARAM_MAX_ITER] = 100;
  * internalSolvers[1]->dparam[SICONOS_DPARAM_TOL] = 1e-12;
 
* internalSolvers[2] = :enumerator:`SICONOS_RELAY_LEMKE`

out :

* dparam[SICONOS_DPARAM_RESIDU]
* dparam[SICONOS_DPARAM_GMP_ERROR_LS]

* iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
* dparam[SICONOS_DPARAM_TOL] = 1e-4
