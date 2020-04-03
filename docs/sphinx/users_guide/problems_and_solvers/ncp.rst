.. index::
   single: Nonlinear Complementarity Problems (NCP)
   
.. contents::

.. _ncp_problem:

Nonlinear Complementarity Problems (NCP)
****************************************

Problem statement
=================

Given a sufficiently smooth function :math:`{F}\colon {{\mathrm{I\!R}}}^{n} \to {{\mathrm{I\!R}}}^{n}` The Nonlinear Complementarity Problem (NCP) is to find two vectors :math:`(z,w \in {{\mathrm{I\!R}}}^{n})` such that:

.. math::

    \begin{align*} w &= F(z) \\ 0 &\le w \perp z \ge 0 \end{align*}


Implementation in numerics
==========================

Structure to define the problem: :class:`NonlinearComplementarityProblem`.

The generic driver for all NCP is :func:`ncp_driver()`.

Solvers list  :enum:`NCP_SOLVER`
 
.. _ncp_solvers:

NCP available solvers
=====================

Newton, Fisher-Burmeister (:enumerator:`SICONOS_NCP_NEWTON_FB_FBLSA`)
---------------------------------------------------------------------

NCP Solver using the FB merit function and a Newton-based method with line-search

driver: :func:`ncp_newton_FBLSA`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;

Newton, min merit function (:enumerator:`SICONOS_NCP_NEWTON_MIN_FBLSA`)
-----------------------------------------------------------------------

NCP Solver using the min merit function (+ the FB as backup) and a Newton-based method with line-search

driver: :func:`ncp_newton_minFBLSA`

parameters: same as :enumerator:`SICONOS_NCP_NEWTON_FB_FBLSA`

Path search algorithm (:enumerator:`SICONOS_NCP_PATHSEARCH`)
------------------------------------------------------------

NCP Solver using a path search algorithm, following the work of D. Ralph.
M. Ferris, and many other collaborators of the latter.

driver: :func:`ncp_newton_pathsearch`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 100
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = NM_LS_MEAN;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 10;
* iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = 5;
* iparam[SICONOS_IPARAM_NMS_WATCHDOG_TYPE] = LINESEARCH;
* iparam[SICONOS_IPARAM_NMS_PROJECTED_GRADIENT_TYPE] = ARCSEARCH;
* iparam[SICONOS_IPARAM_NMS_N_MAX] = 10;
* iparam[SICONOS_IPARAM_PREALLOC] = 0;
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_NMS_DELTA] = 20;
* dparam[SICONOS_DPARAM_NMS_DELTA_VAR] = .8;
* dparam[SICONOS_DPARAM_NMS_SIGMA] = .01;
* dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_WATCHDOG] = 1e-12;
* dparam[SICONOS_DPARAM_NMS_ALPHA_MIN_PGRAD] = 1e-12;
* dparam[SICONOS_DPARAM_NMS_MERIT_INCR] = 1.1

PATH (Ferris) solver (:enumerator:`SICONOS_NCP_PATH`)
-----------------------------------------------------

driver: :func:`ncp_path`

parameters:

* dparam[SICONOS_DPARAM_TOL] = 1e-12
* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
