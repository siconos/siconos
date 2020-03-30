.. index::
   single: Mixed (Non Linear) Complementarity problem (MCP)
   
.. contents::

.. _mcp_problem:

Mixed (Non Linear) Complementarity problem (MCP)
************************************************

Problem statement
=================

Given a sufficiently smooth function :math:`{F} \colon {{\mathrm{I\!R}}}^{n+m} \to {{\mathrm{I\!R}}}^{n+m}` , the Mixed Complementarity problem (MCP) is to find two vectors :math:`(z,w \in {{\mathrm{I\!R}}}^{n+m})` such that:

.. math::

    \begin{align*} w &= \begin{pmatrix}w_e\\w_i\end{pmatrix} = F(z) \\ w_e &=0 \\ 0 &\le w_i \perp z_i \ge 0, \end{align*}

where "i" (resp. "e") stands for inequalities (resp. equalities). The vector :math:`z` is splitted like :math:`w` :

.. math::

    \begin{equation*}z =\begin{pmatrix}z_e\\z_i\end{pmatrix}.\end{equation*}

The vectors :math:`z_i,w_i` are of size ``sizeEqualities`` , the vectors :math:`z_e,w_e` are of size ``sizeInequalities`` and :math:`F` is a non linear function that must be user-defined.

A Mixed Complementarity problem (MCP) is a :ref:`ncp_problem` "augmented" with equality constraints.

Implementation in numerics
==========================

Structure to define the problem: :class:`MixedLinearComplementarityProblem`.

The generic driver for all MCP is :func:`mcp_driver()`.

Solvers list  :enum:`MCP_SOLVER`

.. _mcp_error:

.. _mcp_solvers:

MCP available solvers
=====================

Newton, Fisher-Burmeister (:enumerator:`SICONOS_MCP_NEWTON_FB_FBLSA`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Solver based on Fischer-Burmeister reformulation and line search (VFBLSA in Facchinei--Pang 2003 p. 834)

driver: :func:`mcp_newton_FB_FBLSA()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-10
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;

Newton (min), Fisher-Burmeister (:enumerator:`SICONOS_MCP_NEWTON_MIN_FB_FBLSA`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Solver based on Fischer-Burmeister reformulation and line search.
The descent direction is found using a min reformulation (minFBLSA in Facchinei--Pang 2003 p. 855)

driver: :func:`mcp_newton_min_FB_FBLSA()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-10
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;
Newton, Fisher-Burmeister (:enumerator:`SICONOS_MCP_OLD_FB`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Solver based on Fischer-Burmeister reformulation, old (outdated) version, for the records.

driver: :func:`mcp_old_FischerBurmeister()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-7
