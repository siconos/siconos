.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. index::
   single: Second Order Cone Linear Complementarity Problem (SOCLCP)

.. contents::

.. _soclcp_problem:

Second Order Cone Linear Complementarity Problem (SOCLCP)
*********************************************************

Problem statement
=================

Given

* a symmetric positive semidefinite matrix :math:`{M} \in {{\mathrm{I\!R}}}^{n \times n}`

* a vector :math:`{q} \in {{\mathrm{I\!R}}}^n`

* a vector of coefficients :math:`\mu \in{{\mathrm{I\!R}}}^{n_c}`

the second-order cone linear complementarity problem (SOCLCP), denoted by :math:`\mathrm{SOCLCP}(M,q,\mu)`, is to find two vectors :math:`u\in{{\mathrm{I\!R}}}^n` and :math:`r\in {{\mathrm{I\!R}}}^n` such that

.. math::

    \begin{eqnarray*} \begin{cases}
    u = M r + q \\
    \ C^\star_{\mu} \ni {u} \perp r \in C_{\mu}
    \end{cases} \end{eqnarray*}

where the set :math:`C^{\alpha,\star}_{\mu^\alpha}` is the dual of :math:`C^{\alpha}_{\mu^\alpha}`.

The set C is the second-order cone given by

.. math::

    \begin{eqnarray}
    C_{\mu} = \{ r \} = \prod_{\alpha =1}^{n_c} C^\alpha_{\mu}
    \end{eqnarray}

with

.. math::

    \begin{eqnarray}
    C^\alpha_{\mu} = \{ r \mid \|[r_1, \ldots, r_{n^\alpha}]\| \leq \mu^\alpha * r_0 \} \subset {\mathrm{I\!R}}^{n^\alpha}
    \end{eqnarray}

Implementation in numerics
==========================

Structure to define the problem: :cpp:class:`SecondOrderConeLinearComplementarityProblem`.

The generic driver for all SOCLCP problems is :cpp:func:`soclcp_driver()`.

List of solvers: :cpp:enum:`SOCLCP_SOLVER`

.. _soclcp_solvers:

SOCLCP available solvers
========================

Gauss-Seidel (:cpp:enumerator:`SICONOS_SOCLCP_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

PSOR (Gauss-Seidel with overrelaxation) solver.

driver: :cpp:func:`soclcp_nsgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
* iparam[SICONOS_IPARAM_ERROR_EVALUATION]: error computation method,

    * SICONOS_ERROR_FULL_EVALUATION: full error computation, including v (default)
    * SICONOS_ERROR_LIGHT_EVALUATION: light error computation with incremental values on r, and a check of the absolute error at the end
    * SICONOS_ERROR_LIGHT_EVALUATION_NO_UPDATE: light error computation, without updating v

* iparam[SICONOS_IPARAM_SOCLCP_NSGS_WITH_RELAXATION] = 0: if 1, use over-relaxation
* iparam[SICONOS_IPARAM_NSGS_SHUFFLE]: if 1, shuffle the contact indices in the loop
* dparam[SICONOS_DPARAM_TOL] = 1e-4;
* dparam[SICONOS_DPARAM_SOCLCP_NSGS_RELAXATION] = 1.: relaxation parameter value

internal solver: :cpp:enumerator:`SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration`.


VI, fixed-point (:cpp:enumerator:`SICONOS_SOCLCP_VI_FPP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

VI formulation and fixed-point projection.

driver: :cpp:func:`soclcp_VI_FixedPointProjection()`

parameters: same as :cpp:enumerator:`SICONOS_VI_FPP`, see :ref:`vi_solvers`.


VI, Extra-gradient (:cpp:enumerator:`SICONOS_SOCLCP_VI_EG`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

VI formulation and extra-gradient solver.

driver: :cpp:func:`soclcp_VI_ExtraGradient()`

parameters: same as :cpp:enumerator:`SICONOS_VI_EG`, see :ref:`vi_solvers`.

Projections
"""""""""""

Used as internal solvers of :cpp:enumerator:`SICONOS_SOCLCP_NSGS`.

ids: :cpp:enumerator:`SICONOS_SOCLCP_ProjectionOnConeWithLocalIteration`,
:cpp:enumerator:`SICONOS_SOCLCP_ProjectionOnCone`, :cpp:enumerator:`SICONOS_SOCLCP_ProjectionOnConeWithRegularization`.

drivers:

* :cpp:func:`soclcp_projectionOnCone_solve` for ProjectionOnCone and ProjectionOnConeWithRegularization,
* :cpp:func:`soclcp_projectionOnConeWithLocalIteration` for ProjectionOnConeWithLocalIteration.


parameters:

* iparam[SICONOS_IPARAM_SOCLCP_PROJECTION_CONE_INDEX] (set by soclcp_nsgs)
* dparam[SICONOS_DPARAM_SOCLCP_PROJECTION_RHO] = 0.
