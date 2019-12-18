.. index::
   single: Relay or box-constrained AVI problems
   
.. contents::

.. _lcp_problem:

Relay or box-constrained AVI problems
*************************************

Problem statement
=================

Find :math:`(z,w)` such that:

.. math::

    \begin{equation*}
    \left\lbrace \begin{array}{l}
    w = M z + q\\
    -w \in \mathcal{N}_{K}(z)\\
    \end{array}, \right.
    \end{equation*}

where M is an ( :math:`n \times n` )-matrix, q, z and w are n-dimensional vectors, K is the box defined by :math:`K=\{x\in\mathbb{R}^n\mid lb_i\leq x_i\leq ub_i,i=1,...,n\}` and :math:`\mathcal{N}_K(z)` is the normal cone to :math:`K` at :math:`z` .

Implementation in numerics
==========================

Structure to define the problem: :class:`RelayProblem`.

The generic driver for all Relay problems is :func:`relay_driver()`.

Solvers list  :enum:`RELAY_SOLVER`

.. _relay_solvers:

Relay available solvers
=======================

Direct solvers
--------------

Enumerative solver (:enumerator:`SICONOS_RELAY_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

The relay problem is reformulated as a LCP and solved with enum method, see :ref:`lcp_solvers`.

driver: :func:`relay_enum()`

parameters: same as :enumerator:`SICONOS_LCP_ENUM`.

Lemke solver (:enumerator:`SICONOS_RELAY_LEMKE`)
""""""""""""""""""""""""""""""""""""""""""""""""

The relay problem is reformulated as a LCP and solved with Lemke's method, see :ref:`lcp_solvers`.

driver: :func:`relay_lexicolemke()`

parameters: same as :enumerator:`SICONOS_LCP_ENUM`.


PATH (:enumerator:`SICONOS_RELAY_PATH`)
"""""""""""""""""""""""""""""""""""""""

The relay problem is reformulated as a LCP and solved with the PATH solver

*Works only if Siconos has been built with path support (if PathFerris or PathVI has been found, see :ref:`siconos_install_guide`)**

driver: :func:`relay_path()`

parameters : none.

AVI reformulation
-----------------

AVI, Cao/Ferris solver (:enumerator:`SICONOS_RELAY_AVI_CAOFERRIS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Based on an algorithm by Cao and Ferris for AVI with a polytopic set :math:`K` .

driver:  :func:`relay_avi_caoferris()`

parameters: same as :enumerator:`SICONOS_AVI_CAOFERRIS`, see :ref:`avi_solvers`.


There also exists a test version :enumerator:`SICONOS_RELAY_AVI_CAOFERRIS_TEST` with 

driver:  :func:`relay_avi_caoferris_test()`

parameters: same as :enumerator:`SICONOS_AVI_CAOFERRIS`, see :ref:`avi_solvers`.


Iterative solvers
-----------------

Projected Gauss-Seidel (:enumerator:`SICONOS_RELAY_PGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :func:`relay_pgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
