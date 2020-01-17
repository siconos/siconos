.. index::
   single: Affine Variational Inequalities (AVI)
   
.. contents::

.. _avi_problem:

Affine Variational Inequalities (AVI)
*************************************

Problem statement
=================

The Affine Variational Inequality (AVI) is defined by

Given :math:`q\in\mathbb{R}^n` , :math:`M\in\mathbb{R}^{n\times n}` and a set :math:`K\in\mathbb{R}^n` , find :math:`z\in\mathbb{R}^n` such that:

.. math::

    \begin{equation*}(Mz+q)^T(y -z) \geq 0,\quad \text{ for all } y \in K,\end{equation*}

or equivalently,

.. math::

    \begin{equation*}-(Mz + q) \in \mathcal{N}_K(z)\end{equation*}

where :math:`\mathcal{N}_K` is the normal cone to :math:`K` at :math:`z` .

The AVI is a special case of a Variational Inequality (VI), where the function :math:`F` is affine. For VI solvers, see :ref:`vi_problem`.

From more details on theory and analysis of AVI (and VI in general), we refer to :cite:`Facchinei.2003`

Implementation in numerics
==========================

Structure to define the problem: :class:`AffineVariationalInequalities`.

The generic driver for all VI is :func:`avi_driver()`.

solvers list  :enum:`AVI_SOLVER`

.. _avi_solvers:

AVI Available solvers
=====================

Cao Ferris (:enumerator:`SICONOS_AVI_CAOFERRIS`)
------------------------------------------------

Direct solver for AVI based on pivoting method principle for degenerate problem. Choice of pivot variable is performed via lexicographic ordering

driver: :func:`avi_caoferris()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Pathvi (:enumerator:`SICONOS_AVI_PATHAVI`)
------------------------------------------

Direct solver for VI based on pivoting method principle for degenerate problem.

Ref: "A structure-preserving Pivotal Method for Affine Variational Inequalities" Y. Kim, O. Huber, M.C. Ferris, Math Prog B (2017).

driver: :func:`avi_pathvi()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
