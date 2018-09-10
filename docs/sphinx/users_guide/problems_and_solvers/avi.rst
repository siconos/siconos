.. index:: single: Affine Variational Inequalities (AVI)
.. _avi_problem:

Affine Variational Inequalities (AVI)
=====================================

.. rubric:: Problem:

The Affine Variational Inequality (AVI) is defined by

Given :math:`q\in\mathbb{R}^n` , :math:`M\in\mathbb{R}^{n\times n}` and a set :math:`K\in\mathbb{R}^n` , find :math:`z\in\mathbb{R}^n` such that:

.. math::

    \begin{equation*}(Mz+q)^T(y -z) \geq 0,\quad \text{ for all } y \in K,\end{equation*}

or equivalently,

.. math::

    \begin{equation*}-(Mz + q) \in \mathcal{N}_K(z)\end{equation*}

where :math:`\mathcal{N}_K` is the normal cone to :math:`K` at :math:`z` .

The AVI is a special case of a Variational Inequality (VI), where the function :math:`F` is affine. For VI solvers, see :ref:`Variational Inequality (VI) <vi_problem>` .

From more details on theory and analysis of AVI (and VI in general), we refer to

Facchinei, Francisco; Pang, Jong-Shi (2003), *Finite Dimensional Variational Inequalities and Complementarity Problems* , Vol. 1 & 2, Springer Series in Operations Research, Berlin-Heidelberg-New York: Springer-Verlag.

.. _avi_problem_1aviSolversList:
.. rubric:: Available solvers:

The solvers and their parameters are described in :ref:`Affine Variational Inequalities Solvers <avi_problem_solvers>` .

Use the generic function AVI_driver() to call one the the specific solvers listed below:

* :func:`avi_caoferris()` , direct solver for AVI based on pivoting method principle for degenerate problem.
  
  Choice of pivot variable is performed via lexicographic ordering

(see also the functions/solvers list in ``AVI_Solvers.h`` and numbering in ``AVI_cst.h`` )


.. index:: single: Affine Variational Inequalities Solvers
.. _avi_problem_solvers:

.. rubric:: Affine Variational Inequalities Solvers:

Overview of the available solvers for AVI and their required parameters.

For each solver, the input argument are:

* an :class:`AffineVariationalInequalities`

* the unknown z

* the value of the function F(z)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions struct, which contains iparam and dparam

.. _avi_problem_solvers_1aviCaoFerris:
.. rubric:: The Cao-Ferris algorithm:

Direct solver for AVI based on pivoting method principle for (degenerated) problem.

function: :func:`avi_caoferris()`

parameters:

* iparam[0] (in): max. number of iterations

* iparam[1] (in,out): boolean, 1 if memory has been allocated, 0 otherwise

* iparam[1] (out): number of iterations processed

