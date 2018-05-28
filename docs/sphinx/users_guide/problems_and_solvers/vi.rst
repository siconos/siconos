.. index:: single: Variational Inequality (VI)
.. _doxid-vi_problem:

Variational Inequality (VI)
===========================

.. _doxid-vi_problem_1viIntro:
.. rubric:: Problem statement:

Given

* an integer :math:`n` , the dimension of the ambient space,

* a mapping :math:` F\colon \mathrm{I\!R}^n \rightarrow \mathrm{I\!R}^n`

* a set :math:` {X} \in {{\mathrm{I\!R}}}^n`

the variational inequality problem is to find a vector :math:`z\in{{\mathrm{I\!R}}}^n` ,

.. math::

    \begin{equation*} F(z)^T(y-z) \geq 0,\quad \text{ for all } y \in X \end{equation*}

or equivalently,

.. math::

    \begin{equation*} - F(z) \in \mathcal{N}_X(z) \end{equation*}

where :math:`\mathcal{N}_X` is the normal cone to :math:`X` at :math:`z` .

Reference

Facchinei, Francisco; Pang, Jong-Shi (2003), *Finite Dimensional Variational Inequalities and Complementarity Problems* , Vol. 1 & 2, Springer Series in Operations Research, Berlin-Heidelberg-New York: Springer-Verlag.

The problem is stored and given to the solver in Siconos/Numerics thanks to a C structure VariationalProblem.

.. _doxid-vi_problem_1viSolversList:
.. rubric:: Available solvers for Variational Inequality:

Use the generic function :func:`variationalInequality_driver()` to call one the the specific solvers listed below:

* :func:`variationalInequality_ExtraGradient()` : Extra gradient solver. SolverId: SICONOS_VI_EG = 1000.

* :func:`variationalInequality_FixedPointProjection()` : Fixed-point solver. SolverId: SICONOS_VI_EG = 1001.

* :func:`variationalInequality_HyperplaneProjection()` : Hyperplane Projection based Solver. SolverId: SICONOS_VI_HP_STR = 1002.

* :func:`variationalInequality_box_newton_QiLSA()` : Solver using the merit function proposed by Qi for box-constrained VI. SolverId: SICONOS_VI_BOX_QI_STR = 1003

(see the functions/solvers list in ``VariationalInequality_Solvers.h`` )

.. index:: single: VI problems Solvers
.. _doxid-_v_i_solvers:

.. rubric:: VI problems Solvers:

This page gives an overview of the available solvers for Variational Inequality problems and their required parameters.

This page gives an overview of the available solvers for Variational Inequality problems and their required parameters.

For each solver, the input argument are:

* a VariationalInequality

* the unknowns (x,fx)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

