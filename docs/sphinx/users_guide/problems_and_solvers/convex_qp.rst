.. index:: single: ConvexQP
.. _doxid-convexqp_problem:

ConvexQP
========

.. _doxid-convexqp_problem_1convexQPintro:
.. rubric:: Problem:

Given

* an integer :math:`n` , the dimension of the ambient space,

* a SDP matrix :math:` M \in \mathrm{I\!R}^{n \times n}`

* a vector :math:` q \in \mathrm{I\!R}^n`

* a matrix :math:` A \in \mathrm{I\!R}^{m times n}` of constraints

* a vector :math:` b \in \mathrm{I\!R}^m`

* a convex set :math:` {C} \in {{\mathrm{I\!R}}}^m`

the convex QP problem is to find a vector :math:`z\in{{\mathrm{I\!R}}}^n` ,

.. math::

    \begin{equation*} \begin{array}{lcl} \min & & \frac{1}{2} z^T M z + z^T q \\ s.t & & A z + b \in C \\ \end{array} \end{equation*}

and is most simple example is when :math:` b= 0 A =I` and we obtain

.. math::

    \begin{equation*} \begin{array}{lcl} \min & & \frac{1}{2} z^T M z + Z^T q \\ s.t & & z \in C \\ \end{array} \end{equation*}

Most of the solver returns

* the solution vector :math:` z \in \mathrm{I\!R}^n`

* the vector :math:` u \in \mathrm{I\!R}^m`

* the multiplier :math:` \xi \in \mathrm{I\!R}^m` such that :math:` - \xi \in \partial \Psi_C(u) `

* the vector :math:` w \in \mathrm{I\!R}^n` such that :math:` w =A^T \xi `

In the most simple case, we return

* the solution vector :math:` z = u \in \mathrm{I\!R}^n`

* the vector :math:` w =\xi \in \mathrm{I\!R}^m`


.. index:: single: ConvexQP as VI
.. _doxid-cqp_problem_v_i:


.. rubric:: ConvexQP Solvers:

This page gives an overview of the available solvers for Convex QP problems and their required parameters.

For each solver, the input argument are:

* a ConvexQP problem

* the unknowns (x,fx)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

