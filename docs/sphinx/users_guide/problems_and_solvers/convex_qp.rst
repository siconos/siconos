.. index::
   single: Convex Quadratic Programming problems (ConvexQP)
   
.. contents::

.. _convexqp_problem:

Convex Quadratic Programming (ConvexQP) problems
************************************************

Problem statement
=================


Given

* an integer :math:`n` , the dimension of the ambient space,

* a SDP matrix :math:`M \in \mathrm{I\!R}^{n \times n}`

* a vector :math:`q \in \mathrm{I\!R}^n`

* a matrix :math:`A \in \mathrm{I\!R}^{m times n}` of constraints

* a vector :math:`b \in \mathrm{I\!R}^m`

* a convex set :math:`{C} \in {{\mathrm{I\!R}}}^m`

the convex QP problem is to find a vector :math:`z\in{{\mathrm{I\!R}}}^n` ,

.. math::

    \begin{equation*} \begin{array}{lcl} \min & & \frac{1}{2} z^T M z + z^T q \\ s.t & & A z + b \in C \\ \end{array} \end{equation*}

and is most simple example is when :math:`b= 0 A =I` and we obtain

.. math::

    \begin{equation*} \begin{array}{lcl}
    \min & & \frac{1}{2} z^T M z + Z^T q \\
    s.t & & z \in C \\
    \end{array}\end{equation*}

Most of the solver returns

* the solution vector :math:`z \in \mathrm{I\!R}^n`

* the vector :math:`u \in \mathrm{I\!R}^m`

* the multiplier :math:`\xi \in \mathrm{I\!R}^m` such that :math:`-\xi \in \partial \Psi_C(u)`

* the vector :math:`w \in \mathrm{I\!R}^n` such that :math:`w =A^T \xi`

In the most simple case, we return

* the solution vector :math:`z = u \in \mathrm{I\!R}^n`

* the vector :math:`w =\xi \in \mathrm{I\!R}^m`

Implementation in numerics
==========================

Structure to define the problem: :struct:`ConvexQP`.

No generic driver.

solvers list  :enum:`CONVEXQP_SOLVER`

.. _convex_qp_solvers:

Available solvers
=================

convex QP, projected gradient (:enumerator:`SICONOS_CONVEXQP_PG`)
----------------------------------------------------------------

driver :func:`convexQP_ProjectedGradient()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000
* iparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MAX_ITER] =20
* dparam[SICONOS_CONVEXQP_PGOC_RHO] = -1.e-3 /* rho is variable by default */
* dparam[SICONOS_CONVEXQP_PGOC_RHOMIN] = 1e-9
* dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_MU] =0.9
* dparam[SICONOS_CONVEXQP_PGOC_LINESEARCH_TAU] = 2.0/3.0
  
* dparam[SICONOS_DPARAM_TOL] = 1e-6

convex QP, VI solvers`(:enumerator:`SICONOS_CONVEXQP_VI_FPP` and :enumerator:`SICONOS_CONVEXQP_VI_EG`)
------------------------------------------------------------------------------------------------------

Rewrite QP as Variational Inequality problem.

ids: SICONOS_CONVEXQP_VI_FPP (fixed-point projection) and SICONOS_CONVEXQP_VI_EG (extra-gradient)

driver :func:`convexQP_VI_solver()`

parameters:

same as SICONOS_VI_FPP and SICONOS_VI_EG (see :ref:`vi_solvers`).

convex QP, ADMM (:enumerator:`SICONOS_CONVEXQP_ADMM`)
-----------------------------------------------------

ADMM, alternating direction method of multipliers

driver :func:`convexQP_ADMM()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000
* iparam[SICONOS_CONVEXQP_ADMM_IPARAM_ACCELERATION] = SICONOS_CONVEXQP_ADMM_ACCELERATION_AND_RESTART
* dparam[SICONOS_DPARAM_TOL] = 1e-6
  
* dparam[SICONOS_CONVEXQP_ADMM_RHO] = 1.0
* dparam[SICONOS_CONVEXQP_ADMM_RESTART_ETA] = 0.999

