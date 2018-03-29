.. index:: single: Second Order Cone Linear Complementarity Problem
.. _doxid-soclcp_problem:

Second Order Cone Linear Complementarity Problem
================================================

.. _doxid-soclcp_problem_1soclcpIntro:
.. rubric:: Problem statement:

Given

* a symmetric positive semi definite matrix :math:`{M} \in {{\mathrm{I\!R}}}^{n \times n} `

* a vector :math:` {q} \in {{\mathrm{I\!R}}}^n`

* a vector of coefficients :math:`\mu \in{{\mathrm{I\!R}}}^{n_c}`

the second order cone linear complementarity problem (SOCLCP) is to find two vectors :math:`u\in{{\mathrm{I\!R}}}^n` , and :math:`r\in {{\mathrm{I\!R}}}^n` , denoted by :math:`\mathrm{SOCCLP}(M,q,\mu)` such that

.. math::

    \begin{eqnarray*} \begin{cases} u = M r + q \\ \ C^\star_{\mu} \ni {u} \perp r \in C_{\mu} \end{cases} \end{eqnarray*}

and the set :math:`C^{\alpha,\star}_{\mu^\alpha}` is its dual.

The set C is the second order cone given by

.. math::

    \begin{eqnarray} C_{\mu} = \{ r \} = \prod_{\alpha =1}^{n_c} C^\alpha_{\mu} \end{eqnarray}

with

.. math::

    \begin{eqnarray} C^\alpha_{\mu} = \{ r \mid \|[r_1, \ldots, r_{n^\alpha}]\| \leq \mu^\alpha * r_0 \} \subset {\mathrm{I\!R}}^{n^\alpha} \end{eqnarray}

The problem is stored and given to the solver in numerics thanks to the C structure :ref:`SecondOrderConeLinearComplementarityProblem <doxid-struct_second_order_cone_linear_complementarity_problem>` .

.. _doxid-soclcp_problem_1SOCLCPSolversList:
.. rubric:: Available solvers for SOCCLP:

see ``SOCLCP_cst.h`` for solver ids.

Use the generic function :ref:`soclcp_driver() <doxid-_non_smooth_drivers_8h_1a6fc996a07f24a72c0dbe0dc4b85733b9>` to call one the the specific solvers listed below:

* :ref:`soclcp_nsgs() <doxid-_s_o_c_l_c_p___solvers_8h_1ac8c7dcc8d1f1ee0a12ae16384561029e>` : PSOR (Gauss-Seidel with overrelaxation) solver. SolverId : SICONOS_SOCLCP_NSGS ,

* soclcp_VI_FixedPointProjection() : VI formulation and fixed point projection. SolverId : SICONOS_SOCLCP_VI_FPP ,

* :ref:`soclcp_VI_ExtraGradient() <doxid-_s_o_c_l_c_p___solvers_8h_1a60adc9ddff54149b8188610d4df499a1>` : VI formulation and extra-gradient solver. SolverId : SICONOS_SOCLCP_VI_EG ,

See the related functions/solvers list in ``SOCLCP_Solvers.h`` .

.. index:: single: Second Order Cone Linear Complementarity Problem (SOCLCP) solvers
.. _doxid-_second_order_cone_linear_complementarity_problem_solvers:

.. rubric:: Second Order Cone Linear Complementarity Problem (SOCLCP) solvers:


This page gives an overview of the available solvers for Second Order Cone Linear Complementarity Problem (SOCLCP) and their required parameters.

This page gives an overview of the available solvers for Second Order Cone Linear Complementarity Problem (SOCLCP) and their required parameters.

For each solver, the input argument are:

* a :ref:`SecondOrderConeLinearComplementarityProblem <doxid-struct_second_order_cone_linear_complementarity_problem>`

* the unknowns (r,v)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

.. _doxid-_second_order_cone_linear_complementarity_problem_solvers_1soclcp:
.. rubric:: nsgs Non-Smooth Gauss Seidel Solver:

function: secondOrderConeLinearComplementarity_nsgs(problem, r , v , &info , options); parameters:

