.. index:: single: Linear Complementarity problems (LCP)
.. _doxid-_l_c_problem:

Linear Complementarity problems (LCP)
=====================================

.. _doxid-_l_c_problem_1lcpIntro:
.. rubric:: The problem:

The Linear Complementarity problem (LCP) is defined by

Find :math:`(z,w)` such that:

.. math::

    \begin{equation*} \begin{cases} M \ z + q = w \\ 0 \le w \perp z \ge 0 \end{cases}, \end{equation*}

where :math:` w, z, q` are vectors of size :math:`n` and :math:` M ` is a :math:`n\times n` matrix.

The notation :math:`x \perp y` means that :math:`x^Ty =0` . Inequalities involving vectors are understood to hold component-wise.

From more details on theory and analysis of LCP, we refer to

R.W. Cottle, J.S. Pang, and R.E. Stone. *The Linear Complementarity Problem.* Academic Press, Inc., Boston, MA, 1992.

The problem is stored and given to the solver in numerics thanks to the C structure :ref:`LinearComplementarityProblem <doxid-struct_linear_complementarity_problem>` .

.. _doxid-_l_c_problem_1lcpSolversList:
.. rubric:: Available solvers:

Use the generic functions :ref:`lcp_driver_DenseMatrix() <doxid-_l_c_p___solvers_8h_1a077f9a6851d0cbb0574e4fa7012804ab>` to call one the the specific solvers listed below:

.. _doxid-_l_c_problem_1lcpDirectSolvers:
.. rubric:: Direct solvers:

* :ref:`lcp_lexicolemke() <doxid-_l_c_p___solvers_8h_1a557a22faaecdcab07199133a96384dfc>` , direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

* :ref:`lcp_pivot() <doxid-_l_c_p___solvers_8h_1ad44408b578fc030dcccbb4bdd0ebd188>` , generic solver for pivot-based methods: Bard, Murty and Lemke rules are implemented.

* :ref:`lcp_enum() <doxid-_l_c_p___solvers_8h_1aae09a41219e84f0bf6d4a8f2b8d31bac>` , enumerative solver (brute-force method which tries every possible solution).

* :ref:`lcp_path() <doxid-_l_c_p___solvers_8h_1ae1e4d7c9521eade11bbaf99f4a12dee6>` , interface to the PATH solver

.. _doxid-_l_c_problem_1lcpIterativeSolvers:
.. rubric:: Iterative solvers:

* :ref:`lcp_cpg() <doxid-_l_c_p___solvers_8h_1a38dcc1a3c83b4246f76ef407b4b170bb>` , CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.

* :ref:`lcp_pgs() <doxid-_l_c_p___solvers_8h_1a347e902918529115bfed0a96292e126c>` , PGS is a basic Projected Gauss-Seidel solver for LCP.

* :ref:`lcp_rpgs() <doxid-_l_c_p___solvers_8h_1acd882d2802211d19d8dc9b324cda02c7>` , Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms.

* :ref:`lcp_psor() <doxid-_l_c_p___solvers_8h_1aa26a3cd3c19ad5869fc5f93c31f3c6ec>` , Projected Succesive over relaxation solver for LCP. See Cottle, Pang Stone Chap 5.

* :ref:`lcp_latin() <doxid-_l_c_p___solvers_8h_1a6a4bfa635f9d5a1865cd8b2bf29506dc>` , (LArge Time INcrements) is a basic latin solver for LCP.

* :ref:`lcp_latin_w() <doxid-_l_c_p___solvers_8h_1a098e3401a6c43fe21365c5946ee80553>` , (LArge Time INcrements) is a basic latin solver with relaxation for LCP.

* :ref:`lcp_nsgs_SBM() <doxid-_l_c_p___solvers_8h_1afb791fbc7113d35c90686705df078a2a>` , Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

.. _doxid-_l_c_problem_1lcpEquationBasedSolvers:
.. rubric:: Equation-based solvers:

* :ref:`lcp_newton_min() <doxid-_l_c_p___solvers_8h_1a759c437f42dfae616b69cf80eb595fe5>` , nonsmooth Newton method based on the min formulation of the LCP.

* :ref:`lcp_newton_FB() <doxid-_l_c_p___solvers_8h_1ab44e0b227d2e6f1f96a906f2f023c238>` , uses a nonsmooth newton method based on the Fischer-Bursmeister NCP function.

* :ref:`lcp_newton_minFB() <doxid-_l_c_p___solvers_8h_1a35fe4f067b94917f5b25d210b9274db3>` , nonsmooth Newton method combining the min and FB functions.

.. _doxid-_l_c_problem_1lcpReformulation:
.. rubric:: QP-reformulation:

* :ref:`lcp_qp() <doxid-_l_c_p___solvers_8h_1ac09892b3b83ec2a500231790368041ed>` , quadratic programm formulation

* :ref:`lcp_nsqp() <doxid-_l_c_p___solvers_8h_1af8140adefc613da343d58d751e32c54e>` , quadratic programm formulation for solving an non symmetric LCP

(see also the functions/solvers list in ``LCP_Solvers.h`` and numbering in ``lcp_cst.h`` )

