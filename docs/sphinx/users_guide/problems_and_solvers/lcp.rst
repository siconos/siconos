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

Use the generic functions :func:`lcp_driver_DenseMatrix` to call one the the specific solvers listed below:

.. _doxid-_l_c_problem_1lcpDirectSolvers:
.. rubric:: Direct solvers:

* :func:`lcp_lexicolemke` , direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

* :func:`lcp_pivot` , generic solver for pivot-based methods: Bard, Murty and Lemke rules are implemented.

* :func:`lcp_enum()` , enumerative solver (brute-force method which tries every possible solution).

* :func:`lcp_path()` , interface to the PATH solver

.. _doxid-_l_c_problem_1lcpIterativeSolvers:
.. rubric:: Iterative solvers:

* :func:`lcp_cpg()` , CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.

* :func:`lcp_pgs()` , PGS is a basic Projected Gauss-Seidel solver for LCP.

* :func:`lcp_rpgs()` , Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms.

* :func:`lcp_psor()` , Projected Succesive over relaxation solver for LCP. See Cottle, Pang Stone Chap 5.

* :func:`lcp_latin()` , (LArge Time INcrements) is a basic latin solver for LCP.

* :func:`lcp_latin_w()` , (LArge Time INcrements) is a basic latin solver with relaxation for LCP.

* :func:`lcp_nsgs_SBM()` , Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

.. _doxid-_l_c_problem_1lcpEquationBasedSolvers:
.. rubric:: Equation-based solvers:

* :func:`lcp_newton_min()` , nonsmooth Newton method based on the min formulation of the LCP.

* :func:`lcp_newton_FB()` , uses a nonsmooth newton method based on the Fischer-Bursmeister NCP function.

* :func:`lcp_newton_minFB()` , nonsmooth Newton method combining the min and FB functions.

.. _doxid-_l_c_problem_1lcpReformulation:
.. rubric:: QP-reformulation:

* :func:`lcp_qp()` , quadratic programm formulation

* :func:`lcp_nsqp()` , quadratic programm formulation for solving an non symmetric LCP

(see also the functions/solvers list in ``LCP_Solvers.h`` and numbering in ``lcp_cst.h`` )

