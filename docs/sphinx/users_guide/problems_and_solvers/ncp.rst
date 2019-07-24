.. index:: single: Nonlinear Complementarity Problems (NCP)
.. _doxid-_n_c_problem:

Nonlinear Complementarity Problems (NCP)
========================================

.. _doxid-_n_c_problem_1ncpIntro:
.. rubric:: The problem:

Find :math:`z \in \mathcal{R}^n_+` such that:

.. math::

    \begin{equation*} 0 \le z \perp F(z) \ge 0 \end{equation*}

.. _doxid-_n_c_problem_1ncpSolvers:
.. rubric:: Available solvers/formulations::

* :func:`ncp_newton_FBLSA()` with the FB merit function and a Newton with line-search

* :func:`ncp_newton_minFBLSA()` with the min merit function (with the FB as backup) and a Newton with line-search

* :func:`ncp_pathsearch()` solver using a path search

* NCP_Path() Interface to Path (Ferris)

.. _doxid-_n_c_problem_1ncpProblemIntro:
.. rubric:: Problem Statement:

Given a sufficiently smooth function :math:`{F}\colon {{\mathrm{I\!R}}}^{n} \to {{\mathrm{I\!R}}}^{n}` The Nonlinear Complementarity Problem (NCP) is to find two vectors :math:`(z,w \in {{\mathrm{I\!R}}}^{n})` such that:

.. math::

    \begin{align*} w &= F(z) \\ 0 &\le w \perp z \ge 0 \end{align*}

.. _doxid-_n_c_problem_1ncpSolversList:
.. rubric:: Available solvers::

* ncp_FBLSA(), nonsmooth Newton method based on Fisher-Burmeister function with a line search.

* :func:`ncp_pathsearch()` , a solver based on a path search method

