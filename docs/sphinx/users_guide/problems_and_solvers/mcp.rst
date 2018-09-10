.. index:: single: Mixed (Non Linear) Complementarity problem (MCP)
.. _doxid-_m_c_problem:

Mixed (Non Linear) Complementarity problem (MCP)
================================================

.. _doxid-_m_c_problem_1mcpIntro:
.. rubric:: Problem Statement:

Given a sufficiently smooth function :math:`{F} \colon {{\mathrm{I\!R}}}^{n+m} \to {{\mathrm{I\!R}}}^{n+m} ` , the Mixed Complementarity problem (MCP) is to find two vectors :math:`(z,w \in {{\mathrm{I\!R}}}^{n+m})` such that:



.. math::

    \begin{align*} w &= \begin{pmatrix}w_e\\w_i\end{pmatrix} = F(z) \\ w_e &=0 \\ 0 &\le w_i \perp z_i \ge 0, \end{align*}

where "i" (resp. "e") stands for inequalities (resp. equalities). The vector :math:`z` is splitted like :math:`w` :

.. math::

    \begin{equation*}z =\begin{pmatrix}z_e\\z_i\end{pmatrix}.\end{equation*}

The vectors :math:`z_i,w_i` are of size ``sizeEqualities`` , the vectors :math:`z_e,w_e` are of size ``sizeInequalities`` and :math:`F` is a non linear function that must be user-defined.

A Mixed Complementarity problem (MCP) is a NCP "augmented" with equality constraints.

.. _doxid-_m_c_problem_1mcpSolversList:
.. rubric:: Available solvers ::

* :ref:`mcp_FB() <doxid-mcp__newton___f_b_l_s_a_8h_1afc597b5d99e5486f020864b428ecf810>` , nonsmooth Newton method based on Fisher-Burmeister function.

.. index:: single: MixedComplementarity Problems Solvers
.. _doxid-_m_c_p_solvers:

.. rubric ::MixedComplementarity Problems Solvers:

.. _doxid-_m_c_p_solvers_1mcp_FischerBurmeister:
.. rubric:: semi-smooth Newton/Fisher-Burmeister solver.:

a nonsmooth Newton method based based on the Fischer-Bursmeister convex function

function: :ref:`mcp_FischerBurmeister() <doxid-_m_c_p___solvers_8h_1a3497ef27d73b2a7b75b1ac45c12008a6>`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

