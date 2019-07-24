.. index:: single: Quadratic Programming problems (QP)
.. _doxid-_q_p_solvers:

Quadratic Programming problems (QP)
===================================

.. _doxid-_q_p_solvers_1qpIntro:
.. rubric:: The problem:

Minimize:



.. math::

    \frac{1}{2} x' C x + d' x

subject to:



.. math::

    \begin{eqnarray*} A(j)*x + b(j) = 0 & , & j=1,...,me \\ A(j)*x + b(j) >= 0 & , & j=me+1,...,m \\ xl <= x <= xu \end{eqnarray*}

.. _doxid-_q_p_solvers_1qpSolversList:
.. rubric:: Available solvers:

The qp pack is not yet implemented. The only available function is ql0001() (fortran subroutine)

(see the functions/solvers list in ``QP_Solvers.h`` )

