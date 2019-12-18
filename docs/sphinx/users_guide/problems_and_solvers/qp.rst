.. index::
   single: Quadratic Programming problems (QP)

.. contents::

.. _qp_problem:

Quadratic Programming problems (QP)
***********************************

Problem statement
=================

Minimize:

.. math::

    \frac{1}{2} x' C x + d' x

subject to:

.. math::

    \begin{eqnarray*}
    A(j)*x + b(j) = 0 & , & j=1,...,me \\
    A(j)*x + b(j) >= 0 & , & j=me+1,...,m \\
    xl <= x <= xu \end{eqnarray*}

Available solvers
=================

The qp pack is not yet implemented. The only available function is :func:`ql0001()` (fortran subroutine).


.. seealso::

   :ref:`convexqp_problem`.
