.. index::
   single: Generic mechanical problems

.. contents::

.. _gmp:

Generic mechanical problems
***************************

Complete problem with bilateral equality, complementarity, impact and friction.

Problem statement
=================

Implementation in numerics
==========================

Structure to define the problem: :cpp:class:`GenericMechanicalProblem`.

Solvers list : :cpp:enum:`GENERIC_MECHANICAL_SOLVER`

The generic drivers for generic mechanical problems is :cpp:func:`gmp_driver`.

.. _gmp_error:

Error strategy
==============

Use :cpp:func:`gmp_compute_error`.

Check details in :ref:`fc_error`


.. _gmp_solvers:

GMP solvers
===========
