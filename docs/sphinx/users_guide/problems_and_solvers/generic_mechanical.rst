.. index::
   single: Generic mechanical problems

.. contents::

.. _gmp:

Generic mechanical problems
***************************

A complete problem with bilateral (equality) constraints, complementarity, impact and friction.

Problem statement
=================

Implementation in numerics
==========================

Structure to define the problem: :cpp:class:`GenericMechanicalProblem`.

List of solvers: :cpp:enum:`GENERIC_MECHANICAL_SOLVER`

The generic driver for generic mechanical problems is :cpp:func:`gmp_driver`.

.. _gmp_error:

Error strategy
==============

Use :cpp:func:`gmp_compute_error`.

See the details in :ref:`fc_error`.


.. _gmp_solvers:

GMP solvers
===========
