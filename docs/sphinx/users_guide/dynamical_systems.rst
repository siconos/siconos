.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _dynamical_systems:


Dynamical Systems
=================

:cpp:class:`DynamicalSystem` is the class used in Siconos to describe a set of ordinary differential equations, the first step of any nonsmooth problem description in Siconos.
This base class defines an interface common to all systems. Several derived classes provide specific formulations for different types of problems, as described below.

.. image:: /figures/dynamical_system_classes.*

As usual, a complete description of the interface (members and methods) of these classes can be found in the Doxygen documentation, see for example :cpp:class:`DynamicalSystem`.

Note that :cpp:class:`DynamicalSystem` is an abstract class: no object of this type can be instantiated. It only provides a generic interface for all systems.


Overview
--------

The most general way to write dynamical systems in Siconos is

.. math::

  g(\dot x, x, t, z) = 0

an n-dimensional set of equations, where

* t is the time;
* :math:`x \in R^{n}` is the state;
* :math:`\dot x` is the time derivative of the state;
* :math:`z \in R^{s}` is a vector of arbitrary algebraic variables, a sort of discrete state.
  For example, z may be used to set perturbation parameters, or anything else;
* :math:`g : \mathbb{R}^{n} \times \mathbb{R}^n \times \mathbb{R} \times \mathbb{R}^s \to \mathbb{R}^{n}`.

Under some specific conditions, we can rewrite this as:

.. math::

   \dot x = rhs(x, t, z)

"rhs" means right-hand side.
Note that in this case :math:`\nabla_{\dot x} g` must be invertible.

This class provides members and functions for all types of dynamical systems (i.e. for all derived classes), with specific behaviours depending on the type of system (see the related sections below for details).

*This means that all the members and functions described below are also available in all the derived classes.*

Each system is identified by a number, and the current state of the system is saved as a vector :cpp:func:`DynamicalSystem::x`, with x[0] = :math:`x` and x[1] = :math:`\dot x`.

All the functions and their gradients (:math:`g, rhs, \nabla_x g`...) can be accessed with functions like :cpp:func:`DynamicalSystem::jacobianRhsOver_x` for :math:`\nabla_{x} rhs(x, t, z)`. See the reference for the complete list of members and methods.

The common rules for all members are, 'name' being the variable:

* getName() to get a copy of the content of the object
* name() to get a pointer to the object
* setName(obj) to copy obj into Name
* setNamePtr(objPtr) to link objPtr with Name

Plugins: some members can be connected to user plugin functions that compute them. In that case, the following methods can be used:

* setComputeNameFunction(...) to link name to your own function
* computeName(...) to compute name with your own function

For details about the plugin mechanism, see :ref:`siconos_plugins`.

For instance, the internal forces of Lagrangian systems (see below) can be set in two ways: either the forces are a constant vector, or they are connected to a plugin and can then depend on time, on the state of the system, etc.

First case::

  // we suppose that ds is an existing pointer to a LagrangianDS
  auto myF = std::make_shared<siconos::algebra::SimpleVector>(3);
  // fill myF ...
  ds->setFInt(*myF); // copy the values of myF into fInt
  // OR
  // link fInt to myF: any change in one of them affects the other.
  ds->setFIntPtr(myF);

Second case::

  // we suppose that ds is an existing pointer to a LagrangianDS
  // and that myFunction is a C function implemented in myPlugin.cpp
  ds->setComputeFInt("myPlugin", "myFunction");
  // ...
  ds->computeFint(qdot, q, time);
  // compute fInt value at time for the current state

Note that the signature (i.e. the number and type of arguments) of the function in your plugin must be exactly the same as the one given in kernel/src/plugin/DefaultPlugin.cpp for the corresponding function.


Common interface
----------------

The following functions are (and must be) present in any class derived from DynamicalSystem:

* :cpp:func:`DynamicalSystem::initRhs()`

* :cpp:func:`DynamicalSystem::icomputeRhs()`

* :cpp:func:`DynamicalSystem::computeJacobianRhsOver_x()`

* :cpp:func:`DynamicalSystem::initializeNonSmoothInput()`

* :cpp:func:`DynamicalSystem::swapInMemory()`

* :cpp:func:`DynamicalSystem::display()`

* :cpp:func:`DynamicalSystem::resetAllNonSmoothParts()`

* :cpp:func:`DynamicalSystem::resetNonSmoothPart()`



First-order dynamical systems
-----------------------------

Nonlinear
"""""""""

:cpp:class:`FirstOrderNonLinearDS`

They are described by the following set:

.. math::

   M\dot x(t) &= f(t,x,z) + r \\
   x(t_0)&=x_0

with:

* :math:`M \in \mathbb{R}^{n \times n}`
* f(x,t): the vector field - :math:`f: \mathbb{R}^{n} \times \mathbb{R} \to \mathbb{R}^n`
* r: input due to the nonsmooth behaviour - vector of size n.

* JacobianXF = :math:`\nabla_x f(t,x,z)`, an n x n square matrix, is also a member of the class.

* M is assumed to be invertible (otherwise x[1] = rhs cannot be computed).

* The initial conditions are given by the member x0, a vector of size n. This is the value of x when the simulation starts,
  i.e. after a call to the simulation initialize() function.

* There are plugin functions in this class for f and its Jacobian, jacobianfx.

We have:

.. math::

   rhs &= M^{-1}(f(t,x,z)+r) \\
   \nabla_x rhs &= M^{-1}\nabla_x f(t,x,z)

The other variables are those of the :cpp:class:`DynamicalSystem` class, but some of them are not defined and thus not usable:

* g and its gradients

Linear
""""""

:cpp:class:`FirstOrderLinearDS`

Described by a set of n equations and initial conditions:

.. math::

   \dot x(t) &= A(t,z)x(t)+ b(t,z)+r \\
   x(t_0)&=x_0

With:

* A(t,z): n x n matrix, state-independent but possibly time-dependent.
* b(t,z): vector of size n, possibly time-dependent.

A and b have corresponding plugin functions.
The other variables are those of the :cpp:class:`DynamicalSystem` and FirstOrderNonLinearDS classes, but some of them are not defined and thus not usable:

* g and its gradients
* f and its gradient

And we have:

.. math::

   rhs &= M^{-1}(A(t,z)x(t)+b(t,z)) \\
   \nabla_x rhs&= M^{-1}(A(t,z)

Second-order (Lagrangian) systems
---------------------------------

Nonlinear
"""""""""

:cpp:class:`LagrangianDS`, derived from :cpp:class:`DynamicalSystem`.

Nonlinear second-order Lagrangian systems are described by the following set of nDof equations and initial conditions:

.. math::

   Mass(q,z) \ddot q &= f_L(t,\dot q , q , z) + p \\
   q(t_0) &= q0 \\
   \dot q(t_0) &= velocity0

with:

* Mass(q,z): nDof x nDof inertia matrix.
* q: state of the system - vector of size nDof.
* :math:`\dot q`: time derivative of the state.
* :math:`f_L(t,\dot q , q , z) =  F_{Ext}(t,z) - fGyr(\dot q, q,z) - F_{Int}(t,\dot q , q , z)`
* :math:`fGyr(\dot q, q,z)`: nonlinear, time-independent terms - vector of size nDof.
* :math:`F_{Int}(t,\dot q , q , z)`: time-dependent internal forces - vector of size nDof.
* :math:`F_{Ext}(t,z)`: external forces, time-dependent BUT state-independent - vector of size nDof.
* p: input due to the nonsmooth behaviour - vector of size nDof.

Note that the decomposition of :math:`f_L` only makes the interface more convenient for users; it does not affect the simulation process.

Some gradients are also required:

* jacobianFInt[0] = :math:`\nabla_q F_{Int}(t,q,\dot q,z)` - nDof x nDof matrix.
* jacobianFInt[1] = :math:`\nabla_{\dot q} F_{Int}(t,q,\dot q,z)` - nDof x nDof matrix.
* jacobianfGyr[0] = :math:`\nabla_q fGyr(\dot q, q, z)` - nDof x nDof matrix.
* jacobianfGyr[1] = :math:`\nabla_{\dot q}fGyr(\dot q, q, z)` - nDof x nDof matrix.

We assume that the mass matrix is invertible and that its gradient is zero.

There are plugin functions in this class for :math:`F_{int}, F_{Ext}, M, fGyr` and the four Jacobian matrices.

The other variables are those of the :cpp:class:`DynamicalSystem` class, but some of them are not defined and thus not usable:

* g and its gradients

The links with :cpp:class:`DynamicalSystem` are :math:`n= 2 ndof` and :math:`x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]`.

And we have:

.. math::

   rhs = \left[
   \begin{array}{c}
   \dot q \\
   Mass^{-1}(f_L(t,\dot q , q , z)+p)
   \end{array}\right]

   \nabla_x rhs = \left[
   \begin{array}{cc}
   0 & I \\
   Mass^{-1}\nabla_{q}f_L(t,\dot q , q , z) & Mass^{-1}\nabla_{\dot q}f_L(t,\dot q , q , z)
   \end{array}\right]

I: identity matrix.

Linear and time-invariant
"""""""""""""""""""""""""

Class LagrangianLinearTIDS, derived from LagrangianDS.

.. math::

   Mass \ddot q + C \dot q + K q =  F_{Ext}(t,z) + p

With:

* C: constant nDof x nDof damping matrix
* K: constant nDof x nDof stiffness matrix

The other variables are those of the :cpp:class:`DynamicalSystem` and LagrangianDS classes, but some of them are not defined and thus not usable:

* g and its gradients
* fL, fInt, fGyr and their gradients.

And we have:

.. math::

   rhs = \left[
   \begin{array}{c}
   \dot q \\
   Mass^{-1}(F_{ext}(t,z)- Kq - C\dot q+p)
   \end{array}\right]

   \nabla_x rhs = \left[
   \begin{array}{cc}
   0 & I \\
   -Mass^{-1}K & -Mass^{-1}C
   \end{array}\right]

.. _dsPlugins:

Plugin functions of the dynamical systems
-----------------------------------------

* :cpp:class:`DynamicalSystem`: :math:`g(t,\dot x,x,z), \ \ \nabla_x g(t,\dot x,x,z), \ \ \nabla_{\dot x} g(t,\dot x,x,z)`
* :cpp:class:`FirstOrderNonLinearDS`: :math:`f(t,x,z), \ \ \nabla_x f(t,x,z)`
* :cpp:class:`FirstOrderLinearDS`: A(t,z), b(t,z)
* :cpp:class:`LagrangianDS`: :math:`M(q,z), \ \ fGyr(\dot q,q,z), \ \ F_{Int}(t,\dot q,q ,z), \ \ F_{Ext}(t,z), \ \ \nabla_q F_{Int}(t,\dot q,q,z), \ \ \nabla_{\dot q}F_{Int}(t,\dot q, q, z), \ \ \nabla_q fGyr(\dot q, q, z), \ \ \nabla_{\dot q}fGyr(\dot q, q, z)`.
* :cpp:class:`LagrangianLinearTIDS`: :math:`F_{Ext}(t,z)`
