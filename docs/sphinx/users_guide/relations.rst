.. _relations:

Relations
---------

Relations link the local variables of the Interaction to the global variables of the DynamicalSystems, and thus define the constraints of the system.

:cpp:class:`Relation` is an abstract class that provides a generic interface for all types of relations.
Each relation has a type, which corresponds to the type of dynamical systems it applies to (FirstOrder or Lagrangian), and a subtype (linear, nonlinear, scleronomous...).
Usually, "type+subtype" more or less gives the name of the derived class.
Depending on the subclass, each relation holds plugin functions or operators that define the constraints. They are listed below for each type of relation.


Available classes: :cpp:class:`FirstOrderR`, :cpp:class:`FirstOrderLinearR`, :cpp:class:`FirstOrderLinearTIR`, :cpp:class:`LagrangianR`, :cpp:class:`LagrangianRheonomousR`, :cpp:class:`LagrangianScleronomousR`, :cpp:class:`LagrangianCompliantR`, :cpp:class:`LagrangianLinearR`.

.. image:: /figures/relation_classes.*

First-order relations
^^^^^^^^^^^^^^^^^^^^^

Nonlinear
"""""""""

Class :cpp:class:`FirstOrderR`

.. math::

   output &= y = h(X,t,\lambda,Z)\\
   input &= R = g(X,t,\lambda,Z)

We denote:

.. math::

   \begin{array}{ccc}
   H_0(X,t,\lambda,Z)=\nabla_X h(X,t,\lambda,Z)&, &  H_1(X,t,\lambda,Z)=\nabla_{\lambda} h(X,t,\lambda,Z) \\
   \\
   G_0(X,t,\lambda,Z)=\nabla_X g(X,t,\lambda,Z)&, &  G_1(X,t,\lambda,Z)=\nabla_{\lambda} g(X,t,\lambda,Z)
   \end{array}

:math:`h`, :math:`g` (and their Jacobians with respect to :math:`X` and :math:`\lambda`) are defined with plugin functions.
See the Doxygen documentation of the class :cpp:class:`FirstOrderR` for the list of set/get/compute functions.

Note: for the meaning of :math:`X`, :math:`Z`, :math:`R`, see :ref:`interactions`.

Linear
""""""

Class: :cpp:class:`FirstOrderLinearR`

.. math::

   y &= C(t,Z)X + F(t,Z)Z + D(t,Z) \lambda + e(t,Z) \\
   R &= B(t,Z) \lambda

Plugin functions are available for all operators.

Linear with time-invariant coefficients
"""""""""""""""""""""""""""""""""""""""

Class :cpp:class:`FirstOrderLinearTIR`

.. math::

   y &= CX + FZ + D\lambda + e \\
   R &= B \lambda

Lagrangian (second-order) relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scleronomous
""""""""""""

Class :cpp:class:`LagrangianScleronomousR`

The constraints depend only on the state:

.. math::

   y &= h(Q,Z) \\
   \dot y &= G_0(Q,Z)\dot Q \\
   P &= G_0^t(Q,Z)\lambda

with

.. math::

    G_0(Q,Z) = \nabla_Q h(Q,Z)

Rheonomous
""""""""""

Class :cpp:class:`LagrangianRheonomousR`

The constraints depend on time and state:

.. math::

   y &= h(Q,t,Z)\\
   \dot y &= G_0(Q,t,Z)\dot Q + \frac{\partial h}{\partial t}(Q,t,Z) \\
   P &= G_0^t(Q,t,Z)\lambda

with

.. math::
   G_0(Q,t,Z) = \nabla_Q h(Q,t,Z)  \\
   hdot(Q,t,Z) = \frac{\partial h}{\partial t}(Q,t,Z)

Compliant
"""""""""

Class: :cpp:class:`LagrangianCompliantR`

The constraints depend on the state and on :math:`\lambda`, a function of time for which :math:`\dot\lambda(t)` makes sense.

.. math::

   y &= h(Q,\lambda(t),Z) \\
   \dot y &= G_0(Q,\lambda(t),Z)\dot Q + G_1(Q,\lambda(t),Z)\dot\lambda(t) \\
   P &= G_0^t(Q,\lambda(t),Z)\lambda(t)

with

.. math::

   G_0(Q,\lambda(t),Z) = \nabla_q h(Q,\lambda(t),Z) \\
   G_1(Q,\lambda(t),Z) = \nabla_{\lambda(t)}h(Q,\lambda(t),Z)

Linear with time-invariant coefficients
"""""""""""""""""""""""""""""""""""""""

Class: :cpp:class:`LagrangianLinearR`

Lagrangian linear relations with time-invariant coefficients.

.. math::

   y&= H Q + b + D\lambda +FZ \\
   P &= H^t \lambda

Plugin functions of the relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* FirstOrderR: :math:`h(X,t,\lambda,Z), \ \ g(\lambda,t,Z)`
* FirstOrderLinearR: :math:`C(t,Z), \ \ F(t,Z), \ \ , D(t,Z), \ \ e(t,Z), B(t,Z)`
* LagrangianScleronomousR: :math:`h(Q,Z), \ \ G_0(Q,Z)`
* LagrangianRheonomousR: :math:`h(Q,t,Z), \ \ G_0(Q,t,Z), \ \ hdot(Q,t,Z)`
* LagrangianCompliantR:  :math:`h(Q,\lambda,Z), \ \ G_0(Q,\lambda,Z), \ \ G_1(Q,\lambda,Z)`
