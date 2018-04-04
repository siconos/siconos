.. _relations:

Relations
---------

Relations are used to link local variables of the Interaction and global variables of the DynamicalSystems, and thus define constraints in the systems.

:doxysiconos:`Relation` is an abstract class which provides a generic interface for all types of relations.
Each relation has a type which correspond to the types of dynamical systems they fit with (FirstOrder or Lagrangian), a sub-type, (linear, non linear, scleronomous...).
Usually, "type+subtype" corresponds more or less to the name of the derived class.
Then, depending on the sub-class, each relation holds some plug-in functions or operators used to define the constraints. They are listed below for each available type of relation. 


Available classes: :doxysiconos:`FirstOrderR`, :doxysiconos:`FirstOrderLinearR`, :doxysiconos:`FirstOrderLinearTIR`, :doxysiconos:`LagrangianR`, :doxysiconos:`LagrangianRheonomousR`, :doxysiconos:`LagrangianScleronomousR`, :doxysiconos:`LagrangianCompliantR`, :doxysiconos:`LagrangianLinearR`.

.. image:: /figures/relation_classes.*

First Order Relations
^^^^^^^^^^^^^^^^^^^^^

Non Linear
""""""""""

Class :doxysiconos:`FirstOrderR`

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
   
:math:`h`, :math:`g` (and their jacobian according to :math:`X` and :math:`\lambda`) are defined with some plug-in functions. \n
See the doxygen documentation of the class :doxysiconos:`FirstOrderR` to have a list of the set/get/compute functions.

Note: for the signification of :math:`X`, :math:`Z`, :math:`R` see :ref:`interactions`

Linear
""""""

Class: :doxysiconos:`FirstOrderLinearR`

.. math::
   
   y &= C(t,Z)X + F(t,Z)Z + D(t,Z) \lambda + e(t,Z) \\
   R &= B(t,Z) \lambda

Plug-in functions are available for all operators.

Linear with Time Invariant Coefficients
"""""""""""""""""""""""""""""""""""""""

Class :doxysiconos:`FirstOrderLinearTIR`

.. math::
  
   y &= CX + FZ + D\lambda + e \\
   R &= B \lambda

Lagrangian (second order) Relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Scleronomous
""""""""""""

xClass :doxysiconos:`LagrangianScleronomousR`

The constraints depend only on the state,

.. math::
   
   y &= h(Q,Z) \\
   \dot y &= G_0(Q,Z)\dot Q \\
   P &= G_0^t(Q,Z)\lambda 

with
 
.. math::
    
    G_0(Q,Z) = \nabla_Q h(Q,Z)

Rheonomous
""""""""""

Class :doxysiconos:`LagrangianRheonomousR`

The constraints depend on time and state, 

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

Class: :doxysiconos:`LagrangianCompliantR`

The constraints depends on state and :math:`\lambda`, with a function of time for which :math:`\dot\lambda(t)` makes sense.

.. math::
   
   y &= h(Q,\lambda(t),Z) \\
   \dot y &= G_0(Q,\lambda(t),Z)\dot Q + G_1(Q,\lambda(t),Z)\dot\lambda(t) \\
   P &= G_0^t(Q,\lambda(t),Z)\lambda(t) 

with
 
.. math::

   G_0(Q,\lambda(t),Z) = \nabla_q h(Q,\lambda(t),Z) \\
   G_1(Q,\lambda(t),Z) = \nabla_{\lambda(t)}h(Q,\lambda(t),Z)

Linear and Time Invariant Coefficients
""""""""""""""""""""""""""""""""""""""

Class: :doxysiconos:`LagrangianLinearR`

Lagrangian linear relations with time-invariant coefficients. 

.. math::

   y&= H Q + b + D\lambda +FZ \\
   P &= H^t \lambda 

Relations plug-in functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* FirstOrderR: :math:`h(X,t,\lambda,Z), \ \ g(\lambda,t,Z)`
* FirstOrderLinearR: :math:`C(t,Z), \ \ F(t,Z), \ \ , D(t,Z), \ \ e(t,Z), B(t,Z)`
* LagrangianScleronomousR: :math:`h(Q,Z), \ \ G_0(Q,Z)`
* LagrangianRheonomousR: :math:`h(Q,t,Z), \ \ G_0(Q,t,Z), \ \ hdot(Q,t,Z)`
* LagrangianCompliantR:  :math:`h(Q,\lambda,Z), \ \ G_0(Q,\lambda,Z), \ \ G_1(Q,\lambda,Z)`
