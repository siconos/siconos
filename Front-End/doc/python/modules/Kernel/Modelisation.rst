Modelisation
------------

.. _bouncingball-model:

Usage : example of the modelisation of a bouncing ball
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ball bouncing on the ground may be defined as a linear lagrangian time
invariant dynamical system with one degree of freedom. 

.. math::

      M \ddot q = F_{ext} + p


where :
 - :math:`q` is the state vector, here of size 1. The only coordinate
   corresponds to the height of the center of the ball.
 - :math:`M` is the time invariant mass matrix, here of size 1 * 1.
 - :math:`F_{ext}` are the external forces, here the gravity.
 - :math:`p` is the reaction force due to the nonsmooth interaction with the floor.

We first import the needed classes from `Siconos.Kernel` module:

 - `LagrangianLinearTIDS`, for a linear lagrangian time invariant dynamical system object.
 - `LagrangianLinearTIR`, for a linear lagrangian time invariant relation object.
 - `NewtonImpactNSL`, for a Newton Impact non smooth law object.
 - `Interaction`, to build an object that will tie the relation and the nonsmooth law.
 - `Model`, to build an object that will gather the whole definitions.

.. testcode::

  from Siconos.Kernel import \
      LagrangianLinearTIDS, LagrangianLinearTIR, NewtonImpactNSL, \
      Interaction, Model

To build a `LagrangianLinearTIDS` object, we have to give an initial
position vector, an initial velocity vector and the constant mass
matrix of the object.

The position and velocity are both vectors of size 1 and the mass
matrix is here defined as a one row * one column matrix :math:`\{ 1
\}`:

.. testcode::
  
  position = [1]    # initial position vector
  velocity = [0]    # initial velocity vector
  mass = 1          # mass of the ball
  M = [[mass]]      # mass matrix 
  
  ball = LagrangianLinearTIDS(position, velocity, M)


The gravity is applied to the ball as a constant external force :

.. testcode::

  g = 9.81  
  weight = [- mass * g]    # a vector of size 1 (number of degrees of freedom)
  ball.setFExtPtr(weight)  # set the external force of the lagrangian dynamical system


The ball is constrained to lie above the floor. The relation between
the state space and the constraint space is a linear mapping :math:`y
= C * q` where :math:`y` denotes the constraint vector and :math:`q`
denotes the state vector. :math:`C` is the one row * one column matrix:
:math:`C = \{1\}`

We build an object for this relation with the `LagrangianLinearTIR` class:

.. testcode::

  C = [[1]]           
  relation = LagrangianLinearTIR(C)

The "above floor" constraint is defined by the unilateral law :math:`y
> 0` and a relation between velocities before and after impact. if
:math:`t^{+}` denotes time before impact and :math:`t^{-}` denotes
time after impact then if we add to the unilateral constraint a linear
relation between the velocities :math:`\dot y(t^{+}) = e * \dot
y(t^{-})`, we define a Newton impact nonsmooth law:

.. testcode::

  e = 0.9
  nslaw = NewtonImpactNSL(e)

The relation and the nonsmooth law are tied together in an `Interaction`
object:

.. testcode::

  # the first parameter is here the dimension of the constraint space
  inter = Interaction(1,      
                      nslaw,    
                      relation) 


We finally build a `Model` object to gather the dynamical sytems we
have defined (here just the ball) and link the interactions to them.

.. testcode::

  # the first parameter is the start time
  # the second parameter is the end time
  bouncingBall = Model(0, 10)

  # add the ball to the model's non smooth dynamical system
  bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

  # link the interaction and the dynamical system
  bouncingBall.nonSmoothDynamicalSystem().link(inter, ball)



Modelisation API
^^^^^^^^^^^^^^^^

.. automodule:: Siconos.Kernel
  :members: :eval:`under_directory(['../../../Kernel/src/modelingTools'])`
