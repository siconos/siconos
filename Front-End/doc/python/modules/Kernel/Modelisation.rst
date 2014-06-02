Modelisation
------------

Usage : example of the modelisation of a bouncing ball
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A ball bouncing on the ground may be defined as a linear lagrangian time
invariant dynamical system with one degree of freedom. We first import the
needed classes from `Siconos.Kernel` module:

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



The ball-floor relation is defined as a time invariant and linear lagrangian
relation :math:`y= Cq + e + D\lambda + Fz` 

where :math:`C = {1}` [...]

.. testcode::

  C = [[1]]           
  relation = LagrangianLinearTIR(C)

A nonsmooth law defined as a Newton impact law : 

.. testcode::

  e = 0.9
  nslaw = NewtonImpactNSL(e)



We build an interaction object with the size of the input and the output, the nonsmooth law and the relation

.. testcode::

  # the first parameter is the size of the input and output
  inter = Interaction(1,      
                      nslaw,    
                      relation) 



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
