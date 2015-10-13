Modelisation
------------

.. _bouncingball-model:

Usage :
^^^^^^^
Example of the modelisation of a bouncing ball:
+++++++++++++++++++++++++++++++++++++++++++++++

A ball bouncing on the ground may be defined as a linear Lagrangian time
invariant dynamical system with one degree of freedom. 

.. math::

      M \ddot q = F_{ext} + p


where:
 - :math:`q` is the state vector. Here, we have :math:`q \in \mathbb{R}^{1}` as the only required coordinate
   corresponds to the height of the center of the ball.
 - :math:`M` is the time invariant mass matrix. In this example :math:`M` is a 1 x 1 matrix.
 - :math:`F_{ext}` contains the external forces. Here, the gravity is applied.
 - :math:`p` is the reaction force due to the nonsmooth interaction with the floor.

We first import the needed classes from `siconos.kernel` module:

 - `LagrangianLinearTIDS`, for a Linear Lagrangian Time Invariant Dynamical System object.
 - `LagrangianLinearTIR`, for a Linear Lagrangian Time Invariant Relation object.
 - `NewtonImpactNSL`, for a Newton Impact NonSmooth Law object.
 - `Interaction`, to build an object that glues the relation and the nonsmooth law.
 - `Model`, to build an object that gathers all the modeling and simulation objects.

.. testcode::

  from siconos.kernel import \
      LagrangianLinearTIDS, LagrangianLinearTIR, NewtonImpactNSL, \
      Interaction, Model

To build a `LagrangianLinearTIDS` object, we have to give an initial
position vector, an initial velocity vector and the constant mass
matrix of the object.

The position and velocity are both vectors of size 1. The mass
matrix is defined as a 1 x 1 matrix.

.. testcode::
  
  position = [1]    # initial position vector
  velocity = [0]    # initial velocity vector
  mass = 1          # mass of the ball
  M = [[mass]]      # mass matrix 
  
  ball = LagrangianLinearTIDS(position, velocity, M)

The gravity is expressed in the coordinates chosen for the ball. It is then
applied as a constant external force.

.. testcode::

  g = 9.81  
  weight = [- mass * g]    # a vector of size 1
  ball.setFExtPtr(weight)  # apply gravity

The ball is constrained to lie above the floor. The relation between
the state space and the constraint space is a linear mapping :math:`y
= C * q` where :math:`y` denotes the constraint vector and :math:`q`
denotes the state vector. :math:`C` is a 1 x 1 matrix: 
:math:`C = \{1\}`

We build an object for this relation with the `LagrangianLinearTIR` class:

.. testcode::

  C = [[1]]           
  relation = LagrangianLinearTIR(C)

The "above floor" constraint is unilateral and defined by :math:`y
> 0` and a relation between velocities before and after impact. Let
:math:`t^{-}` denotes the time instant before impact and :math:`t^{+}` denotes
time instant after impact. Then if we add to the unilateral constraint a linear
relation between the velocities :math:`\dot y(t^{+}) = e * \dot
y(t^{-})`, we define a Newton impact nonsmooth law:

.. testcode::

  e = 0.9
  nslaw = NewtonImpactNSL(e)

The relation and the nonsmooth law are tied together in an `Interaction`
object:

.. testcode::

  inter = Interaction(nslaw,    
                      relation) 

We finally build a `Model` object to gather the dynamical systems we
have defined (here just the ball) and link the interactions to them.

.. testcode::

  # the first parameter is the start time
  # the second parameter is the end time
  bouncingBall = Model(0, 10)

  # add the ball to the model data
  bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)

  # link the interaction and the dynamical system
  bouncingBall.nonSmoothDynamicalSystem().link(inter, ball)


.. _diodes-bridge-model:

Example of the modelisation of a diodes bridge:
+++++++++++++++++++++++++++++++++++++++++++++++

This is an example of an electrical circuit involving:

  - a linear dynamical system consisting of an LC oscillator (1 µF ,
    10 mH)
  - a non smooth system (a 1000 Ohm resistor supplied through a 4
    diodes bridge) in parallel with the oscillator

.. image:: /figures/electronics/DiodeBridge/SchemaDiodeBridge.*

Expected behavior:

The initial state (:math:`v_c = 10` V , :math:`i_L = 0`) of the oscillator provides
an initial energy. The period is :math:`2 \pi \sqrt(LC) \equiv 0,628` ms.

The non smooth system is a full wave rectifier: each phase (positive
and negative) of the oscillation allows current to flow through the
resistor in a constant direction, resulting in an energy loss: the
oscillation damps.

State variables :
  - the voltage across the capacitor (or inductor)
  - the current through the inductor

Since there is only one dynamical system, the interaction is defined
by :

 - complementarity laws between diodes current and voltage. Depending
   on the diode position in the bridge, :math:`y` stands for the reverse
   voltage across the diode or for the diode current (see figure in
   the template file)
 - a linear time invariant relation between the state variables and :math:`y`
   and :math:`\lambda` (derived from Kirchhoff laws)

The oscillator is a time-invariant linear dynamical system, and using
the Kirchhoff current and voltage laws and branch constitutive
equations, its dynamics is written as:

.. math::
   
   \dot x = A.x + r,
   

where :math:`x = \left[\begin{array}{c}  v_L\\ i_L \end{array}\right]`, :math:`\lambda = \left[\begin{array}{c} -v_{DR1}\\ -v_{DF2}\\ i_{DF1}\\ i_{DR2} \end{array}\right]`, :math:`A=\left[\begin{array}{cc} 0 & \frac{-1}{C}\\ \frac{1}{L} & 0 \end{array}\right]` and :math:`r= \left[\begin{array}{cccc} 0 & 0 & \frac{-1}{C} & \frac{1}{C}\\ 0 & 0 & 0 & 0 \end{array}\right]\lambda`


We first import the needed classes for the construction of the model:

.. testcode::

   from siconos.kernel import FirstOrderLinearDS, FirstOrderLinearTIR, \
       ComplementarityConditionNSL, Interaction,\
       Model

We define the model parameters:

.. testcode::

   Lvalue = 1e-2    # inductance
   Cvalue = 1e-6    # capacitance
   Rvalue = 1e3     # resistance
   Vinit = 10.0     # initial voltage

A `FirstOrderLinearDS` object is built with the initial state and the linear
operator :math:`A`:

.. testcode::
    
   initstate = [Vinit, 0]

   A = [[0,          -1.0/Cvalue],
        [1.0/Lvalue, 0          ]]

   diodesBridge = FirstOrderLinearDS(initstate, A)


The linear relations between voltage and current given by:

:math:`y = Cx + D\lambda`, where :math:`y=\left[ \begin{array}{c} i_{DR1}\\ i_{DF2}\\ -v_{DF1}\\ -v_{DR2} \end{array} \right]`, :math:`C = \left[ \begin{array}{cc} 0 & 0\\ 0 & 0\\ -1 & 0\\ 1 & 0 \end{array} \right]` :math:`D = \left[ \begin{array}{cccc} \frac{1}{R} & \frac{1}{R} & -1 & 0\\ \frac{1}{R} & \frac{1}{R} & 0 & -1\\ 1 & 0 & 0 & 0\\ 0 & 1 & 0 & 0 \end{array} \right]` and :math:`\lambda =\left[ \begin{array}{c} -v_{DR1}\\ -v_{DF2}\\ i_{DF1}\\ i_{DR2} \end{array} \right]`

are encoded as a `FirstOrderLinearTIR` object:

.. testcode::

   C = [[0.,   0.],
        [0,    0.],
        [-1.,  0.],
        [1.,   0.]]

   D = [[1./Rvalue, 1./Rvalue, -1.,  0.],
        [1./Rvalue, 1./Rvalue,  0., -1.],
        [1.,        0.,         0.,  0.],
        [0.,        1.,         0.,  0.]]

   B = [[0.,        0., -1./Cvalue, 1./Cvalue],
        [0.,        0.,  0.,        0.       ]]

   relation = FirstOrderLinearTIR(C, B)
   relation.setDPtr(D)

The behavior of each diode of the bridge, supposed to be ideal,
may be described with a complementarity condition between current and
reverse voltage (variables (:math:`y`, :math:`\lambda`) ). Depending on the
diode position in the bridge, y stands for the reverse voltage across
the diode or for the diode current. Then, the complementarity
conditions, results of the ideal diodes characteristics are given
by:

:math:`\begin{array}{l} 0 \leq -v_{DR1} \, \perp \, i_{DR1} \geq 0\\0 \leq -v_{DF2} \, \perp \, i_{DF2} \geq 0\\0 \leq i_{DF1} \, \perp \, -v_{DF1} \geq 0\\0 \leq i_{DR2} \, \perp \, -v_{DR2} \geq 0\end{array} \ \ \ \ \ \ or \ \ \ \ \ \  0 \leq y \, \perp \, \lambda \geq 0`

These conditions are defined with a `ComplementarityConditionNSL` object of size 4:

.. testcode::
   
   nslaw = ComplementarityConditionNSL(4)

Then the relation and the nonsmooth law are tied in an `Interaction` object:

.. testcode::

   interaction = Interaction(nslaw, relation)

Finally a `Model` object is built :

.. testcode::

   # the first parameter is the start time
   # the second parameter is the end time
   DiodesBridgeModel = Model(0, 10)

.. testcode::

   # add the dynamical system to the model data
   DiodesBridgeModel.nonSmoothDynamicalSystem().insertDynamicalSystem(diodesBridge)

   # link the interaction and the dynamical system
   DiodesBridgeModel.nonSmoothDynamicalSystem().link(interaction, diodesBridge)

Modelisation API
^^^^^^^^^^^^^^^^

.. automodule:: siconos.kernel
  :members: :eval:`under_directory(['src/modelingTools'])`
