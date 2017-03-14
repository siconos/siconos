# Siconos
|  **`ubuntu-14-04-gcc-atlas-lapack`**   |
|----------------------------------------|
|[![Build Status](https://travis-ci.org/siconos/siconos.svg?branch=master)](https://travis-ci.org/siconos/siconos)|

A software package for the modeling and simulation of nonsmooth dynamical systems in C++ and in Python.

Siconos is an open-source scientific software primarily targeted at modeling and simulating nonsmooth dynamical systems:

  * _Mechanical systems_ (rigid or solid) with unilateral contact and Coulomb friction and impact (Nonsmooth mechanics, 
contact dynamics, multibody systems dynamics or granular materials). 
  * _Switched Electrical Circuit_ such as electrical circuits with ideal and piecewise linear components: power converter, rectifier, Phase-Locked Loop (PLL) or Analog-to-Digital converter.
  * _Sliding mode control_ systems.
  * _Biology_ Gene regulatory networks.
 
Other applications are found in Systems and Control (hybrid systems, differential inclusions,
optimal control with state constraints), Optimization (Complementarity systems and Variational inequalities), 
Fluid Mechanics, Computer Graphics, ...

Read more about Siconos at the [Siconos homepage](http://siconos.gforge.inria.fr)
# Main components

Each component can be used either from a low-level language like C/C++ or from Python.

## siconos/numerics (C)

Collection of low-level algorithms for solving optimization problems arising in the simulation of nonsmooth dynamical systems:

  * Complementarity problems ([LCP](https://en.wikipedia.org/wiki/Linear_complementarity_problem), [MLCP](https://en.wikipedia.org/wiki/Mixed_linear_complementarity_problem), [NCP](https://en.wikipedia.org/wiki/Nonlinear_complementarity_problem))
  * Friction-contact problems (2D or 3D)
  * Second-order cone programming (SOCP)
  * Primal or Dual Relay problems
  * Finite dimensional [Variational Inequality](https://en.wikipedia.org/wiki/Variational_inequality) (AVI and VI)

## siconos/kernel (C++)

Library for the modeling and simulation of nonsmooth dynamical systems.

  * Dynamical systems formalism: first order systems, Lagrangian and Newton-Euler formulations
  * Numerical integration techniques: Event-detecting (event-driven) and Event-Capturing (time-stepping) schemes
  * Nonsmooth laws: complementarity, Relay, normal cone inclusion, Friction Contact, Newton impact, multiple impact law.

## siconos/mechanics (C++)

Component for the simulation of mechanical systems in interaction with their environment:
* Contact detection procedure between simple primitives (homemade) and meshes [bullet3](https://github.com/bulletphysics/bullet3)
* Contact detection between Brep representation based on [oce. Open CASCADE Community Edition](https://github.com/tpaviot/oce) and on [pythonOCC](https://github.com/tpaviot/pythonocc) 3D CAD/CAM package for python 

## siconos/control (C++)

Library to add a controller to a simulation. For now almost all the implemented control schemes are based on sliding modes with an implicit discretization.

## siconos/io (C++)

This component can be used to 
* serialize almost any simulation using [boost::serialization](http://www.boost.org/doc/libs/1_60_0/libs/serialization/doc/index.html)
* generate mechanical examples from HDF5 and to write HDF5 in view of visualization through [vtk](http://www.vtk.org)

# License

Siconos is currently distributed under Apache Licenses (v2).

### The archetypal example: "The bouncing ball"
```python
from siconos.kernel import LagrangianLinearTIDS, NewtonImpactNSL,\
LagrangianLinearTIR, Interaction, Model, MoreauJeanOSI,\
TimeDiscretisation, LCP, TimeStepping
from numpy import eye, empty

t0 = 0       # start time
T = 10       # end time
h = 0.005    # time step
r = 0.1      # ball radius
g = 9.81     # gravity
m = 1        # ball mass
e = 0.9      # restitution coeficient
theta = 0.5  # theta scheme

# the dynamical system
x = [1, 0, 0]    # initial position
v = [0, 0, 0]    # initial velocity
mass = eye(3)  # mass matrix
mass[2, 2] = 2. / 5 * r * r
ball = LagrangianLinearTIDS(x, v, mass)
weight = [-m * g, 0, 0] 
ball.setFExtPtr(weight) #set external forces
# Interaction ball-floor
H = [[1, 0, 0]]
nslaw = NewtonImpactNSL(e)
relation = LagrangianLinearTIR(H)
inter = Interaction(nslaw, relation)
# Model
bouncingBall = Model(t0, T)
# add the dynamical system to the non smooth dynamical system
bouncingBall.nonSmoothDynamicalSystem().insertDynamicalSystem(ball)
# link the interaction and the dynamical system
bouncingBall.nonSmoothDynamicalSystem().link(inter, ball)
# Simulation
# (1) OneStepIntegrators
OSI = MoreauJeanOSI(theta)
# (2) Time discretisation 
t = TimeDiscretisation(t0, h)
# (3) one step non smooth problem
osnspb = LCP()
# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(t, OSI, osnspb)
# end of model definition

# computation
bouncingBall.setSimulation(s)
bouncingBall.initialize() # simulation initialization
N = (T - t0) / h # the number of time steps
# time loop
while s.hasNextEvent():
    s.computeOneStep()
    s.nextStep()
