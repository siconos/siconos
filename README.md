# Siconos

A software for modeling and simulation of nonsmooth dynamical systems in C++ and in Python.

Siconos is an open-source scientific software primarily targeted at modeling and simulating nonsmooth dynamical systems:

  * _Mechanical systems_ (rigid or solid) with unilateral contact and Coulomb friction and impact (Nonsmooth mechanics, 
contact dynamics, multibody systems dynamics or granular materials). 
  * _Switched Electrical Circuit_ such as electrical circuits with ideal and piecewise linear components: power converter, rectifier, Phase-Locked Loop (PLL) or Analog-to-Digital converter.
  * _Sliding mode control_ systems.
  * _Biology_ Gene regulatory network. 

Other applications are found in Systems and Control (hybrid systems, differential inclusions,
optimal control with state constraints), Optimization (Complementarity systems and Variational inequalities), 
Fluid Mechanics, Computer graphics, ...

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

Library for the modeling and simulation of the nonsmooth dynamical systems.

  * Dynamical systems formalism: first order systems, Lagrangian and Newton-Euler formulations
  * Numerical integration techniques: Event-detecting (event-driven) and Event-Capturing (time-stepping) schemes
  * Nonsmooth laws: complementarity, Relay, normal cone inclusion, Friction Contact, Newton impact, multiple impact law.

## siconos/mechanics (C++)

Component for the simualtion of mechanical systems in interaction with their environment:
* Contact detection procedure between simple primitives (homemade) and meshes [bullet3](https://github.com/bulletphysics/bullet3)
* Contact detection between Brep representation based on [oce. Open CASCADE Community Edition](https://github.com/tpaviot/oce) and on [pythonOCC](https://github.com/tpaviot/pythonocc) 3D CAD/CAM package for python 

## siconos/control (C++)

Library to add a controller to a simulation. For now the almost all implemented control scheme are based on sliding modes with an implicit discretization.

## siconos/io (C++)

This component can be used to 
* serialize almost any simulation using boost::archive
* generate mechanical examples from hdf5 and to write hdf5 in view of vizualization through [vtk](http://www.vtk.org)

# License

Siconos is currently distributed under GPL license (v2).
