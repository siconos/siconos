# Siconos

[![pipeline status](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/badges/master/pipeline.svg)](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/commits/master)  [![Version](https://img.shields.io/github/release/siconos/siconos.svg)](https://github.com/siconos/siconos/releases/latest)  [![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/siconos/siconos/blob/master/COPYING)

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

## Usage

### Through your web browser

Perfect to understand and test the software or run lite simulations, demos, teaching ...

Prerequisite: a web browser and a network connection.

Pros: no installation prerequisites, nothing to do

Cons: limited resources. No long/important simulations.

Just click there [![Binder](https://plmbinder.math.cnrs.fr/binder/badge_logo.svg)](https://plmbinder.math.cnrs.fr/binder/v2/git/https%3A%2F%2Fgricad-gitlab.univ-grenoble-alpes.fr%2Fnonsmooth%2Fsiconos-tutorials.git/HEAD)

### Jupyter Lab environment with siconos ready to use and a set of end-user examples:

Prerequisite: Docker installed and usable on your computer.


To use last Siconos release, run the command line
```bash
docker run --rm -p 8888:8888 -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-master:latest
```
You will get something like
```
...
To access the server, open this file in a browser:
...
        http://127.0.0.1:8888/lab?token=b8a131ed7aed720c8fe1d7fae034cd2b669dcf686126c598
```
Copy this last line into your browser, and there you will be able to start a terminal or a notebook session.


You can also start a simple terminal to use Siconos through the command line with

```bash
docker run --rm --entrypoint /bin/bash  -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-master:latest
```



### From source

Assuming you have cloned the project into <siconos-sources>, to build and install the libraries and the python interface ::

* Create a user options file. Some templates are provided in <siconos-sources>/config_samples.
* Run

```
cmake -S <path-to-siconos-sources> -B build-siconos -DUSER_OPTIONS_FILE=<your-options-file>
cmake --build build-siconos -j 4 # or the number of cores available on your computer.
ctest --test-dir build-siconos   # run tests, optional
cmake --install build-siconos
```

More details in [Siconos download and install guide](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/install_guide/index.html).

<!-- ## Docker images -->

<!-- #### Docker images with siconos ready to use: -->

<!-- * latest version (development) -->

<!-- ``` -->
<!-- docker run -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconos-master:latest -->
<!-- ``` -->

<!-- * A specific (release) version X.Y: -->

<!-- ``` -->
<!-- docker run -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconos-release-X.Y:latest -->
<!-- ``` -->

<!-- #### Jupyter Lab environment with siconos ready to use and a set of end-user examples: -->

<!-- * latest version (development) -->

<!-- ``` -->
<!-- docker run -p 8888:8888 -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-master -->
<!-- ``` -->


<!-- * A specific (release) version X.Y: -->

<!-- ``` -->
<!-- docker run -p 8888:8888 -ti gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-release-X.Y -->
<!-- ``` -->



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
LagrangianLinearTIR, Interaction, NonSmoothDynamicalSystem, MoreauJeanOSI,\
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
bouncingBall = NonSmoothDynamicalSystem(t0, T)
# add the dynamical system to the non smooth dynamical system
bouncingBall.insertDynamicalSystem(ball)
# link the interaction and the dynamical system
bouncingBall.link(inter, ball)
# Simulation
# (1) OneStepIntegrators
OSI = MoreauJeanOSI(theta)
# (2) Time discretisation 
t = TimeDiscretisation(t0, h)
# (3) one step non smooth problem
osnspb = LCP()
# (4) Simulation setup with (1) (2) (3)
s = TimeStepping(bouncingBall, t, OSI, osnspb)
# end of model definition

# computation
N = (T - t0) / h # the number of time steps
# time loop
while s.hasNextEvent():
    s.computeOneStep()
    s.nextStep()
