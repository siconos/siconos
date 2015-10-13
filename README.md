Siconos
================
A software for modeling and simulation of nonsmooth dynamical systems in C++ and in Python.

Siconos is an open-source scientific software primarily targeted at modeling and simulating nonsmooth dynamical systems :
 * Mechanical systems (rigid or solid) with unilateral contact and Coulomb friction and impact (Nonsmooth mechanics, 
contact dynamics, multibody systems dynamics or granular materials). 
 * Switched Electrical Circuit such as electrical circuits with ideal and piecewise linear components: power converter, rectifier, Phase-Locked Loop (PLL) or Analog-to-Digital converter.
 * Sliding mode control systems.
 * Biology (Gene regulatory network). 
 Other applications are found in Systems and Control (hybrid systems, differential inclusions,
optimal control with state constraints), Optimization (Complementarity systems and Variational inequalities), 
Fluid Mechanics, Computer graphics, ...

Siconos contains
   * Siconos/Numerics. Collection of low-level algorithms for solving basic Algebra and optimization problem arising in the simulation of nonsmooth dynamical systems:
     * Complementarity problems (LCP, MLCP, NCP)
     * Friction-contact problems (2D or 3D)
     * Second-order cone programming (SOCP)
     * Primal or Dual Relay problems
     * Variational inequality (VI)
   * Siconos/Kernel.  Modeling and simulation of the nonsmooth dynamical systems. it contains :
     * Dynamical systems classes : first order and Lagrangian systems, Newton-Euler systems
     * Nonsmooth laws : complementarity, Relay, FrictionContact, impact
     * Numerical integration techniques : Event-detecting (event-driven) and Event-Capturing schemes
   * Siconos/Mechanics
   * Siconos/Control
   * Siconos/IO
