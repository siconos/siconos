Siconos
================
A software for modeling and simulation of nonsmooth dynamical systems

SICONOS is an Open Source scientific software primarily targeted at modeling and simulating nonsmooth dynamical systems :

Mechanical systems (Rigid body or solid) with Unilateral contact and Coulomb friction as we find in Non-smooth mechanics, 
Contact dynamics or Granular material. Switched Electrical Circuit such as electrical circuits with ideal and piecewise linear
components: Power converter, Rectifier, Phase-locked loop (PLL) or Analog-to-digital converter.
Sliding mode control systems. Other applications are found in Systems and Control (hybrid systems, differential inclusions,
optimal control with state constraints), Optimization (Complementarity systems and Variational inequalities), 
Biology (Gene regulatory network), 
Fluid Mechanics, Computer graphics, ....

The software is based on 4 main components :
   * Siconos/Numerics (C API). Collection of low-level algorithms for solving basic Algebra and optimization problem arising in the simulation of nonsmooth dynamical systems:
     * Linear complementarity problems (LCP)
     * Mixed linear complementarity problems (MLCP)
     * Nonlinear complementarity problems (NCP)
     * Quadratic programming problems (QP)
     * Friction-contact problems (2D or 3D)
     * (Second-order cone programming (SOCP))
     * Primal or Dual Relay problems
   * Siconos/Kernel. API C++ used to model and simulate the NonSmooth dynamical systems. it contains :
     * Dynamical systems classes : first order and Lagrangian systems, Newton-Euler systems
     * Nonsmooth laws : complementarity, Relay, FrictionContact, impact
   * Siconos/Mechanics
   * Siconos/Control
