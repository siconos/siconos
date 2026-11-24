.. _architecture:

############
Architecture
############

This page gives developers an overview of how Siconos is built: its components,
what each one contains, how they depend on each other, and what happens during
one time step of a simulation. For modelling and numerical details, see the
:ref:`siconos_users_guide`.

Overview
========

Siconos is a library for modelling and simulating **nonsmooth dynamical
systems**: mechanical systems with impacts and friction, electrical circuits with
ideal diodes or switches, sliding-mode control, and similar systems. It is split
into components. Each component builds one shared library, plus Python bindings
written with pybind11.

.. figure:: /figures/siconos_components.*
   :width: 600
   :alt: Siconos components and their dependencies

   Siconos components and their dependencies.

Components
==========

+----------------+----------+-----------------------------------------------------------+
| Component      | Language | Content                                                   |
+================+==========+===========================================================+
| ``externals``  | C, C++,  | Third-party code shipped with Siconos: BLAS/LAPACK        |
|                | Fortran  | wrappers, ODEPACK/Hairer integrators, lp_solve, PATH SDK, |
|                |          | boost add-ons, …                                          |
+----------------+----------+-----------------------------------------------------------+
| ``numerics``   | C        | Formulation and solvers of the one-step nonsmooth         |
|                |          | problems (LCP, MLCP, friction contact, …), sparse and     |
|                |          | dense matrix storage, generic solver options.             |
+----------------+----------+-----------------------------------------------------------+
| ``kernel``     | C++20    | Modelling of nonsmooth dynamical systems, time            |
|                |          | integration and simulation strategies. Linear algebra     |
|                |          | is built on Eigen.                                        |
+----------------+----------+-----------------------------------------------------------+
| ``control``    | C++      | Sensors, observers, controllers and actuators (sliding    |
|                |          | mode control, twisting, PID…) used together with kernel   |
|                |          | simulations.                                              |
+----------------+----------+-----------------------------------------------------------+
| ``mechanics``  | C++      | Rigid bodies and contact detection (Bullet or native),    |
|                |          | joints, OpenCASCADE contact shapes, finite elements,      |
|                |          | cohesive zone models.                                     |
+----------------+----------+-----------------------------------------------------------+
| ``mechanisms`` | C++      | Optional. Multibody systems defined from CAD files        |
|                |          | (MBTB/CADMBTB, needs OpenCASCADE).                        |
+----------------+----------+-----------------------------------------------------------+
| ``io``         | C++,     | Serialization (boost), HDF5 files for mechanical          |
|                | Python   | simulations (``mechanics_hdf5.py``, the                   |
|                |          | ``mechanics_run.py`` runner) and the ``vview`` viewer.    |
+----------------+----------+-----------------------------------------------------------+
| ``vkernel``    | C++      | Optional, experimental. Struct-of-arrays storage.         |
+----------------+----------+-----------------------------------------------------------+

Each component depends only on the components above it in the table
(``externals`` → ``numerics`` → ``kernel`` → ``control``, ``mechanics`` →
``io``). ``numerics`` can be used on its own as a C library. The kernel only uses
it through the one-step nonsmooth problem classes (see `Coupling between kernel
and numerics`_).

Kernel
======

The kernel source tree (``kernel/src``) is split into the following C++
namespaces:

``siconos::algebra`` (``algebra/``)
   ``SiconosVector``, ``BlockVector``, ``SiconosMatrix`` (dense, Eigen
   ``MatrixXd``), sparse matrices and ``SiconosMemory`` (state history). Eigen
   does **not** initialize storage. Always call ``setZero()`` (or use a zeroing
   constructor) before reading a newly allocated matrix.

``siconos::modeling`` (``modeling/``)
   Building blocks of the model:

   - ``NonSmoothDynamicalSystem`` (NSDS): owns the dynamical systems and
     interactions, stored as a graph (DS = vertices, interactions = edges);
   - ``DynamicalSystem`` and its subclasses;
   - ``Interaction`` = one ``Relation`` + one ``NonSmoothLaw``, linking one or
     two dynamical systems;
   - ``Relation`` and ``NonSmoothLaw`` hierarchies;
   - ``BoundaryCondition``.

``siconos::integrators`` (``simulation/``)
   One-step integrators (OSI): time discretisation of the smooth dynamics of
   the dynamical systems.

``siconos::nonsmooth_formulations`` (``simulation/``)
   One-step nonsmooth problems (OSNS): the problem solved at each step (or
   event) to compute the multipliers λ of the interactions.

``siconos::simulation`` (``simulation/``)
   ``Simulation`` strategies (time-stepping, event-driven), ``TimeDiscretisation``,
   ``EventsManager``, ``Topology`` and ``InteractionManager`` (dynamic
   creation of interactions, for example by contact detection).

``utils/``
   Exceptions (``siconos::exception``), graphs, quaternions, serialization
   helpers.

The Python modules ``siconos.modeling``, ``siconos.integrators``,
``siconos.nonsmooth_formulations`` and ``siconos.simulation`` follow these
namespaces.

Main class hierarchies
----------------------

The tree below shows the main classes only. The Doxygen reference lists all of
them.

.. code-block:: text

   DynamicalSystem
   ├── FirstOrderNonLinearDS
   │   └── FirstOrderLinearDS
   ├── SecondOrderDS
   │   ├── LagrangianDS
   │   │   ├── LagrangianLinearTIDS
   │   │       └── LagrangianLinearDiagonalDS
   │   ├── LagrangianSparseDS ── LagrangianSparseLinearTIDS
   │   └── NewtonEulerDS            (rigid bodies, quaternion parametrization)

   Relation
   ├── FirstOrderR
   │   ├── FirstOrderNonLinearR, FirstOrderType1R, FirstOrderType2R
   │   └── FirstOrderLinearR, FirstOrderLinearTIR
   ├── LagrangianR
   │   ├── LagrangianScleronomousR (── Lagrangian2d*R), LagrangianRheonomousR
   │   ├── LagrangianCompliantR, LagrangianCompliantLinearTIR
   │   └── LagrangianLinearTIR
   ├── LagrangianSparseR            (for LagrangianSparseDS)
   └── NewtonEulerR
       ├── NewtonEulerNormalVelocityR
       │   └── NewtonEulerVelocityR, NewtonEulerVelocityAngularVelocityR
       │       └── (mechanics) ContactR
       └── (mechanics) NewtonEulerJointR: PivotJointR, PrismaticJointR, …

   NonSmoothLaw
   ├── ComplementarityConditionNSL, MixedComplementarityConditionNSL,
   │   EqualityConditionNSL, RelayNSL, NormalConeNSL
   ├── NewtonImpactNSL, MultipleImpactNSL
   ├── NewtonImpactFrictionNSL, NewtonImpactRollingFrictionNSL,
   │   FremondImpactFrictionNSL
   └── CohesiveZoneModelNIFNSL, MohrCoulombPlasticityNSL

   Simulation
   ├── TimeStepping
   │   ├── TimeSteppingCombinedProjection
   │   └── TimeSteppingDirectProjection
   ├── TimeSteppingD1Minus          (with D1MinusLinearOSI)
   └── EventDriven

   OneStepIntegrator
   ├── EulerMoreauOSI              (first-order systems, θ-method)
   ├── MoreauJeanOSI               (Moreau–Jean scheme, local formulation)
   │   ├── MoreauJeanGOSI          (global formulation)
   │   └── MoreauJeanCombinedProjectionOSI, MoreauJeanDirectProjectionOSI
   ├── MoreauJeanBilbaoOSI
   ├── SchatzmanPaoliOSI, D1MinusLinearOSI, NewMarkAlphaOSI
   ├── ZeroOrderHoldOSI            (linear systems, used by control)
   └── LsodarOSI, Hem5OSI          (event-driven integrators)

   OneStepNSProblem
   ├── LinearOSNS
   │   ├── LCP, MLCP (── MLCPProjectOnConstraints), Relay, AVI, Equality
   │   ├── FrictionContact ── CohesiveFrictionContact
   │   ├── RollingFrictionContact
   │   ├── GlobalFrictionContact ── GlobalRollingFrictionContact
   │   ├── GenericMechanical
   │   └── MultipleImpact
   └── QP

Relations and dynamical systems are dispatched through visitors
(``DynamicalSystemVisitor.hpp``, ``RelationVisitor.hpp``). A new concrete class
must be added to the relevant visitor, otherwise ``accept`` falls back to the
base class overload, which throws.

Simulation strategies
---------------------

Siconos has two families of simulation strategies (see :ref:`simulation_nsds`):

**Event-capturing (time-stepping)**: ``TimeStepping`` and its variants
   Fixed time step. Impulses (not forces) are computed at each step, so impacts
   and accumulations of impacts need no event detection. Used with all OSIs
   except ``LsodarOSI``/``Hem5OSI``. Nonlinear systems are solved with a Newton
   loop (``TimeSteppingType``: ``LINEAR``, ``LINEAR_IMPLICIT``, ``NONLINEAR``,
   ``NONLINEAR_FULL``). ``TimeSteppingCombinedProjection`` and
   ``TimeSteppingDirectProjection`` also project positions onto the
   constraints.

**Event-driven**: ``EventDriven``
   The smooth dynamics are integrated with an ODE/DAE solver (``LsodarOSI``,
   ``Hem5OSI``) up to the next event, which is located with root finding. At
   events, an impact problem (velocity level) and a contact problem
   (acceleration level) are solved. This strategy is accurate when there are
   few events, but cannot go through accumulations of impacts (Zeno behaviour).

Numerics
========

``numerics`` is a C library. Each problem type has a problem structure, a
driver, and a set of solvers selected through ``SolverOptions``
(``core/SolverOptions.h``): a solver id, integer and double parameters, and
optional internal solvers.

+--------------------------------+---------------------------------------------+------------------------------+
| Problem                        | Description                                 | Driver                       |
+================================+=============================================+==============================+
| LCP                            | Linear complementarity problem              | ``linearComplementarity_     |
|                                |                                             | driver``                     |
+--------------------------------+---------------------------------------------+------------------------------+
| MLCP                           | Mixed LCP (equalities + complementarity)    | ``mlcp_driver``              |
+--------------------------------+---------------------------------------------+------------------------------+
| MCP, NCP                       | (Mixed) nonlinear complementarity problems  | ``mcp_driver``,              |
|                                |                                             | ``ncp_driver``               |
+--------------------------------+---------------------------------------------+------------------------------+
| AVI, VI                        | (Affine) variational inequalities           | ``avi_driver``,              |
|                                |                                             | ``variationalInequality_     |
|                                |                                             | driver``                     |
+--------------------------------+---------------------------------------------+------------------------------+
| Relay                          | Relay (box-constrained) problem             | ``relay_driver``             |
+--------------------------------+---------------------------------------------+------------------------------+
| ConvexQP                       | Convex quadratic program                    | ``convexQP_*`` solvers       |
+--------------------------------+---------------------------------------------+------------------------------+
| SOCLCP                         | Second-order cone LCP                       | ``soclcp_driver``            |
+--------------------------------+---------------------------------------------+------------------------------+
| FrictionContact (2D/3D)        | Coulomb friction, reduced (local)           | ``fc2d_driver``,             |
|                                | formulation                                 | ``fc3d_driver``              |
+--------------------------------+---------------------------------------------+------------------------------+
| GlobalFrictionContact (2D/3D)  | Coulomb friction, global formulation        | ``gfc2d_driver``,            |
|                                | (velocities and reactions)                  | ``gfc3d_driver``             |
+--------------------------------+---------------------------------------------+------------------------------+
| RollingFrictionContact         | Coulomb friction with rolling resistance    | ``rolling_friction_2d/3d_    |
|                                | (local and global)                          | driver``, ``global_rolling_  |
|                                |                                             | friction_3d_driver``         |
+--------------------------------+---------------------------------------------+------------------------------+
| CohesiveFrictionContact        | Friction contact with cohesion              | ``cohesive_friction_3d_      |
|                                |                                             | driver``                     |
+--------------------------------+---------------------------------------------+------------------------------+
| Plasticity (2D)                | Plasticity (von Mises, Mohr–Coulomb)        | ``plasticity_2d_driver``     |
+--------------------------------+---------------------------------------------+------------------------------+
| GenericMechanical              | Mix of the above, one block per interaction | ``gmp_driver``               |
+--------------------------------+---------------------------------------------+------------------------------+

The exact driver signatures are in ``numerics/src/NonSmoothDrivers.h``.

The solvers belong to a few families, which are shared between problem types
when possible:

- **pivoting methods**: Lemke, enumeration, simplex (LCP, MLCP, Relay);
- **splitting methods**: projected/nonsmooth Gauss–Seidel (PGS, NSGS), PSOR.
  For friction problems, NSGS solves one contact at a time with a local solver
  (projection, local Newton…). The common logic is in
  ``NonSmoothSolvers/nsgs_*``;
- **projection and fixed-point methods**: fixed-point projection,
  extra-gradient, hyperplane projection (VI), De Saxcé and Panagiotopoulos
  fixed points, proximal point, ADMM;
- **semi-smooth Newton methods**: Alart–Curnier, Jean–Moreau (natural map),
  Fischer–Burmeister, with line search or Levenberg–Marquardt stabilisation;
- **interior point methods** (global and rolling friction contact);
- **external solvers**: PATH, GAMS (optional).

Matrices are ``NumericsMatrix`` (``linalg/NumericsMatrix.h``) with three
storage types: ``NM_DENSE`` (column-major), ``NM_SPARSE_BLOCK`` (sparse block
matrix ``SparseBlockStructuredMatrix``, for one block per interaction) and
``NM_SPARSE`` (CSparse compressed column/triplet). Most solvers work on all
three. Look at the ``NM_*`` functions before writing code for one storage
type only.

See :ref:`problems_and_solvers` for each problem and the list of its solvers.

Coupling between kernel and numerics
------------------------------------

Each ``OneStepNSProblem`` of the kernel builds the matching numerics problem
structure and calls its driver:

- ``LinearOSNS::computeDiagonalInteractionBlock`` and
  ``computeInteractionBlock`` assemble the operator (Delassus matrix
  W = H M⁻¹ Hᵀ for the local formulations) block by block, one block per pair
  of interactions. ``OSNSMatrix`` stores it as a ``NumericsMatrix``;
- ``computeqBlock`` assembles, for each interaction, the free-motion term q
  from the OSI's free state;
- ``compute`` fills the numerics problem (operator, q, friction
  coefficients…), calls the driver (for example ``linearComplementarity_driver``
  in ``LCP.cpp``, ``fc3d_driver`` in ``FrictionContact.cpp``) and copies the
  solution into the interactions' y and λ;
- the ``SolverOptions`` are owned by the OSNS and can be changed from C++ or
  Python (``osnspb.numericsSolverOptions()``).

One time step
=============

The sequence below is the time-stepping loop of
``kernel/src/simulation/TimeStepping.cpp``. Event-driven simulations reuse the
same building blocks between events.

1. **Initialization** (``Simulation::initialize``, first call): link the OSIs
   with the dynamical systems, create the interactions' work vectors, compute
   the index sets, allocate the OSNS matrices. Each later step calls it again
   to take into account systems or interactions added to the NSDS (for
   example by an ``InteractionManager``).
2. ``computeOneStep()`` → ``advanceToEvent()`` → ``newtonSolve()``. Each Newton
   iteration:

   a. ``prepareNewtonIteration``: OSIs compute their iteration matrices
      (for example W = M + hθC + h²θ²K for Moreau–Jean with a linear DS);
   b. ``computeFreeState``: state without the nonsmooth contribution;
   c. ``computeOneStepNSProblem``: assemble and solve the OSNS for λ;
   d. ``updateInput``: r (or p) = Hᵀλ, the nonsmooth input to the DS;
   e. ``updateState``: OSIs compute the new state from the free state and r;
   f. ``updateOutput``: y = h(x, t), then ``updateIndexSets``;
   g. ``newtonCheckConvergence`` (nonlinear case only; linear systems take
      one iteration).

3. ``nextStep()``: save the state in memory, move the time discretisation and
   the events forward.

Mechanics, control and io
=========================

**mechanics** (``mechanics/src``)
   - ``collision/``: ``RigidBodyDS`` (Newton–Euler), contact shapes,
     ``SiconosCollisionManager`` (an ``InteractionManager``) with a Bullet
     backend (``bullet/``) and a native backend (``native/``);
     ``ContactR`` and variants;
   - ``joints/``: ``NewtonEulerJointR`` and pivot, prismatic, cylindrical,
     fixed, coupler, knee joints, joint stops and joint friction;
   - ``occ/``: OpenCASCADE bodies and contact shapes (optional);
   - ``fem/``: finite elements (native solid/beam elements, cables);
   - ``czm/``: cohesive zone models.

**control** (``control/src``)
   ``Sensor`` → ``Observer`` → ``Actuator``/controller, scheduled by a
   ``ControlManager`` through events of the ``EventsManager``.
   ``ControlZOHSimulation`` and ``ControlLsodarSimulation`` drive a plant with
   its controllers.

**io** (``io/src``, ``io/python``)
   Boost serialization of a whole simulation (``SiconosRestart``),
   ``MechanicsIO`` (contact points and forces), and the Python HDF5 layer:
   scenes are written with ``MechanicsHdf5Runner``, simulated with
   ``mechanics_run.py`` and shown with ``siconos_vview``.

Python bindings
===============

Bindings are written with pybind11 (``*_pybind11.cpp`` files in
``<component>/python``) and installed as the ``siconos`` package:

.. code-block:: text

   siconos.numerics                numerics problems, drivers, SolverOptions
   siconos.modeling                DS, relations, nonsmooth laws, NSDS
   siconos.integrators             one-step integrators
   siconos.nonsmooth_formulations  one-step nonsmooth problems
   siconos.simulation              TimeStepping, EventDriven, TimeDiscretisation
   siconos.mechanics.{collision, joints, occ, fem, czm}
   siconos.io                      HDF5 runner, vview, tools

pybind11 converts arguments to temporary objects. A C++ object must not keep a
reference to an argument: it must store a copy, or the binding must use
``py::keep_alive``. User functions (external forces, nonlinear relations) are
``std::function`` objects, so they can be Python callables. The former plugin
mechanism (C functions loaded by name) has been removed.

Bouncing ball example (full version in the tutorials,
``examples/mechanics/BouncingBall/BouncingBallTS.py``):

.. code-block:: python

   import numpy as np
   import siconos.modeling as sm
   import siconos.integrators as si
   import siconos.nonsmooth_formulations as snf
   import siconos.simulation as ss

   t0, T, h = 0.0, 10.0, 0.005
   r, m, g, e = 0.1, 1.0, 9.81, 0.9

   # Dynamical system: a ball (3 dof), with gravity
   mass = np.eye(3, order="F")
   mass[2, 2] = 2.0 / 5 * r * r
   ball = sm.LagrangianLinearTIDS(np.array([1.0, 0, 0]), np.zeros(3), mass, sm.alias_t)
   ball.setConstantFext(np.array([-m * g, 0, 0]), sm.alias_t)

   # Interaction: unilateral contact with the ground + Newton impact law
   relation = sm.LagrangianLinearTIR(np.array([[1.0, 0, 0]], order="F"))
   inter = sm.Interaction(sm.NewtonImpactNSL(e), relation)

   nsds = sm.NonSmoothDynamicalSystem(t0, T)
   nsds.insertDynamicalSystem(ball)
   nsds.link(inter, ball)

   # Simulation: Moreau-Jean scheme + LCP at each step
   sim = ss.TimeStepping(nsds, ss.TimeDiscretisation(t0, h),
                         si.MoreauJeanOSI(0.5), snf.LCP())

   while sim.hasNextEvent():
       sim.computeOneStep()
       # read ball.q(), ball.velocity(), inter.lambda_python(1) ...
       sim.nextStep()

Directory structure
===================

.. code-block:: text

   siconos/
   ├── externals/        third-party code (blas_lapack wrappers, netlib, hairer, PATH_SDK…)
   ├── numerics/
   │   ├── src/
   │   │   ├── LCP/ MLCP/ MCP/ NCP/ AVI/ VI/ QP/ Relay/ SOCP/
   │   │   ├── FrictionContact/ RollingFrictionContact/ CohesiveFrictionContact/
   │   │   ├── GenericMechanical/ Plasticity/
   │   │   ├── NonSmoothSolvers/   generic NSGS and nonsmooth Newton
   │   │   ├── linalg/             NumericsMatrix, sparse block, CSparse
   │   │   ├── core/               SolverOptions, error handling
   │   │   ├── optimization/ tools/ utils/ io/
   │   │   └── NonSmoothDrivers.h
   │   └── python/
   ├── kernel/
   │   ├── src/
   │   │   ├── algebra/            vectors, matrices (Eigen)
   │   │   ├── modeling/           DS, relations, nonsmooth laws, NSDS
   │   │   ├── simulation/         OSI, OSNS, simulation strategies
   │   │   └── utils/
   │   └── python/
   ├── control/          src/{Controller,Observer,Sensor,Simulation,utils}, python/
   ├── mechanics/        src/{collision,joints,occ,fem,czm}, python/
   ├── mechanisms/       src/{MBTB,CADMBTB}   (optional)
   ├── io/               src/{serialization,mechanics}, python/
   ├── vkernel/          experimental
   ├── python/           top-level python package (siconos/__init__.py, …)
   ├── cmake/            build system modules
   ├── config_samples/   option files for -DUSER_OPTIONS_FILE
   ├── ci_gitlab/        CI configuration and Docker images
   └── docs/             Sphinx + Doxygen documentation

Unit tests are in a ``test/`` (or ``tests/``) directory next to the sources they
test: cppunit in the kernel, C tests generated from solver/problem collections in
numerics, pytest in the ``python/tests`` directories.

Build system
============

Siconos is built with CMake:

.. code-block:: bash

   cmake -S <sources> -B <build> -DUSER_OPTIONS_FILE=<options file>
   cmake --build <build> -j
   cmake --install <build>
   ctest --test-dir <build>

Option files are in ``config_samples/`` (for example ``siconos_release.cmake``,
or ``siconos_asan.cmake`` for a debug build with sanitizers). Components are
enabled with ``WITH_<component>`` options, for example
``-DWITH_mechanisms=ON``. Dependencies between components are resolved
automatically. Optional dependencies also have ``WITH_*`` options:
``WITH_BULLET``, ``WITH_OpenCASCADE``, ``WITH_MUMPS``, ``WITH_FCLIB``,
``WITH_SANITIZER``….

The install location is derived from the active Python environment (virtualenv
or conda), so the C++ libraries and the Python package are installed together.
See :ref:`siconos_install_guide`.

Further reading
===============

- :ref:`modeling_nsds` and :ref:`simulation_nsds`: modelling and simulation
  concepts;
- :ref:`time_integrators` and :ref:`osns_problems`: OSIs and OSNS in detail;
- :ref:`problems_and_solvers`: numerics problems and solvers;
- :ref:`siconos_algebra`: matrix storage in the kernel and in numerics;
- :ref:`siconos_devel_guide`: other developer documentation.
