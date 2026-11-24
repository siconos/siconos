.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _diode_bridge_example:



Tutorial: a four-diode bridge wave rectifier
============================================

.. highlight:: python

.. note:: The up-to-date version of this tutorial is available as a Python notebook in the `siconos-tutorials repository <https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/-/blob/main/siconos-notebooks/st01_diode_bridge.ipynb>`_.

Preamble
--------

This tutorial is dedicated to the simulation of the system shown in :ref:`fig_diode_bridge` below. We describe step by step how to model it
as a nonsmooth dynamical system and how to integrate it in time.

.. _fig_diode_bridge:

.. figure:: /figures/electronics/DiodeBridge/diodeBridge.*
   :align: center

   fig 1: Diode bridge

An LC oscillator, initialized with a given voltage
across the capacitor and zero current through the inductor, supplies
energy to a load resistor through a full-wave rectifier
consisting of a bridge of four ideal diodes. Both half-waves of the oscillating
voltage across the LC are applied to the resistor, with the current
always flowing in the same direction. The energy is dissipated in the
resistor, which results in a damped oscillation.

The behaviour of the diode is shown in :ref:`fig_diode`: the left-hand sketch shows the ideal diode characteristic and the right-hand sketch the usual exponential characteristic given by Shockley's law.

.. _fig_diode:

.. figure:: /figures/electronics/diodeNonSmooth.*
   :align: center

   fig 2: diode characteristics

The diodes, assumed to be ideal, lead to *complementarity* between voltage and current, which introduces nonsmoothness into the system. This is detailed later.


Siconos driver file
-------------------

To run a simulation in Siconos, you write a "driver" file, either in C++ (this file must be compiled and linked with the Siconos libraries)
or in Python, and execute it. In both cases, the whole process is handled by the 'siconos' script::

  siconos DiodeBridge.cpp

compiles, links and executes the program, while::

  siconos DiodeBridge.py

executes your Python script.

Siconos can of course be used in an interactive Python session or notebook. This is probably the easiest way to follow this tutorial.

*Remarks:*

    * *this example is available in the Siconos examples, in the Electronics/DiodeBridge directory, in C++, in Python and as a notebook.
      For the latter, run*::

	ipython notebook DiodeBridge.ipynb

      *to interactively run this tutorial.*

    * *In this tutorial, we assume that Siconos is properly installed on your system, as explained in* :ref:`siconos_install_guide`.

Let us start with a short description of the three main steps always required to run a simulation.

#. First, describe the system as a nonsmooth dynamical system, i.e.:

   * define the set of ordinary differential equations (the **Dynamical Systems**) that represent the dynamics;
   * define the 'nonsmooth' part of the system, through **nonsmooth laws** and **relations** between variables that may constrain the state.

#. Then choose a simulation strategy, i.e. how the nonsmooth system is integrated over a time step: which discretisation and integrators for the dynamics (**one-step integrators**), which formulation and solvers for the **nonsmooth problem**, and so on.
#. Finally, run the simulation and post-process the results.



.. contents::
   :backlinks: entry

Building a nonsmooth dynamical system
-------------------------------------

Modeling the dynamics
^^^^^^^^^^^^^^^^^^^^^

The oscillator is a linear dynamical system with time-invariant coefficients. Using Kirchhoff's current and voltage laws and the
branch constitutive equations, the dynamics of the system read

.. math::

   \left[\begin{array}{c}
   \dot v_L\\
   \dot i_L
   \end{array}\right]=
   \left[\begin{array}{cc}
   0 & \frac{-1}{C}\\
   \frac{1}{L} & 0
   \end{array}\right].
   \left[\begin{array}{c}
   v_L\\
   i_L
   \end{array}\right]
   +
   \left[\begin{array}{cccc}
   0 & 0 & \frac{-1}{C} & \frac{1}{C}\\
   0 & 0 & 0 & 0
   \end{array}\right].
   \left[\begin{array}{c}
   -v_{DR1}\\
   -v_{DF2}\\
   i_{DF1}\\
   i_{DR2}
   \end{array}\right]

and if we denote

.. math::

   x = \left[\begin{array}{c}
   v_L\\
   i_L
   \end{array}\right], \
   \lambda = \left[\begin{array}{c}
   -v_{DR1}\\
   -v_{DF2}\\
   i_{DF1}\\
   i_{DR2}
   \end{array}\right], A=\left[\begin{array}{cc}
   0 & \frac{-1}{C}\\
   \frac{1}{L} & 0
   \end{array}\right], r= \left[\begin{array}{cccc}
   0 & 0 & \frac{-1}{C} & \frac{1}{C}\\
   0 & 0 & 0 & 0
   \end{array}\right].\lambda

we get a first-order linear system

.. math::

   \dot x = A x + r

with the unknowns :math:`x` and :math:`r`.

To represent this kind of ordinary differential equation, Siconos has the class :cpp:class:`siconos::modeling::FirstOrderLinearDS`,
which inherits from :cpp:class:`siconos::modeling::DynamicalSystem`. See :ref:`dynamical_systems` for a complete review of the dynamical system formalisms available in the software.

::

   # import siconos package
   import siconos.modeling as sm
   # numpy for vectors and matrices
   import numpy as np

   # dynamical system parameters
   Lvalue = 1e-2    # inductance
   Cvalue = 1e-6    # capacitance
   Rvalue = 1e3     # resistance
   Vinit = 10.0     # initial voltage
   x0 = [Vinit, 0.] # initial state
   # A matrix of the linear oscillator
   A = np.zeros((2, 2), dtype=np.float64)
   A.flat[...] = [0., -1.0/Cvalue, 1.0/Lvalue, 0.]

   # build the dynamical system
   ds = sk.FirstOrderLinearDS(x0, A)

A few remarks:

* in Python, you can use either lists or numpy arrays to build the vectors or matrices passed to Siconos methods;
* help on Siconos objects is available with the standard Python help function. For example, to find how
  the system can be built::

    help(sk.FirstOrderLinearDS)

  or in the :ref:`cpp_api` or :ref:`python_api` documentation.

Modeling the interactions
^^^^^^^^^^^^^^^^^^^^^^^^^

Now the nonsmooth part of the system must be defined, namely the nonsmooth laws and the constraints between the variables.
In Siconos, a **nonsmooth law** together with a **relation** between one or two dynamical systems is called an **Interaction**
(see :ref:`interactions`).
A set of dynamical systems and the interactions between them thus make up the complete nonsmooth dynamical system.

For the oscillator of :ref:`fig_diode_bridge`, there are linear relations (constraints) between the voltages and currents of the diodes, given by

.. math::

   \left[ \begin{array}{c}
   i_{DR1}\\
   i_{DF2}\\
   -v_{DF1}\\
   -v_{DR2}
   \end{array} \right]
   =
   \left[ \begin{array}{cc}
   0 & 0\\
   0 & 0\\
   -1 & 0\\
   1 & 0
   \end{array} \right]
   \cdot
   \left[ \begin{array}{c}
   v_L\\
   i_L
   \end{array} \right]
   +
   \left[ \begin{array}{cccc}
   \frac{1}{R} & \frac{1}{R} & -1 & 0\\
   \frac{1}{R} & \frac{1}{R} & 0 & -1\\
   1 & 0 & 0 & 0\\
   0 & 1 & 0 & 0
   \end{array} \right]
   \cdot
   \left[ \begin{array}{c}
   -v_{DR1}\\
   -v_{DF2}\\
   i_{DF1}\\
   i_{DR2}
   \end{array} \right]

with

.. math::

   y=\left[ \begin{array}{c}i_{DR1}\\i_{DF2}\\-v_{DF1}\\-v_{DR2}\end{array} \right], C= \left[ \begin{array}{cc}0 & 0\\ 0 & 0\\ -1 & 0\\ 1 & 0\end{array} \right], D=\left[ \begin{array}{cccc}\frac{1}{R} & \frac{1}{R} & -1 & 0\\ \frac{1}{R} & \frac{1}{R} & 0 & -1\\1 & 0 & 0 & 0\\0 & 1 & 0 & 0\end{array} \right]

and recalling that

.. math::

   \lambda = \left[\begin{array}{c}
   -v_{DR1}\\
   -v_{DF2}\\
   i_{DF1}\\
   i_{DR2}
   \end{array}\right],
   B = \left[\begin{array}{cccc}
   0 & 0 & \frac{-1}{C} & \frac{1}{C}\\
   0 & 0 & 0 & 0
   \end{array}\right], \ r = B\lambda


this is equivalent to the linear relation between :math:`(x, r)` and :math:`(y, \lambda)`:

.. math::

   y = Cx + D\lambda, \ \ \ r = B\lambda

To represent this kind of algebraic equation, Siconos has the class :cpp:class:`siconos::modeling::FirstOrderLinearTIR` (TIR stands for time-invariant relation),
which inherits from :cpp:class:`siconos::modeling::Relation`. See :ref:`relations` for a complete review of the relation formalisms available in the software.

::

   # --- build an interaction ---
   interaction_size = 4  # number of constraints

   # B, C, D matrices of the relation
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

   # set relation type
   relation= sk.FirstOrderLinearTIR(C, B)
   relation.setConstantD(D)

   # set nonsmooth law
   nonsmooth_law = sk.ComplementarityConditionNSL(interaction_size)

   # nslaw + relation == interaction
   interaction = sk.Interaction(nonsmooth_law, relation)

Note that the complete form of a :cpp:class:`siconos::modeling::FirstOrderLinearTIR` is

.. math::

   y &= Cx + D\lambda + Fz + e \\
   r &= B\lambda

All the components not set at construction are zero (here F and e).

Each diode of the bridge is assumed to be ideal, with the behaviour shown in the left-hand sketch of :ref:`fig_diode`.
Such a behaviour can be described by a **complementarity condition** between the current and the reverse voltage.

Complementarity between two variables :math:`y \in R^m, \lambda \in R^m` reads

.. math::

   y \geq 0, \ \lambda \geq 0, \ \text{ and } \ y > 0 \ \Rightarrow \ \lambda = 0, \ \ \lambda > 0 \ \Rightarrow \ y = 0

or, using the ":math:`\perp`" symbol,

.. math::

   0 \leq y\, \perp \, \lambda \geq 0\\

The inequalities hold component-wise.

Back to our circuit: the complementarity conditions given by the ideal diode characteristics are:

.. math::

   \begin{array}{l}
   0 \leq -v_{DR1} \, \perp \, i_{DR1} \geq 0\\
   0 \leq -v_{DF2} \, \perp \, i_{DF2} \geq 0\\
   0 \leq i_{DF1} \, \perp \, -v_{DF1} \geq 0\\
   0 \leq i_{DR2} \, \perp \, -v_{DR2} \geq 0\\
   \end{array} \ \ \ \ \ \ \text{ or } \ \ \ \ \ \  0 \leq y \, \perp \, \lambda \geq 0

with :math:`y` and :math:`\lambda` defined above. Note that, depending on the position of the diode in the bridge, :math:`y_i` is either the reverse voltage across
the diode or the diode current.

To represent such a nonsmooth law, Siconos has the class :cpp:class:`ComplementarityConditionNSL` (NSL appears in the name of every class defining a nonsmooth law),
which inherits from :cpp:class:`NonSmoothLaw`. See :ref:`non_smooth_laws` for a complete review of the laws available in the software.

::

   nonsmooth_law = sk.ComplementarityConditionNSL(interaction_size)

The interaction can now be fully defined::

  interaction = sk.Interaction(nonsmooth_law, relation)

Note that, for the moment, this interaction only describes relations and laws: it is not connected to any dynamical system yet.

The modeling part is almost complete, since only one dynamical system and one interaction are needed to describe the problem.
They must be gathered into a specific object, the :cpp:class:`NonSmoothDynamicalSystem`.
Building this object is simple:
set the time window of the simulation, insert the dynamical systems and link them to the right interactions.

::

   # dynamical systems and interactions must be gathered into a NonSmoothDynamicalSystem
   t0 = 0. # initial time
   T = 5.0e-3 # duration of the simulation
   DiodeBridge = sk.NonSmoothDynamicalSystem(t0, T)
   # add the dynamical system to the nonsmooth dynamical system
   DiodeBridge.insertDynamicalSystem(ds)

   # link the interaction and the dynamical system
   DiodeBridge.link(interaction, ds)


Describing the simulation of the nonsmooth dynamical system
-----------------------------------------------------------

You now need to define how the nonsmooth dynamical system is integrated over time. This is the role of the simulation, which defines:

* how the dynamical systems are discretized and integrated over a time step;
* how the nonsmooth problem is formulated and solved.

Two strategies are available: event-capturing (a.k.a. time-stepping) schemes and event-driven schemes. See :ref:`simulation_nsds` for details,
or :cite:`Acary.Brogliato2008` for even more details.

For the diode bridge example, we use an event-capturing strategy, with an Euler-Moreau integrator and an LCP (Linear Complementarity Problem)
formulation.

Let us start with the 'one-step integrator', i.e. the discretisation and integration of the dynamics over a time step, between
times :math:`t_i` and :math:`t_{i+1}`. The integration of the equation over the time step is based on a :math:`\theta`-method. The process is detailed in :ref:`event_capturing` and, for first-order systems, leads to

.. math::

   x_{i+1} &= x^{free}_{i}  +  hW^{-1}r_{i+1} \\
   W &= (M - h\theta A) \\
   x^{free}_i &=  x_i + h W^{-1}(A x_i + b)

This corresponds to the :cpp:class:`EulerMoreauOSI` integrator, which inherits from :class:`OneStepIntegrator`. See :ref:`time_integrators` for a complete review of the integrators available in the software.

::

  theta = 0.5
  osi = sk.EulerMoreauOSI(theta)

Note that each dynamical system of the model must be associated with one and only one one-step integrator.

Next, depending on the simulation strategy and the time integration, a one-step nonsmooth problem must be formulated, see :ref:`simulation_nsds` and :ref:`osns_problems`.

With the following discretization of the relations and nonsmooth law defined above

.. math::

   y_{i+1} &= C(t_{i+1})x_{i+1} + D(t_{i+1})\lambda_{i+1} \\
   R_{i+1} &= B(t_{i+1})\lambda_{i+1}\\
   0 \leq y_{i+1}\ &\perp  \lambda_{i+1} \geq 0  \\

we get

.. math::

   y_{i+1} &= q + M\lambda_{i+1} \\
   0 \leq y_{i+1}\ &\perp  \lambda_{i+1} \geq 0  \\

with

.. math::

   q = Cx^{free}_{i+1},    M = hCW^{-1}B + D

This is a Linear Complementarity Problem, represented in Siconos by the :cpp:class:`LCP` class, which inherits from :class:`OneStepNSProblem`.
As usual, see :ref:`osns_problems` for a complete review of the nonsmooth problem formulations available in Siconos.
Each formulation must be associated with a solver, picked from the list given in :ref:`lcp_solvers`::

  import siconos.numerics as sn
  #  Nonsmooth problem
  osnspb = sk.LCP(sn.solver_ids.SICONOS_LCP_NSQP)

Note that the solvers come from Siconos numerics and are identified by an id. The correspondence between ids and solvers is given in :ref:`lcp_solvers`.

The last step is the creation of the simulation, with its time discretisation::

  # simulation and time discretisation
  time_step =  1.0e-6
  td = sk.TimeDiscretisation(t0, time_step)
  simu = sk.TimeStepping(DiodeBridge, td, osi, osnspb)


Running the simulation
----------------------

The easiest way to run your simulation is to call::

  s->run()

But then you only have access to the values computed at the last
time step, which may not be enough.

Here, :math:`x, y` and :math:`\lambda` are needed at each time step for post-processing. Here is how to get them and save them in a numpy array::

  N = int((T - t0) / simu.timeStep()) + 1
  data_plot = np.zeros((N, 8))
  y = interaction.y(0)
  lamb = interaction.lambda_(0)
  x = ds.x()
  k = 0
  data_plot[k, 1] = x[0] #  inductor voltage
  data_plot[k, 2] = x[1] # inductor current
  data_plot[k, 3] = y[0] # diode R1 current
  data_plot[k, 4] = - lamb[0] # diode R1 voltage
  data_plot[k, 5] = - lamb[1] # diode F2 voltage
  data_plot[k, 6] = lamb[2] # diode F1 current
  data_plot[k, 7] = y[0] + lamb[2] # resistor current
  while simu.hasNextEvent():
      k += 1
      simu.computeOneStep()
      data_plot[k, 0] = simu.nextTime()
      data_plot[k, 1] = x[0]
      data_plot[k, 2] = x[1]
      data_plot[k, 3] = y[0]
      data_plot[k, 4] = - lamb[0]
      data_plot[k, 5] = - lamb[1]
      data_plot[k, 6] = lamb[2]
      data_plot[k, 7] = y[0] + lamb[2]
      simu.nextStep()


* :cpp:func:`Simulation::hasNextEvent()` is true as long as there are events to process, i.e. until T is reached;
* :cpp:func:`Simulation::nextStep()` mainly increments the time step, saves the current state and prepares the initial values of the next step;
* :cpp:func:`Simulation::computeOneStep()` performs the computation over the current time step. With Moreau's time-stepping scheme, it first integrates the dynamics to
  get the so-called free state, i.e. without nonsmooth effects, then formulates and solves an LCP, and finally integrates the dynamics again with the LCP results.

The results can now be post-processed, with matplotlib for example::

  import matplotlib.pyplot as plt
  plt.subplot(411)
  plt.title('inductor voltage')
  plt.plot(data_plot[0:k - 1, 0], data_plot[0:k - 1, 1])
  plt.grid()
  plt.subplot(412)
  plt.title('inductor current')
  plt.plot(data_plot[0:k - 1, 0], data_plot[0:k - 1, 2])
  plt.grid()
  plt.subplot(413)
  plt.title('diode R1 (blue) and F2 (green) voltage')
  plt.plot(data_plot[0:k - 1, 0], -data_plot[0:k - 1, 4])
  plt.plot(data_plot[0:k - 1, 0], data_plot[0:k - 1, 5])
  plt.grid()
  plt.subplot(414)
  plt.title('resistor current')
  plt.plot(data_plot[0:k - 1, 0], data_plot[0:k - 1, 7])
  plt.grid()


.. image:: /figures/electronics/DiodeBridge/diodeBridgeResult.*
   :align: center
