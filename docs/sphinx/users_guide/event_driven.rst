.. SPDX-FileCopyrightText: INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _event_driven:

Event-driven schemes
====================

General principle
-----------------

Event-driven schemes decompose the dynamics in time into modes, time intervals where the dynamics are smooth, and discrete events, times where the dynamics are nonsmooth. Numerically, the event-driven scheme uses this decomposition to

* detect and solve the nonsmooth dynamics at events, with a reinitialization rule for the state;
* integrate the smooth dynamics between two events with an ODE solver with root finding, possibly with bilateral constraints on the state.

Event-driven implementation
---------------------------

Integration of the smooth dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Between events, the dynamics are integrated with the Lsodar algorithm, through the function EventDriven::integrate.

Considering two given functions f(x,t) and g(x,t), a call to::

  integrate(tinit, tend, tout, iout)

integrates f(x,t) between tinit and tend and searches for the roots of the function g(x,t). If a root is found, the integration stops and the time reached is saved in tout.

The in-out parameter iout is an indicator that must be set to 1 at the first call. On output, it is 2 if no root was found, and 3 otherwise.

Thus, for an EventDriven simulation:

* the dynamics of all the DynamicalSystems involved are rewritten as :math:`\dot x = f(x,t)`;
* the relations are used to define the functions g(x,t).

The definition of f and g depends on the type of system, as described below.

Events
^^^^^^

In Siconos, an Event object only holds a type (see below) and a long int, the time at which the event occurs. It also defines a process function, whose action depends on the type of event.

The possible types (derived classes) are:

* TimeDiscretisationEvent: an event at one of the user-defined time discretisation points. These events are created and scheduled when the EventDriven simulation is initialized. The process function:

  * updates (computes) the output of all the Interactions involved;
  * saves the current values (states of the DynamicalSystems and input/output of the Interactions) in memory vectors. The last saved values become the initial values of the next integration.

* NonSmoothEvent: "points" where the dynamics are nonsmooth and need a special treatment. These events are detected by a root-finding algorithm and correspond to the violation of some constraints (relations). The process function roughly does the following (the full process depends on the type of system and is described below):

  * updates (computes) the output of all the Interactions involved;
  * updates the index sets;
  * formulates and solves one or more nonsmooth problems;
  * saves the current values (states of the DynamicalSystems and input/output of the Interactions) in memory vectors. The last saved values become the initial values of the next integration.

The Events manager
^^^^^^^^^^^^^^^^^^

All the events are handled by a dedicated object, the EventsManager. It belongs to the EventDriven class and holds two eventsContainers (sets of Events):

* pastEvents: all the Events that have already been processed;
* unProcessedEvents: the future events, already scheduled but not yet processed.

We also call "currentEvent" the last processed event, which is the starting point of the current integration, and "nextEvent" the event that follows "currentEvent".

The events manager is initialized with time-discretisation events, from the user time discretisation. Then, each time a new event is detected (added by the user, or when a root is found during integration), it is scheduled in the manager.
The manager also has a :cpp:func:`EventsManager::processEvents` function, which moves currentEvent to the set of past events, processes nextEvent and prepares the next step.

Other useful functions are :cpp:func:`EventsManager::startingTime` and :cpp:func:`EventsManager::nextTime`.

The simulation process
----------------------

A general integration step of an EventDriven simulation is thus:
"While there are events in the unProcessedEvents set, integrate the smooth dynamics between the current and the next event, then process the behaviour at the event".
In code:

.. code-block:: c++

  // Assume a NonSmoothDynamicalSystem nsds and
  // a TimeDiscretisation td ...
  auto s = std::make_shared<siconos::simulation::EventDriven>(nsds, td);

  // ...

  // We get the events manager
  auto eventsManager = s->eventsManager();

  // while there are events ...
  while(eventsManager->hasNextEvent())
  {
  // integrate between current and next event
  s->advanceToEvent();
  // solve the non-smooth dynamics, if necessary ...
  eventsManager->processEvents();
  }

  // Or in one step:
  while(eventsManager->hasNextEvent())
  {
  s->computeOneStep();
  }

.. _event_driven_lagrange:

Event-driven algorithm for Lagrangian systems
---------------------------------------------

Currently, the only event-driven algorithm available in Siconos is for Lagrangian dynamical systems, subject to perfect unilateral constraints, with Newton impact laws.

Because of the unilateral constraints, the evolution of the system may be nonsmooth. The velocity may jump and the "acceleration" may not be defined everywhere. The generalized coordinates, assumed to be absolutely continuous, are:

.. math::

   q(t) = q(t_0) +\int_{t_0}^t v^+(t)dt \ with \ v = \dot q

The indices "+" and "-" denote the right and left values of a variable at a discontinuity.

The equations of motion are written in terms of a measure differential equation:

.. math::

   M(q)dv +  F_{int}(t, q,  v^+)dt=F_{ext}(t) +  dr

r being the generalized force due to the unilateral constraints.
Using the Lebesgue decomposition theorem and its variants, the differential measures dv and dr are decomposed as:

.. math::

   dv = \gamma dt + (v^+-v^-)\sum_i\delta_{t_i} + dv_s \\
   dr = fdt + \sum_ip_i\delta_{t_i}+dr_s

The first term of the decomposition is the smooth part, with :math:`\gamma =\ddot q`, the acceleration in the usual sense. The second term is the behaviour at the times of discontinuity (:math:`\delta_{t_i}`: Dirac measure), and the last term, a singular measure, is neglected.

With these decompositions, the nonsmooth dynamics can be split into "impact equations", which correspond to the nonsmooth events, and "smooth dynamics". These equations are completed by the constraints, written at different kinematic levels, as shown in the following paragraphs.

The impact equations
^^^^^^^^^^^^^^^^^^^^

The impact equations can be written at the times :math:`t_i` of discontinuity:

.. math::

   M(q(t_i))(v^{+}(t_i)- v^{-}(t_i)) = p_i,

:math:`p_i` is an impulse.

This equation is solved at the time of impact together with an impact law, here a Newton impact law:

.. math::
   M(q(t_i))(v^{+}(t_i)- v^{-}(t_i)) = p_i, \\
   \dot y^{+}(t_i) = \nabla_q h(q(t_i)) v^{+}(t_i) \\
   \dot y^{-}(t_i) = \nabla_q h(q(t_i)) v^{-}(t_i) \\
   p_i =   \nabla_q^T h(q(t_i)) P_{N,i}\\
   0\leq  \dot y^{+}(t_i)+ e \dot y^{-}(t_i) \perp P_{N,i} \geq 0

If the matrix :math:`M(q(t_i))` is invertible, this problem can be reduced to the local unknowns :math:`\dot y^{+}(t_i),P_{N,i}`, which gives the following Linear Complementarity Problem at the time :math:`t_i` of a discontinuity of v:

.. math::
   \dot y^{+}(t_i) =  \nabla_q h(q(t_i)) (M(q(t_i)))^{-1} \nabla_q^T h(q(t_i))   P_{N,i} + \dot y^{-}(t_i) \\
   0\leq  \dot y^{+}(t_i)+ e \dot y^{-}(t_i) \perp P_{N,i} \geq 0

This system is called the "LCP at velocity level" below.

The smooth dynamics
^^^^^^^^^^^^^^^^^^^

The smooth dynamics, valid almost everywhere for the Lebesgue measure :math:`dt`, are governed by the following equation:

.. math::

   M(q) \ddot q^+ +  F_{int}(t, q,  v^+)&= F_{ext}(t) +  f^+ \quad (dt-a.e.)

where we assume that :math:`f^+=f^-=f\, (dt-a.e.)`.

The following smooth systems must then be solved:

.. math::

   M(q(t)) \ddot q^{+}(t) + F_{int}(t, q, v^+)= F_{ext}(t) + f^{+}(t)\\
   y = h(q(t)) \\
   f^+ =  \nabla_q h(q(t))^T F^+(t) \\
   0 \leq y \perp F^+(t) \geq 0

To solve these systems at each time, i.e. to know the configuration after each event and to integrate it numerically, it is useful to write the complementarity laws at different kinematic levels. We also introduce the following index sets:

* :math:`I_0` is the set of all potential UnitaryRelations (UR);
* :math:`I_1 = \{ ur_\alpha\in I_{0} , y_{\alpha} = 0 \}` (if the UR is in :math:`I_1`, there is contact);
* :math:`I_2 = \{ ur_\alpha\in I_{1} , \dot y_{\alpha} = 0 \}` (if the UR is in :math:`I_2`, the contact persists: no take-off).

This leads to the following form of the **bilateral smooth dynamics**:

.. math::

   M(q) \ddot q^{+} + F_{int}(t, q, v)= F_{ext} +  \nabla_q h(q)^T F^+\\ \\
   \ddot y^+ = \nabla_q h(q) \ddot q^+ + \dot{ \nabla_q h(q)} v^+    \\ \\
   F^{+,\alpha} = 0,   \quad \forall \alpha \in I_0-I_2 \\ \\
   \ddot y^{+,\alpha} = 0  \quad \forall \alpha \in I_2

which, if M(q) is invertible, can be reduced to the variables :math:`\ddot y^+` and :math:`F^+` for :math:`\alpha \in I_2`:

.. math::

   \ddot y^{+,\alpha} = \nabla_q h(q) M^{-1}(q)(- F_{int}(t, q, v^+)+ F_{ext}(t)  ) +  \dot{ \nabla_q h(q)} v^+  +\nabla_q h(q) M^{-1}  \nabla_q h(q(t))^T F^{+,\alpha}(t)  \\ \\
   0 \leq \ddot y^{+,\alpha} \perp F^{+,\alpha} \geq 0

This system is called the **"LCP at acceleration level"** below.

The algorithm
-------------

The event-driven algorithm is then as follows.

Knowing :math:`y, \dot y` and :math:`I_1, I_2` at the beginning of the time step :math:`[t_k, t_{k+1}]`:

#. **Integrate the bilateral smooth dynamics** up to an event, given by the root finding of the following function:

   .. math::

      y^\alpha =0,\quad \forall \alpha \in I_0 - I_2 \\
      or \\
      F^{+,\alpha} = 0, \quad \forall \alpha \in I_2

   This gives :math:`y, \dot y` at the new point and updates the index sets :math:`I_1` and :math:`I_2`.

#. If :math:`I_1 - I_2 \neq \emptyset`, impacts occur:

   * formulate and solve the **"LCP at velocity level"**;
   * update the index sets :math:`I_1` and :math:`I_2` and check that :math:`I_1 - I_2 =\emptyset`.

#. If :math:`I_2\neq \emptyset`:

   * formulate and solve the **"LCP at acceleration level"**;
   * for each :math:`\alpha \in I_2`: if :math:`\ddot y_{\alpha} >0, F_{\alpha} = 0`, remove :math:`\alpha` from :math:`I_2` and :math:`I_1`;
     if :math:`\ddot y_{\alpha} =0, F_{\alpha}=0`, the case is undetermined.

#. Go to the next time step.

Implementation in Siconos
^^^^^^^^^^^^^^^^^^^^^^^^^

As explained in :ref:`dynamical_systems`, the dynamics of Lagrangian systems are written in Siconos as:

.. math::
   M(q) \ddot q + fGyr(\dot q, q) + F_{Int}(\dot q , q , t) &= F_{Ext}(t) + p \\

In the following, the :math:`fGyr` term is omitted and considered as included in :math:`F_{Int}`.
The Lagrangian relations are (see :ref:`relations`):

.. math::

   y &= h(Q) \\
   \dot y &= \nabla_q h(Q)\dot Q \\
   P &= \nabla_q h(Q)^t\lambda

Q (resp. P) being the collection of all the q (resp. p) of the Dynamical Systems involved in the Interaction.

As seen in the previous section, the kinematic level is essential. Siconos represents it with the
"[i]" notation: for each Unitary Relation, y[i] is the i-th time derivative of y.
In the same way, :math:`\lambda[i]` is the variable linked to y[i] by a nonsmooth law (usually a complementarity condition).
Finally, each :math:`\lambda[i]` has a corresponding p[i].
To make things clearer, let us rewrite the systems defined above with the Siconos notation:

* **Bilateral smooth dynamics**:

.. math::

   M(q) \ddot q + F_{int}(t, q, \dot q)= F_{ext} +  \nabla_q h(q)^T \lambda[2] \\ \\
   y[2] = \nabla_q h(q) \ddot q + \dot{ \nabla_q h(q)} \dot q    \\ \\
   \lambda[2]_{\alpha} = 0,   \quad \forall \alpha \in I_0-I_2 \\ \\
   y[2]_{\alpha} = 0  \quad \forall \alpha \in I_2

with root finding of:

.. math::

   g(x,t) = y[0]_\alpha,\quad \forall \alpha \in I_0 - I_2 \\
   or \\
   g(x,t) = \lambda[2]_\alpha, \quad \forall \alpha \in I_2

* **"LCP at velocity level"**:

.. math::

   y[1]^{+} =  \nabla_q h(q(t_i)) (M(q(t_i)))^{-1} \nabla_q^T h(q(t_i))\lambda[1] + y[1]^{-} \\
   0\leq y[1]^{+} + e y[1]^{-}  \perp \lambda[1] \geq 0

* **"LCP at acceleration level"**:

.. math::

   y[2]_{\alpha} = \nabla_q h(q) M^{-1}(q)(- F_{int}(t, q, \dot q)+ F_{ext}(t)  ) +  \dot{ \nabla_q h(q)} \dot q  +\nabla_q h(q) M^{-1}  \nabla_q h(q(t))^T \lambda[2]_{\alpha}  \\ \\
   0 \leq y[2]_{\alpha}\perp \lambda[2]_{\alpha} \geq 0

To build an EventDriven simulation, you therefore need two OneStepNSProblems, one at the velocity level and one at the acceleration level.
Here is a typical construction of the simulation::

  EventDriven* s = new EventDriven(ball);
  // -- Time discretisation --
  TimeDiscretisation * t = new TimeDiscretisation(timeStep,s);
  // -- OneStepIntegrators --
  OneStepIntegrator * OSI = new Lsodar(setOfDS,s);
  // -- OneStepNsProblem --
  OneStepNSProblem * impact = new LCP(s, "impact",solverName,101, 0.0001,"max",0.6);
  OneStepNSProblem * acceleration = new LCP(s, "acceleration",solverName,101, 0.0001,"max",0.6);

In Siconos, the algorithm described above becomes:

#. Integration of the bilateral smooth dynamics.
   To integrate these systems with Lsodar, we need to define f(x,t) and g(x,t).
   To compute f(x,t), we:

   * formulate and solve an "LCP at acceleration level" to compute :math:`(y[2],\lambda[2])`;
   * collect and rewrite the dynamics of all the Dynamical Systems as a first-order system, including the result of the LCP.

   The function g(x,t) is given by:

   .. math::

      g(x,t) &= y[0], \quad \forall \alpha \in I_0 - I_2 \\
      \\
      g(x,t) &= \lambda[2], \quad \forall \alpha \in I_2

   Corresponding code::

     s->advanceToEvent()
     // This calls Lsodar->integrate and schedules new nonsmooth events if necessary

   The next steps are done during the call to eventsManager->processEvents(); they are detailed below.

#. Compute y[0] and y[1] and update the index sets::

     simulation->updateOutput(0, 1);
     simulation->updateIndexSets();

#. If :math:`I_1 - I_2 \neq \emptyset`, formulate and solve an LCP at velocity level::

     simulation->computeOneStepNSProblem("impact");

#. Compute p[1], the post-impact velocity, y[1] and indexSet[2]::

     simulation->update(1);

#. If :math:`I_2 \neq \emptyset`, formulate and solve an LCP at acceleration level, and update the index sets with some conditions::

     simulation->computeOneStepNSProblem("acceleration");
     simulation->updateIndexSetsWithDoubleCondition();

#. Next time step::

     simulation->nextStep();
