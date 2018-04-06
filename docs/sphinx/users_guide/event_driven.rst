.. _event_driven:

Event-Driven schemes
====================

General Principle
-----------------

The principle of the event-driven is based on the time-decomposition of the dynamics in modes, time-intervals where the dynamics is smooth, and discrete events, times where the dynamics are nonsmooth. From the numerical point of view, the event-driven scheme use the decomposition in time of the dynamics in order to

* detect and solve the non smooth dynamics at events with a reinitialization rule of the state,
* integrate the smooth dynamics between two events with any ODE solvers with root-findings and possibly bilateral constraints on the state.

Event Driven implementation
---------------------------

Integration of the smooth dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Between events, the dynamics are integrated thanks to Lsodar algorithm and with the function EventDriven::integrate.

Considering two given functions f(x,t) and g(x,t), a call to::

  integrate(tinit, tend, tout, iout)

results in the integration of the function f(x,t) between tinit and tend, and search for roots of the function g(x,t). If roots are found, integration stops, and last time is saved in tout.

The in-out parameter iout is an indicator that must be set to 1 at first call. If no root was found, it is equal to 2 if so to 3.

Thus for an EventDriven simulation:

* the dynamics of all the concerned DynamicalSystems is rewritten as :math:`\dot = f(x,t)`
* the relations are used to defined the g(x,t) functions

Once again the proper definition of f and g depends on the system type as described below.

Events
^^^^^^

In Siconos the object Event just handle a type (see below) and a long int which corresponds to the time of occurrence of the event. A process function is also defined, which action depends on the event type.

The possible types (derived classes) are:

* TimeDiscretisationEvent: event that corresponds to user-defined time discretisation points. These events are created and schedule when the EventDriven simulation is initialized. Process function call results in :
  * update (compute) output for all the concerned Interactions
  * save current values (DynamicalSystems states and Interactions input/output) in memory vectors. Last saved values become initial values for next integration.

* NonSmoothEvent: "points" where the dynamics are non smooth and which required a special treatment. These events are detected thanks to a roots-finding algorithm, and corresponds to violation of some given constraints (relation). The action of the process function is roughly (the full process depends on the system type and is described in \ref docSimuEDDetails):
  * update (compute) output for all the concerned Interactions
  * update the index sets 
  * formalize and solve one or more non-smooth problems
  * save current values (DynamicalSystems states and Interactions input/output) in memory vectors. Last saved values become initial values for next integration.

The Events manager
^^^^^^^^^^^^^^^^^^

To handle all the events, a specific object is built: the EventsManager. It belongs to the EventDriven class and holds two eventsContainers (sets of Events):
- pastEvents: for all the Events that has already been treated
- unProcessedEvents: for the future events, already scheduled but not treated
We also denote "currentEvent" the last processed event, which corresponds to the initial point of the current integration, and "nextEvent" the event following "currentEvent". 

The events manager is initialized with time-discretisation events, from the user time-discretisation. Then, each time a new event is detected (added by user or when a root is found during integration) it is scheduled in the manager.
The manager has also a \link EventsManager::processEvents processEvents \endlink function, which moves currentEvent to past events set, processes nextEvent and prepare the next step.

Other useful functions are: 
- \link EventsManager::startingTime() startingTime \endlink, \link EventsManager::nextTime() nextTime \endlink

The Simulation process
----------------------

Thus, a general step  of integration for EventDriven will looks like:
"While there are some events in the unProcessedEvents set, integrate the smooth dynamics between current and next event and then process with the behavior at event".
Or:

code:: c++
  
  SP::EventDriven s(new EventDriven(myModel));

  s->initialize();

  // We get the events manager
  SP::EventsManager eventsManager = s->eventsManager();
  
  // while there are some events ...
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

Event Driven algorithm for Lagrangian systems
---------------------------------------------

At the time, the only available event-driven algorithm in Siconos is for Lagrangian dynamical systems, subjected to perfect unilateral constraints and with the Newton impact rules.

Because of the unilateral constraints, the evolution of the considered system may be non-smooth. Some jumps can occur in the velocity and the "acceleration" may not be defined everywhere. The generalized coordinates, assumed to be absolutely continuous are:

.. math::
   
   q(t) = q(t_0) +\int_{t_0}^t v^+(t)dt \ with \ v = \dot q

We will index with "+" and "-" right and left values of the variable at discontinuity.

The equations of motion are written in terms of a measure differential equation:

.. math::
   
   M(q)dv +  F_{int}(t, q,  v^+)dt=F_{ext}(t) +  dr

r being the generalized force due the unilateral constraints. 
Using the Lebesgue decomposition theorem and its variants, the differential measure dv and dr are decomposed in:

.. math::
   
   dv = \gamma dt + (v^+-v^-)\sum_i\delta_{t_i} + dv_s \\
   dr = fdt + \sum_ip_i\delta_{t_i}+dr_s

First term of the decomposition corresponds to the smooth part, with :math:`\gamma =\ddot q`, the acceleration in the usual sense. The second term corresponds to the behavior at times of discontinuities, ( :math:`\delta_{t_i}`: Dirac), and the last term, a singular measure, will be neglected.

Thanks to these decompositions, the non-smooth Dynamics can be split into "impact equations", that will correspond to the non-smooth events, and some "smooth Dynamics". These equations are completed by the constraints, formulated at different kinematics levels, as shown in the following paragraphs.

The impact equations
^^^^^^^^^^^^^^^^^^^^

The impact equations can be written at the time :math:`t_i` of discontinuities:

.. math::
   
   M(q(t_i))(v^{+}(t_i)- v^{-}(t_i)) = p_i,
   
:math:`p_i` is like an impulsion.
      
This equation will be solved at the time of impact together with an impact law. That is for a Newton impact law

.. math::
   M(q(t_i))(v^{+}(t_i)- v^{-}(t_i)) = p_i, \\
   \dot y^{+}(t_i) = \nabla_q h(q(t_i)) v^{+}(t_i) \\
   \dot y^{-}(t_i) = \nabla_q h(q(t_i)) v^{-}(t_i) \\
   p_i =   \nabla_q^T h(q(t_i)) P_{N,i}\\
   0\leq  \dot y^{+}(t_i)+ e \dot y^{-}(t_i) \perp P_{N,i} \geq 0

This problem can be reduced on the local unknowns :math:`\dot y^{+}(t_i),P_{N,i}` if the matrix :math:`M(q(t_i))` is assumed to be invertible, leading to the following Linear Complementarity Problem at time :math:`t_i` of discontinuities of v:

.. math::
   \dot y^{+}(t_i) =  \nabla_q h(q(t_i)) (M(q(t_i)))^{-1} \nabla_q^T h(q(t_i))   P_{N,i} + \dot y^{-}(t_i) \\ 
   0\leq  \dot y^{+}(t_i)+ e \dot y^{-}(t_i) \perp P_{N,i} \geq 0

Later this system will be identified as "LCP at velocity level". 

The smooth Dynamics
^^^^^^^^^^^^^^^^^^^

The smooth dynamics which is valid almost everywhere for the Lebesgue measure :math:`dt` is governed by  the following equation:

.. math::

   M(q) \ddot q^+ +  F_{int}(t, q,  v^+)&= F_{ext}(t) +  f^+ \quad (dt-a.e.)

where we assume that :math:`f^+=f^-=f\, (dt-a.e.)`.

The following smooth systems are then to be solved:

.. math::
   
   M(q(t)) \ddot q^{+}(t) + F_{int}(t, q, v^+)= F_{ext}(t) + f^{+}(t)\\
   y = h(q(t)) \\
   f^+ =  \nabla_q h(q(t))^T F^+(t) \\
   0 \leq y \perp F^+(t) \geq 0

To solve these systems, at each time, i.e. to known the configuration after each events and to integrate it numerically, it is useful to express the complementarity laws at different kinematics level. We also introduce the pre-defined index sets (about index sets, see \ref docSimuIndexSets):\n

:math:`I_0` is the set of all the potential UnitaryRelations (UR).
:math:`I_1 = \{ ur_\alpha\in I_{0} , y_{\alpha} = 0 \}` (or if the UR is in :math:`I_1` then contact occurs).
:math:`I_2 = \{ ur_\alpha\in I_{1} , \dot y_{\alpha} = 0 \}` (or if the UR is in :math:`I_2`, contact remains, no take off).

This results in the new writing of the <b>Bilateral Smooth Dynamics</b>: 

.. math::

   M(q) \ddot q^{+} + F_{int}(t, q, v)= F_{ext} +  \nabla_q h(q)^T F^+\\ \\
   \ddot y^+ = \nabla_q h(q) \ddot q^+ + \dot{ \nabla_q h(q)} v^+    \\ \\
   F^{+,\alpha} = 0,   \quad \forall \alpha \in I_0-I_2 \\ \\
   \ddot y^{+,\alpha} = 0  \quad \forall \alpha \in I_2

which can be reduced on variable :math:`\ddot y^+` and :math:`F^+`, if M(q) is invertible, when :math:`\alpha \in I_2`:

.. math::

   \ddot y^{+,\alpha} = \nabla_q h(q) M^{-1}(q)(- F_{int}(t, q, v^+)+ F_{ext}(t)  ) +  \dot{ \nabla_q h(q)} v^+  +\nabla_q h(q) M^{-1}  \nabla_q h(q(t))^T F^{+,\alpha}(t)  \\ \\
   0 \leq \ddot y^{+,\alpha} \perp F^{+,\alpha} \geq 0 

Later this system will be identified as <b>"LCP at acceleration level"</b>. 

The algorithm
-------------

Finally, the event-driven algorithm will be:

knowing the value of :math:`y, \dot y` and :math:`I_1, I_2` at the beginning of the time step :math:`[t_k, t_{k+1}]`:

-# <b> Integration of the Bilateral Smooth Dynamics </b> up to an event given by the root-finding of the following function :

.. math::
   y^\alpha =0,\quad \forall \alpha \in I_0 - I_2 \\
   or \\
   F^{+,\alpha} = 0, \quad \forall \alpha \in I_2

This results in the computation of :math:`y, \dot y` at this new point and to an update of the index sets :math:`I_1` and :math:`I_2`.

-# if :math:`I_1 - I_2 \neq \emptyset` then Impacts occur: 
    - Formalize and solve the <b>"LCP at velocity level"</b>
    - Update the index sets :math:`I_1` and :math:`I_2` and check that  :math:`I_1 - I_2 =\emptyset`
   endif

-# if :math:`I_2\neq \emptyset` then 
    - Formalize and solve the <b>"LCP at acceleration level"</b>
    - for :math:`\alpha \in I_2` do
      if :math:`\ddot y_{\alpha} >0, F_{\alpha} = 0` remove :math:`\alpha` from :math:`I_2` and :math:`I_1`
      else if :math:`\ddot y_{\alpha} =0, F_{\alpha}=0` then undetermined case.
      endif\n
     endfor\n 
    endif\n

-# go to the next time step.

Implementation in Siconos
^^^^^^^^^^^^^^^^^^^^^^^^^

According to \ref doc_lagds, in Siconos, the Dynamics of Lagrangian systems is written as:

.. math::
   M(q) \ddot q + fGyr(\dot q, q) + F_{Int}(\dot q , q , t) &= F_{Ext}(t) + p \\

Next,:math:`fGyr` term will be forget and considered as included in :math:`F_{Int}`.
And Lagrangian relations are (see \ref docRelationLag): 

.. math::

   y &= h(Q) \\
   \dot y &= \nabla_q h(Q)\dot Q \\
   P &= \nabla_q h(Q)^t\lambda 

Q (resp. P) being a collection of all the q (resp. p) of the Dynamical Systems involved in the Interaction.

As we have seen in the previous section, the notion of kinematics level is really important. We introduce this in Siconos thanks to 
"[i]" notation. More precisely, for each Unitary Relation, we define y[i] as the derivative number i of variable y, according to time.
In the same way, we denote :math:`\lambda[i]` the variable that is linked with y[i] through a Non-Smooth law (usually a complementarity). 
Finally to each :math:`\lambda[i]` corresponds a p[i].
To make things clearer, let us rewrite the previous defined systems with Siconos notations: 

- <b>Bilateral Smooth Dynamics</b>:

.. math::

   M(q) \ddot q + F_{int}(t, q, \dot q)= F_{ext} +  \nabla_q h(q)^T \lambda[2] \\ \\
   y[2] = \nabla_q h(q) \ddot q + \dot{ \nabla_q h(q)} \dot q    \\ \\
   \lambda[2]_{\alpha} = 0,   \quad \forall \alpha \in I_0-I_2 \\ \\
   y[2]_{\alpha} = 0  \quad \forall \alpha \in I_2

with roots finding of:

.. math::
   
   g(x,t) = y[0]_\alpha,\quad \forall \alpha \in I_0 - I_2 \\
   or \\
   g(x,t) = \lambda[2]_\alpha, \quad \forall \alpha \in I_2

- <b>"LCP at velocity level"</b>

.. math::

   y[1]^{+} =  \nabla_q h(q(t_i)) (M(q(t_i)))^{-1} \nabla_q^T h(q(t_i))\lambda[1] + y[1]^{-} \\ 
   0\leq y[1]^{+} + e y[1]^{-}  \perp \lambda[1] \geq 0

- <b>"LCP at acceleration level"</b>

.. math::
   
   y[2]_{\alpha} = \nabla_q h(q) M^{-1}(q)(- F_{int}(t, q, \dot q)+ F_{ext}(t)  ) +  \dot{ \nabla_q h(q)} \dot q  +\nabla_q h(q) M^{-1}  \nabla_q h(q(t))^T \lambda[2]_{\alpha}  \\ \\
   0 \leq y[2]_{\alpha}\perp \lambda[2]_{\alpha} \geq 0 

Then, to build an EventDriven simulation, it is necessary to define two OneStepNSProblems, one at velocity and one at acceleration level.
So here is a classical code for simulation construction::

  EventDriven* s = new EventDriven(ball);
  // -- Time discretisation --
  TimeDiscretisation * t = new TimeDiscretisation(timeStep,s);
  // -- OneStepIntegrators --
  OneStepIntegrator * OSI = new Lsodar(setOfDS,s); 
  // -- OneStepNsProblem --
  OneStepNSProblem * impact = new LCP(s, "impact",solverName,101, 0.0001,"max",0.6);
  OneStepNSProblem * acceleration = new LCP(s, "acceleration",solverName,101, 0.0001,"max",0.6);

Finally, the algorithm described earlier is:

-# Integration of the Bilateral Smooth Dynamics:
To integrate these systems thanks to lsodar, we need to define f(x,t) and g(x,t).
To compute f(x,t), we:
  - formalize and solve a "LCP at acceleration level" to compute :math:`(y[2],lambda[2])`
  - collect and rewrite the Dynamics of all the Dynamical Systems as a first order system, including the result of the LCP computation.
The function g(x,t) is given by:

.. math::

   g(x,t) &= y[0], \quad \forall \alpha \in I_0 - I_2 \\
   \\
   g(x,t) &= \lambda[2], \quad \forall \alpha \in I_2
   
Corresponding code::

  s->advanceToEvent()
  // This results in a call to Lsodar->integrate and to schedule of new non-smooth events if necessary
  
The next steps are done during call to eventsManager->processEvents(), but they will be detailed below.
-# Compute y[0] and y[1] and update the index sets::

  simulation->updateOutput(0, 1);
  simulation->updateIndexSets();

-# if :math:`I_1 - I_2 \neq \emptyset`, formalize and solve a LCP at velocity level::

  simulation->computeOneStepNSProblem("impact"); 

-# compute p[1], post-impact velocity, y[1] and indexSet[2]::

  simulation->update(1);
  
-# if :math:`I_2 \neq \emptyset`, formalize and solve a LCP at acceleration level, and update index sets with some conditions::

  simulation->computeOneStepNSProblem("acceleration");
  simulation->updateIndexSetsWithDoubleCondition();

-# next time step::

  simulation->nextStep();
