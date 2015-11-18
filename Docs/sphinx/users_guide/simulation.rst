.. _simulation_nsds:

Nonsmooth dynamical systems simulation
======================================

.. toctree::
   :maxdepth: 3

	      
A "simulation" is the description of how a nonsmooth dynamical system will be computed in Siconos.
At the time, two strategies are available:

* event-capturing (time-stepping) schemes
* event-driven schemes
  
Key points of the simulation are :
* the time discretisation
* the time integration of dynamical systems
* the formulation and resolution of one or more nonsmooth problems
  

Time discretisation
-------------------

The discretisation scheme is characterized by a vector of size nSteps+1, denoted tk, that handles the values at each time-step, nSteps being the number of time steps. At the time, only a constant time-step (denoted h) time discretisation is available.

A TimeDiscretisation must be linked to one and only one Model (required parameter of any constructor). This Model provides the initial time value and possibly the final time value, which is an optional parameter. Thus depending on the constructor you use, you may need to give the time-step, the final time ... \n
Just take care to avoid redundant or conflicting information. See the constructors list in the TimeDiscretisation class documentation for full details. 

Example::

  SP::Model m(new Model(t0)); // only initial time is provided
  double h = ...;
  unsigned int nSteps = ...;
  SP::TimeDiscretisation td(new TimeDiscretisation(h, nSteps,m));
  // tk and final time will be automatically computed. 
  // 
  // Then a simulation is created and associated to m through td. 
  SP::Simulation s(new TimeStepping(td));


Time integration of the dynamics
--------------------------------

Dynamical systems integration over a time-step or between two events must be defined thanks to
'one-step integrators'.

* Euler-Moreau (:doxysiconos:`EulerMoreauOSI`)

For first-order dynamical systems, in an 'event-capturing' simulation strategy.
  
.. math::

   M x_{k+1} &=& M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1}) + h(1-\gamma)r(t_k) \\
   y_{k+1} &=&  h(t_{k+1},x_{k+1},\lambda _{k+1}) \\
   r_{k+1} &=& g(x_{k+1},\lambda_{k+1},t_{k+1})\\

with a nonsmooth law linking :math:`y_{k+1}` and :math:`\lambda_{k+1}`,
and :math:`\theta \in [0,1], \gamma \in [0,1]`.

Another variant can also be used (FullThetaGamma scheme)

.. math::

   M x_{k+1} = M x_{k} +h f(x_{k+\theta},t_{k+1}) + h r(t_{k+\gamma}) \\[2mm]
   y_{k+\gamma} =  h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
   r_{k+\gamma} = g(x_{k+\gamma},\lambda_{k+\gamma},t_{k+\gamma})\\[2mm]
   \mbox{nslaw} ( y_{k+\gamma} , \lambda_{k+\gamma})

* Moreau-Jean (:doxysiconos:`MoreauJeanOSI`)

For mechanical (second-order) systems, in an 'event-capturing' simulation strategy.

.. math::
   
   M (v_{k+1}-v_k) + h K q_{k+\theta} + h C v_{k+\theta} - h F_{k+\theta} = p_{k+1} = G P_{k+1}\\ 
   q_{k+1} = q_{k} + h v_{k+\theta}, \\
   U_{k+1} = G^\top\, v_{k+1}, \\
   \begin{array}{lcl}
   0 \leq U^\alpha_{k+1} + e  U^\alpha_{k} \perp P^\alpha_{k+1}  \geq 0,& \quad&\alpha \in \mathcal I_1, \\
   P^\alpha_{k+1}  =0,&\quad& \alpha \in \mathcal I \setminus \mathcal I_1,\end{array}

with  :math:`\theta \in [0,1]`. The index set :math:`\mathcal I_1` is the discrete equivalent
to the rule that allows us to apply the Signorini  condition at the velocity level.
Numerically, this set is defined as

.. math::

   \mathcal I_1 = \{\alpha \in \mathcal I \mid G^\top (q_{k} + h v_{k}) + w \leq 0\text{ and } U_k \leq 0 \}.

* Schatzman-Paoli (:doxysiconos:`SchatzmanPaoliOSI`)
* zero-order  (:doxysiconos:`ZeroOrderHoldOSI`) 
* Lsodar (:doxysiconos:`LsodarOSI`)

For 'event-driven' simulation strategy. Integrator based on LSODAR (https://computation.llnl.gov/casc/odepack/) rootfinding routine :
"Lsodar solves problems dy/dt = f with full or banded Jacobian and automatic method selection, and at the same time, it finds the roots of any of a set of given functions of the form g(t,y). This is often useful for finding stop conditions or points at which switches are to be made in the function f". 
In Siconos, Lsodar is used for event-driven algorithm, to integrate the dynamics with stops at new non-smooth events (violation of a constraint)

* Hem5 (:doxysiconos:`Hem5OSI`)
For 'event-driven' simulation strategy. Based on Ernst Hairer HEM5 integrator (http://www.unige.ch/~hairer/software.html)

* Newmark (:doxysiconos:`NewmarkAlphaOSI`)



Nonsmooth problems formulation and solve
----------------------------------------

When dynamical systems and their interactions have been properly defined inside a model and its nonsmooth dynamical system,
a proper formulation for the latter must be chosen, associated to a nonsmooth solver.

Linear nonsmooth problems
^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n}
where :math:`w \in R^{n}, z \in R^{n}` are the unknowns.

* Linear Complementarity Problems (:doxysiconos:`LCP`)

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n} \\
   w \geq 0, z \geq 0,  z^{T} w =0

* Mixed Linear Complementarity Problems (:doxysiconos:`MLCP`)

.. math::
  
   0 =  Au + Cv + a\\
   z =  Du + Bv + b\\
   v \geq 0, z \geq 0,  z^{T} v =0

where

.. math::
   u \in R^{n},  v \in R^{m}, z \in R^{m} \ the \ unknowns, \\
   a \in R^{n}, b \in R^{m} \\
   A \in R^{n \times n }, B \in R^{m \times m }\\
   C \in R^{n \times m }, D \in R^{m \times n }
   
* :doxysiconos:`Relay`
* :doxysiconos:`Equality`
* :doxysiconos:`AVI`
* 2D or 3D friction contact problem :doxysiconos:`FrictionContact`

.. math::

  velocity =  q + M reaction \\
  velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0

and a Coulomb friction law.
With :math:`velocity \in R^{n}, reaction \in R^{n}` the unknowns,
and :math:`M \in R^{n \times n }, q \in R^{n}`

  
* multiple-impact problem (:doxysiconos:`OSNSMultipleImpact`)
* primal friction contact problems (:doxysiconos:`GlobalFrictionContact`)

.. math::

   M velocity =  q +  reaction \\
   localVelocity = H^T velocity + tildeLocalVelocity\\
   reaction = H localReaction \\

and :math:`localVelocity,localReaction` belongs to the Coulomb friction law with unilateral contact.

With :math:`velocity \in R^{n}, reaction \in R^{n}, localVelocity \in R^{m}, localReaction \in R^{m}` the unknowns,
:math:`M \in R^{n \times n }, q \in R^{n}`. :math:`tildeLocalVelocity \in R^{m}` is the modified local velocity (:math:`e U_{N,k}`), :math:`M \in R^{n \times n }, q \in R^{n}, H \in R^{n \times m }`.
   
* Generic mechanical problem (:doxysiconos:`GenericMechanical`)
  
Complete problem with bilateral equality, complementarity, impact and friction.
