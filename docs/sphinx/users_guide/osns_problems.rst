.. _osns_problems:

Nonsmooth problems formulation and solve
========================================

When dynamical systems and their interactions have been properly defined inside a model and its nonsmooth dynamical system,
a proper formulation for the latter must be chosen, associated to a nonsmooth solver.

Linear nonsmooth problems
-------------------------

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n}

where :math:`w \in R^{n}, z \in R^{n}` are the unknowns.

* Linear Complementarity Problems (:class:`LCP`)

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n} \\
   w \geq 0, z \geq 0,  z^{T} w =0

* Mixed Linear Complementarity Problems (:class:`MLCP`)

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
   
* :class:`Relay`
* :class:`Equality`
* :class:`AVI`
* 2D or 3D friction contact problem :class:`FrictionContact`

.. math::

  velocity =  q + M reaction \\
  velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0

and a Coulomb friction law.
With :math:`velocity \in R^{n}, reaction \in R^{n}` the unknowns,
and :math:`M \in R^{n \times n }, q \in R^{n}`

  
* multiple-impact problem (:class:`OSNSMultipleImpact`)
* primal friction contact problems (:class:`GlobalFrictionContact`)

.. math::

   M velocity =  q +  reaction \\
   localVelocity = H^T velocity + tildeLocalVelocity\\
   reaction = H localReaction \\

and :math:`localVelocity,localReaction` belongs to the Coulomb friction law with unilateral contact.

With :math:`velocity \in R^{n}, reaction \in R^{n}, localVelocity \in R^{m}, localReaction \in R^{m}` the unknowns,
:math:`M \in R^{n \times n }, q \in R^{n}`. :math:`tildeLocalVelocity \in R^{m}` is the modified local velocity (:math:`e U_{N,k}`), :math:`M \in R^{n \times n }, q \in R^{n}, H \in R^{n \times m }`.
   
* Generic mechanical problem (:class:`GenericMechanical`)
  
Complete problem with bilateral equality, complementarity, impact and friction.


The Simulation process
----------------------

As for Event-Driven, we introduce level index sets, with level = 0 for first order systems and level=1 for second order systems (this is related to the relative degrees but we won't get into details about that here).

:math:`I_0` is the set of all the potential UnitaryRelations (UR).
For second order systems:
:math:`I_1 = \{ ur_\alpha\in I_{0} , y^p_{\alpha} = 0 \}`. 
Thus, the LCP is built only for unitary relations that belongs to :math:`I_level`, level=0 for first order and level=1 for second order systems. 

Then, the steps of a Moreau's Time-Stepping simulation will be:

Knowing all values at the beginning of the time step :math:`[t_i,t_{i+1}]`,

-# compute the free solutions
-# for :math:`ur \in I_level` formalize and solve a LCP 
-# update the state (according to the possibly LCP results)
-# go to next time step

::
   
   SP::TimeStepping s(new TimeStepping(myModel));  
   SP::TimeDiscretisation t(new TimeDiscretisation(timeStep,s));

   s->initialize();

   int N = t->getNSteps(); // Number of time steps
   
   // --- Time loop ---
   while(k < N)// for each time step ...
   {
   // compute xFree, or qFree,vFree
   s->computeFreeStep();
   // Formalize and solve a LCP
   computeOneStepNSProblem("timeStepping");
   // Update state, using last computed values
   s->update(level); // 
   // transfer of state i+1 into state i and time incrementation
   s->nextStep();
   }

Note that all time-independent operators are computed during simulation initialisation.



Customize simulation behavior
-----------------------------

Each time :function:`ComputeOneStepNS()` function, i.e. the numerics solver, is called, it returns an int, giving some information about the convergence of the solver:

* output = 0 => solver succeeded,
* else, the meaning of output depends on the solver called (see :ref:`numerics_solvers`).
  
By default, when the convergence is not achieved, an exception is throwed and the process stops.
Change this behavior is possible by defining a specific function of the form::

  //
  // your inputFile.cpp
  //
  void myF(int info, SP::Simulation s)
  {
  // do what you need ...
  }
  
  int main(int argc, char* argv[])
  {
  // 
  // ...
  SP::TimeStepping your_simulation = ...
  your_simulation->setCheckSolverFunction(&myF);

Then after each call to your_simulation->computeOneStepNS(...), the function myF will be called.
That may be usefull to change the solver type, the tolerance or whatever is needed. 

