.. _osns_problems:

Nonsmooth problems formulation and solve
========================================

When dynamical systems and their interactions have been properly defined inside a model and its nonsmooth dynamical system,
a proper formulation for the latter must be chosen, associated to a nonsmooth solver.

For details regarding the available formulations, see :ref:`problems_and_solvers`.


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

Each time :func:`ComputeOneStepNS()` function, i.e. the numerics solver, is called, it returns an int, giving some information about the convergence of the solver:

* output = 0 => solver succeeded,
* else, the meaning of output depends on the solver called (see :ref:`problems_and_solvers`).
  
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

