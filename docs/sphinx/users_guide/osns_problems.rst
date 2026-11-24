.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _osns_problems:

Formulation and solution of nonsmooth problems
==============================================

Once the dynamical systems and their interactions have been defined in a model and its nonsmooth dynamical system,
a formulation of the latter must be chosen, together with a nonsmooth solver.

For details regarding the available formulations, see :ref:`problems_and_solvers`.


The simulation process
----------------------

As for event-driven simulations, we introduce index sets for each level, with level = 0 for first-order systems and level = 1 for second-order systems (this is related to the relative degree, but we will not go into details here).

:math:`I_0` is the set of all potential UnitaryRelations (UR).
For second-order systems:
:math:`I_1 = \{ ur_\alpha\in I_{0} , y^p_{\alpha} = 0 \}`.
The LCP is therefore built only for the unitary relations that belong to :math:`I_{level}`, with level = 0 for first-order and level = 1 for second-order systems.

The steps of a Moreau time-stepping simulation are then:

Knowing all values at the beginning of the time step :math:`[t_i,t_{i+1}]`,

#. compute the free solutions;
#. for :math:`ur \in I_{level}`, formulate and solve an LCP;
#. update the state (using the LCP results, if any);
#. go to the next time step.

::

   auto s = std::make_shared<siconos::simulation::TimeStepping>(myModel);
   auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(timeStep,s);

   s->initialize();

   int N = t->getNSteps(); // Number of time steps

   // --- Time loop ---
   while(k < N)// for each time step ...
   {
   // compute xFree, or qFree,vFree
   s->computeFreeStep();
   // Formulate and solve an LCP
   computeOneStepNSProblem("timeStepping");
   // Update state, using last computed values
   s->update(level); //
   // transfer of state i+1 into state i and time incrementation
   s->nextStep();
   }

Note that all time-independent operators are computed during the initialisation of the simulation.



Customize the simulation behaviour
----------------------------------

Each call to the :cpp:func:`ComputeOneStepNS()` function, i.e. to the numerics solver, returns an int that gives information about the convergence of the solver:

* output = 0 => the solver succeeded,
* otherwise, the meaning of the output depends on the solver (see :ref:`problems_and_solvers`).

By default, when the solver does not converge, an exception is thrown and the process stops.
This behaviour can be changed by defining a function of the form::

  //
  // your inputFile.cpp
  //
  void myF(int info, std::shared_ptr<siconos::simulation::Simulation> s)
  {
  // do what you need ...
  }

  int main(int argc, char* argv[])
  {
  //
  // ...
  auto your_simulation = ...
  your_simulation->setCheckSolverFunction(&myF);

Then, after each call to your_simulation->computeOneStepNS(...), the function myF is called.
This is useful to change the solver type, the tolerance or anything else.
