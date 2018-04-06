.. _fc_solvers:

Friction-Contact solvers
========================

This page gives an overview of the available solvers for friction-contact (2D and 3D) problems and their required parameters.

For each solver, the input argument are:

* a FrictionContactProblem
* the unknowns (reaction,velocity)
* info, the termination value (0: convergence, >0 problem which depends on the solver)
* a :class:`SolverOptions` structure, which handles iparam and dparam

2D solvers
----------

CPG
^^^

function: :function:`pfc_2D_cpg()`

parameters:

* iparam[0] (in), the maximum number of iterations allowed,
* iparam[1] (out), the number of iterations performed by the algorithm.
* dparam[0] (in), the tolerance required,
* dparam[1] (out), the residu.

3D solvers
----------

Non-Smooth Gauss Seidel
^^^^^^^^^^^^^^^^^^^^^^^

function: :function:`fc3d_nsgs()`

parameters:



