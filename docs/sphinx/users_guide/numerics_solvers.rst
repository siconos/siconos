.. _numerics_solver:

Solvers definition (numerics)
=============================

To define a non-smooth problem in Numerics, the structure :class:`SolverOptions` is used. It handles the name of the solver and its input-output parameters.

:class:`SolverOptions` main components are:
   * a name
   * two lists of input-output parameters (int: iparam, double: dparam) and their sizes

Check each type of formulation of the problem to find which solvers are available and what are the required parameters,
see for example :ref:`lcp_solvers` or :ref:`fc_solvers`.

As an example, consider a Linear Complementarity Problem :
M is a NumericsMatrix and can be saved as a double* or as a SparseBlockStructuredMatrix.
One needs to define a SolverOptions, say "options", by choosing one solver among those given in :ref:`lcp_solvers` and set::

  int nbSolvers = 1;
  SolverOptions options;
  strcpy(options.solverName,"PGS");
  int iparam[2] ={maxIter, 0};
  double dparam[2] = {tolerance,0.0};
  options.iSize = 2;
  options.dSize = 2;
  options.iparam = iparam;
  options.dparam = dparam;
  options.isSet = 1;

And then call the driver::

  int info = lcp_driver(myProblem, z,w, &options, nbSolvers, &global_options);

which will result in the resolution of the LCP defined in myProblem thanks to a PGS solver.

On the other side if M is saved as a SparseBlockStructuredMatrix, with N rows of blocks, one needs to used a
"block-solver" with possibly one or more specific local solver dedicated to each local problem.
In that case options must be a vector of SolverOptions, with:

* options[0] the definition for the global "block" solver
* options[i], i>0, the solvers used for each local problem.

Example with a LCP::

  // First define a vector of options
  int nbSolvers = 3;
  SolverOptions options[nbSolvers];

  // The global solver:
  strcpy(options[0].solverName,"GaussSeidel_SBM");
  int iparam[2] ={maxIter, 0};
  double dparam[2] = {tolerance,0.0};
  options[0].iSize = 2;
  options[0].dSize = 2;
  options[0].iparam = iparam;
  options[0].dparam = dparam;
  options[0].isSet = 1;
  
  // The local solvers:
  strcpy(options[1].solverName,"PGS");
  int iparam[2] ={maxIter, 0};
  double dparam[2] = {tolerance,0.0};
  options[1].iSize = 2;
  options[1].dSize = 2;
  options[1].iparam = iparam;
  options[1].dparam = dparam;
  options[1].isSet = 1;
  strcpy(options[2].solverName,"Lemke");
  int iparam[2] ={maxIter,0};
  double dparam[2] = {tolerance,0.0};
  options[2].iSize = 2;
  options[2].dSize = 2;
  options[2].iparam = iparam;
  options[2].dparam = dparam;
  options[2].isSet = 1;

The call of the driver remains the same::

  int info = lcp_driver(myProblem, z,w, options,nbSolvers, &global_options);

In this case, if the matrix M has N rows of blocks, the global problem will be solved thanks to the Gauss-Seidel block solver, with the first local problem (first row) solved thanks to a PGS and the others with a Lemke.
Note that options[i+1] is used for row i of M, while i<nbSolvers-1 and options[nbSolvers-1] for row i when i>=nbSolvers.


.. toctree::
   :maxdepth: 4

   lcp_solvers
   fc_solvers
