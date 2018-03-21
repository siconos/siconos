.. _lcp_solvers:

LCP solvers
===========

This page gives an overview of the available solvers for LCP in numerics component and their required parameters.

For each solver, the input argument are:

* a LinearComplementarityProblem structure
* the unknowns z and w
* info, the termination value (0: convergence, >0 problem which depends on the solver)
* a :doxysiconos:`SolverOptions` structure, which handles solver parameters (iparam and dparam)

Remark: when the filterOn parameter (from :doxysiconos:`SolverOptions`) is different from 0, lcp_compute_error() is called at the end of the
process to check the validity of the solution. This function needs a tolerance value and returns an error.
In that case, tolerance is dparam[0] and error output dparam[1]. Thus, in the following solvers, when dparam[0,1] are omitted, that means that they are not required inputs, and that if filter is on, some default values will be used.

lexicographic Lemke
-------------------

Direct solver for LCP based on pivoting method principle for degenerated problem.

 function: :doxysiconos:`lcp_lexicolemke()`
 
 parameters:

 * iparam[0] (in) : max. number of iterations
 * iparam[1] (out): number of iterations processed

QP Solver
----------

quadratic programm formulation for solving a LCP with a symmetric matrix M.

The QP we solve is

  Minimize: :math:`z^T (M z + q)` subject to :math:`Mz  + q  \geq  0`

  which is the classical reformulation that can be found
  in Cottle, Pang and Stone (2009).

  If the symmetry condition is not fulfilled, use the NSQP Solver

function: :doxysiconos:`lcp_qp()`

parameters:

* dparam[0] (in): tolerance

NSQP Solver
-----------

non symmetric (and not nonsmooth as one could have thought in a plateform dedicated to nonsmooth problems)
quadratic programm formulation for solving an LCP with a non symmetric matrix.

function: :doxysiconos:`:doxysiconos:`lcp_nsqp()`

parameters:

* dparam[0] (in): tolerance

CPG Solver
----------

Conjugated Projected Gradient solver for LCP based on quadratic minimization.
Reference: "Conjugate gradient type algorithms for frictional multi-contact problems: applications to granular materials",
M. Renouf, P. Alart. doi:10.1016/j.cma.2004.07.009

function: :doxysiconos:`lcp_cpg()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error

PGS Solver
----------

Projected Gauss-Seidel solver

function: :doxysiconos:`lcp_pgs()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error

RPGS Solver
-----------

Regularized Projected Gauss-Seidel, solver for LCP, able to handle with matrices with null diagonal terms

function: :doxysiconos:`lcp_rpgs()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error
* dparam[2] (in): rho

PSOR Solver
-----------

Projected Succesive over relaxation solver for LCP. See Cottle, Pang and Stone (2009), Chap 5 

function: :doxysiconos:`lcp_psor()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error
* dparam[2] (in): relaxation parameter

NewtonMin Solver
----------------

a nonsmooth Newton method based on the min formulation of the LCP

function: :doxysiconos:`lcp_newton_min()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* iparam[2] (in): if > 0, keep the work vector (reduce the number of memory allocation if the same type of problem is solved multiple times)
* iparam[3] (in): if > 0. use a non-monotone linear search
* iparam[4] (in): if a non-monotone linear search is used, specify the number of merit values to remember
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error

NewtonFB Solver
---------------

a nonsmooth Newton method based based on the Fischer-Burmeister NCP function.
It uses a variant of line search algorithm (VFBLSA in Facchinei-Pang 2003).

function: :doxysiconos:`lcp_newton_FB()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* iparam[2] (in): if > 0, keep the work vector (reduce the number of memory allocation if the same type of problem is solved multiple times)
* iparam[3] (in): if > 0. use a non-monotone linear search
* iparam[4] (in): if a non-monotone linear search is used, specify the number of merit values to remember
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error

Newton min + FB Solver
----------------------

a nonsmooth Newton method based based on the minFBLSA algorithm : the descent direction is given
by a min reformulation but the linesearch is done with Fischer-Burmeister (and if needed the gradient direction).

function: :doxysiconos:`lcp_newton_minFB()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* iparam[2] (in): if > 0, keep the work vector (reduce the number of memory allocation if the same type of problem is solved multiple times)
* iparam[3] (in): if > 0. use a non-monotone linear search
* iparam[4] (in): if a non-monotone linear search is used, specify the number of merit values to remember
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error

Path (Ferris) Solver
--------------------

This solver uses the external PATH solver

function: :doxysiconos:`lcp_path()`

parameters:

* dparam[0] (in): tolerance

Enumerative Solver
------------------

A brute-force method to find the solution of the LCP

function: :doxysiconos:`lcp_enum()`

parameters:

* iparam[0] (in): search for multiple solutions if 1
* iparam[1] (out): key of the solution
* iparam[1] (out): number of solutions
* iparam[3] (in):  starting key values (seed)
* iparam[4] (in):  use DGELS (1) or DGESV (0).
* dparam[0] (in): tolerance

Latin Solver
------------

LArge Time INcrements solver

function: :doxysiconos:`lcp_latin()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error
* dparam[2] (in): latin parameter

Latin_w Solver
--------------

LArge Time INcrements solver with relaxation

function: :doxysiconos:`lcp_latin_w()`

parameters:

* iparam[0] (in): maximum number of iterations allowed
* iparam[1] (out): number of iterations processed
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error
* dparam[2] (in): latin parameter
* dparam[3] (in): relaxation parameter

Block solver (Gauss Seidel)
---------------------------

Gauss-Seidel for Sparse-Block matrices. \n
Matrix M of the LCP must be a SparseBlockStructuredMatrix. \n
This solver first build a local problem for each row of blocks and then call any of the other solvers through lcp_driver()`.

function: :doxysiconos:`lcp_nsgs_SBM()`

parameters:

* iparam[0] (in): maximum number of iterations allowed for GS process
* iparam[1] (out): number of GS iterations processed
* iparam[2] (out): sum of all local number of iterations (if it has sense for the local solver)
* dparam[0] (in): tolerance
* dparam[1] (out): resulting error
* dparam[2] (in): sum of all local error values

