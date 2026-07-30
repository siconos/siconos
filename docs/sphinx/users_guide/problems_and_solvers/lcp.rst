.. index::
   single: Linear Complementarity Problem (LCP)
   
.. contents::

.. _lcp_problem:

Linear Complementarity Problems (LCP)
*************************************

Problem statement
=================

The Linear Complementarity problem (LCP) is defined by

Find :math:`(z,w)` such that:

.. math::

    \begin{equation*} \begin{cases}
    M \ z + q = w \\
    0 \le w \perp z \ge 0
    \end{cases} \end{equation*}


For more details on theory and analysis of LCP, we refer to :cite:`Cottle.1992`.

Implementation in numerics
==========================

Structure to define the problem: :cpp:class:`LinearComplementarityProblem`.

The generic driver for all LCPs is :cpp:func:`linearComplementarity_driver()`
which switches either to :cpp:func:`lcp_driver_DenseMatrix` or :cpp:func:`lcp_driver_SparseMatrix` according
to the kind of storage used in the problem.

Solvers list  :cpp:enum:`LCP_SOLVER`

.. _lcp_error:

Error computation
=================

The criterion is based on :

.. math::

   error = \frac{1}{\|q\| }\sqrt(\sum_{i} (z_i - (z_i - (Mz+q)_i)^{+})^2)

   error = \frac{1}{\|q\| }\sum_{i} [ (z_i*(Mz+q)_i)^{+} + (z_i)^{-} + {(Mz+q)_i}^{-} ]
   
with :math:`x^{+} = max(0,x)` and :math:`x^{-} = max(0,-x)`.



* :cpp:func:`lcp_compute_error` returns 0  if  :math:`error \leq tolerance`, else 1.
* A call to this function updates the content of the input vector w with :math:`Mz + q`.

  
.. _lcp_solvers:

LCP available solvers
=====================

Direct solvers
--------------

Lemke (:cpp:enumerator:`SICONOS_LCP_LEMKE`)
"""""""""""""""""""""""""""""""""""""""""""

direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

driver :cpp:func:`lcp_lexicolemke`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[2] = 0.0 
* dparam[3] = 0.0

Pivot based methods
"""""""""""""""""""

:cpp:enumerator:`SICONOS_LCP_PIVOT`, :cpp:enumerator:`SICONOS_LCP_BARD`,
:cpp:enumerator:`SICONOS_LCP_MURTY`, :cpp:enumerator:`SICONOS_LCP_PATHSEARCH`,
:cpp:enumerator:`SICONOS_LCP_PIVOT_LUMOD`
            
generic solver for pivot-based methods: Bard, Murty and Lemke rules are implemented.

drivers:

* :cpp:func:`lcp_pivot` for :cpp:enumerator:`SICONOS_LCP_PIVOT`, :cpp:enumerator:`SICONOS_LCP_BARD`, and :cpp:enumerator:`SICONOS_LCP_MURTY`,
* :cpp:func:`lcp_pathsearch` for :cpp:enumerator:`SICONOS_LCP_PATHSEARCH`,
* :cpp:func:`lcp_pivot_lumod` for :cpp:enumerator:`SICONOS_LCP_PIVOT_LUMOD`.


parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] =

  * SICONOS_LCP_PIVOT_BARD for SICONOS_LCP_BARD
  * SICONOS_LCP_PIVOT_LEAST_INDEX for SICONOS_LCP_MURTY
  * SICONOS_LCP_PIVOT_LEMKE for SICONOS_LCP_PIVOT, SICONOS_LCP_PATHSEARCH and SICONOS_LCP_PIVOT_LUMOD

* iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = 0;
* dparam[SICONOS_DPARAM_TOL] = 100 * epsilon (machine precision)



Enumerative solver (:cpp:enumerator:`SICONOS_LCP_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

Brute-force method which tries every possible solution.

driver: :cpp:func:`lcp_enum()`

parameters:

* iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 0 (0 : use dgesv, 1: use dgels)
* iparam[SICONOS_LCP_IPARAM_ENUM_CURRENT_ENUM] (out): key of the solution
* iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = 0, starting key values
* iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS] = 0,  search for multiple solutions if 1
* iparam[SICONOS_LCP_IPARAM_ENUM_NUMBER_OF_SOLUTIONS] (out): number of solutions
* dparam[SICONOS_DPARAM_TOL] = 1e-6

PATH (:cpp:enumerator:`SICONOS_LCP_PATH`)
"""""""""""""""""""""""""""""""""""""""""

*Works only if Siconos has been built with path support (if PathFerris or PathVI has been found, see :ref:`siconos_install_guide`)*

driver: :cpp:func:`lcp_path()`

parameters
* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

 
Iterative solvers
-----------------

Conjugated Projected Gradient (:cpp:enumerator:`SICONOS_LCP_CPG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Conjugated Projected Gradient solver for LCP based on quadratic minimization.
Reference: "Conjugate gradient type algorithms for frictional multi-contact problems: applications to granular materials",
M. Renouf, P. Alart. doi:10.1016/j.cma.2004.07.009

driver: :cpp:func:`lcp_cpg()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

Projected Gauss-Seidel (:cpp:enumerator:`SICONOS_LCP_PGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :cpp:func:`lcp_pgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

Regularized Projected Gauss-Seidel (:cpp:enumerator:`SICONOS_LCP_RPGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms.

driver: :cpp:func:`lcp_rpgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_RHO] = 1.0

PSOR (:cpp:enumerator:`SICONOS_LCP_PSOR`)
"""""""""""""""""""""""""""""""""""""""""

Projected Successive over relaxation solver for LCP. See :cite:`Cottle.1992`, Chap 5.

driver: :cpp:func:`lcp_psor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_RHO] = 0.1

Latin method (:cpp:enumerator:`SICONOS_LCP_LATIN` and :cpp:enumerator:`SICONOS_LCP_LATIN_W`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Latin stands for LArge Time INcrements.
'w' version is the Latin solver with relaxation.

drivers: :cpp:func:`lcp_latin()` and :cpp:func:`lcp_latin_w()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_LATIN_PARAMETER] = 0.3
* dparam[SICONOS_LCP_DPARAM_RHO] = 1.0 (only useful for solver with relaxation)

Sparse-block Gauss-Seidel (:cpp:enumerator:`SICONOS_LCP_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

Matrix M of the LCP must be a SparseBlockStructuredMatrix.

This solver first build a local problem for each row of blocks and then call any of the other solvers through lcp_driver().

driver: :cpp:func:`lcp_nsgs_SBM()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_LCP_IPARAM_NSGS_ITERATIONS_SUM] (out): sum of all local number of iterations (if it has sense for the local solver)
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_NSGS_LOCAL_ERROR_SUM] (in): sum of all local error values

internal solver : :cpp:enumerator:`SICONOS_LCP_PSOR`

Equation-based solvers
----------------------

Nonsmooth Newton, min formulation (:cpp:enumerator:`SICONOS_LCP_NEWTONMIN`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton method based on the min formulation of the LCP.

.. math::

   0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)

   H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0

References: Alart & Curnier 1990, Pang 1990


driver: :cpp:func:`lcp_newton_min()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6


Nonsmooth Newton, Fisher-Burmeister (:cpp:enumerator:`SICONOS_LCP_NEWTON_FB_FBLSA`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton method based on the Fischer-Bursmeister NCP function.
It uses a variant of line search algorithm (VFBLSA in Facchinei-Pang 2003).

.. math::

   0 \le z \perp w \ge 0 \Longrightarrow \phi(z,w)=\sqrt{z^2+w^2}-(z+w)=0

   \Phi(z) = \left[ \begin{array}{c}  \phi(z_1,w_1) \\ \phi(z_1,w_1) \\ \vdots \\  \phi(z_n,w_n)  \end{array}\right] =0

References: Alart & Curnier 1990, Pang 1990

driver: :cpp:func:`lcp_newton_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0 (in): if > 0. use a non-monotone linear search
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0 (in): if a non-monotone linear search is used, specify the number of merit values to remember
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-10
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;

Nonsmooth Newton, Fisher-Burmeister (:cpp:enumerator:`SICONOS_LCP_NEWTON_MIN_FBLSA`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A nonsmooth Newton method based based on the minFBLSA algorithm : the descent direction is given
by a min reformulation but the linesearch is done with Fischer-Burmeister (and if needed the gradient direction).

driver: :cpp:func:`lcp_newton_minFB()`

parameters: same as :cpp:enumerator:`SICONOS_LCP_NEWTON_FB_FBLSA`.

GAMS solver (:cpp:enumerator:`SICONOS_LCP_GAMS`)
""""""""""""""""""""""""""""""""""""""""""""""""

Optimization solvers from `GAMS <https://www.gams.com/optimization-solvers/>`_.

*Works only if Siconos has been built with GAMS support (see :ref:`siconos_install_guide`)**

driver: :cpp:func:`lcp_gams()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

QP-reformulation
----------------

quadratic program formulation (:cpp:enumerator:`SICONOS_LCP_QP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Quadratic program formulation for solving a LCP

The QP we solve is

Minimize: :math:`z^T (M z + q)` subject to :math:`Mz  + q  \geq  0`

which is the classical reformulation that can be found
in Cottle, Pang and Stone (2009).

If the symmetry condition is not fulfilled, use the NSQP Solver.

driver: :cpp:func:`lcp_qp()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6

quadratic program formulation (:cpp:enumerator:`SICONOS_LCP_NSQP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Non symmetric (and not nonsmooth as one could have thought in a plateform dedicated to nonsmooth problems) quadratic program formulation for solving an LCP with a non symmetric matrix.

driver: :cpp:func:`lcp_nsqp()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6

AVI reformulation
-----------------

AVI with Cao/Ferris solver (:cpp:enumerator:`SICONOS_LCP_AVI_CAOFERRIS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reformulates the LCP as an :ref:`avi_problem`, then uses the solver by Cao and
Ferris.

driver: :cpp:func:`lcp_avi_caoferris()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
