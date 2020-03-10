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

Structure to define the problem: :class:`LinearComplementarityProblem`.

The generic driver for all LCPs is :func:`linearComplementarity_driver()`
which switches either to :func:`lcp_driver_DenseMatrix` or :func:`lcp_driver_SparseMatrix` according
to the kind of storage used in the problem.

Solvers list  :enum:`LCP_SOLVER`

.. _lcp_error:

Error computation
=================

The criterion is based on :

.. math::

   error = \frac{1}{\|q\| }\sqrt(\sum_{i} (z_i - (z_i - (Mz+q)_i)^{+})^2)

   error = \frac{1}{\|q\| }\sum_{i} [ (z_i*(Mz+q)_i)^{+} + (z_i)^{-} + {(Mz+q)_i}^{-} ]
   
with :math:`x^{+} = max(0,x)` and :math:`x^{-} = max(0,-x)`.



* :func:`lcp_compute_error` returns 0  if  :math:`error \leq tolerance`, else 1.
* A call to this function updates the content of the input vector w with :math:`Mz + q`.

  
.. _lcp_solvers:

LCP available solvers
=====================

Direct solvers
--------------

Lemke (:enumerator:`SICONOS_LCP_LEMKE`)
"""""""""""""""""""""""""""""""""""""""

direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

driver :func:`lcp_lexicolemke`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[2] = 0.0 
* dparam[3] = 0.0

Pivot based methods
"""""""""""""""""""

:enumerator:`SICONOS_LCP_PIVOT`, :enumerator:`SICONOS_LCP_BARD`,
:enumerator:`SICONOS_LCP_MURTY`, :enumerator:`SICONOS_LCP_PATHSEARCH`,
:enumerator:`SICONOS_LCP_PIVOT_LUMOD`
            
generic solver for pivot-based methods: Bard, Murty and Lemke rules are implemented.

drivers:

* :func:`lcp_pivot` for :enumerator:`SICONOS_LCP_PIVOT`, :enumerator:`SICONOS_LCP_BARD`,
and :enumerator:`SICONOS_LCP_MURTY`,
* :func:`lcp_pathsearch` for :enumerator:`SICONOS_LCP_PATHSEARCH`,
* :func:`lcp_pivot_lumod`for :enumerator:`SICONOS_LCP_PIVOT_LUMOD`.


parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] =

  * SICONOS_LCP_PIVOT_BARD for SICONOS_LCP_BARD
  * SICONOS_LCP_PIVOT_LEAST_INDEX for SICONOS_LCP_MURTY
  * SICONOS_LCP_PIVOT_LEMKE for SICONOS_LCP_PIVOT, SICONOS_LCP_PATHSEARCH and SICONOS_LCP_PIVOT_LUMOD

* iparam[SICONOS_IPARAM_PATHSEARCH_STACKSIZE] = 0;
* dparam[SICONOS_DPARAM_TOL] = 100 * epsilon (machine precision)



Enumerative solver (:enumerator:`SICONOS_LCP_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""

Brute-force method which tries every possible solution.

driver: :func:`lcp_enum()`

parameters:

* iparam[SICONOS_LCP_IPARAM_ENUM_USE_DGELS] = 0 (0 : use dgesv, 1: use dgels)
* iparam[SICONOS_LCP_IPARAM_ENUM_SEED] = 0, starting key values
* iparam[SICONOS_LCP_IPARAM_ENUM_MULTIPLE_SOLUTIONS] = 0;
* dparam[SICONOS_DPARAM_TOL] = 1e-6

PATH (:enumerator:`SICONOS_LCP_PATH`)
"""""""""""""""""""""""""""""""""""""

*Works only if Siconos has been built with path support (if PathFerris or PathVI has been found, see :ref:`siconos_install_guide`)*

driver: :func:`lcp_path()`

parameters
* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

 
Iterative solvers
-----------------

Conjugated Projected Gradient (:enumerator:`SICONOS_LCP_CPG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Solver based on quadratic minimization.

driver: :func:`lcp_cpg()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

Projected Gauss-Seidel (:enumerator:`SICONOS_LCP_PGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :func:`lcp_pgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

Regularized Projected Gauss-Seidel (:enumerator:`SICONOS_LCP_RPGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Regularized Projected Gauss-Seidel, is a solver for LCP, able to handle matrices with null diagonal terms.

driver: :func:`lcp_rpgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_RHO] = 1.0

PSOR (:enumerator:`SICONOS_LCP_PSOR`)
"""""""""""""""""""""""""""""""""""""

Projected Succesive over relaxation solver for LCP. See :cite:`Cottle.1992`, Chap 5.

driver: :func:`lcp_psor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_RHO] = 0.1

Latin method (:enumerator:`SICONOS_LCP_LATIN` and :enumerator:`SICONOS_LCP_LATIN_W`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Latin stands for LArge Time INcrements.
'w' version is the Latin solver with relaxation.

drivers: :func:`lcp_latin()` and :func:`lcp_latin_w()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_LCP_DPARAM_LATIN_PARAMETER] = 0.3
* dparam[SICONOS_LCP_DPARAM_RHO] = 1.0 (only useful for solver with relaxation)

Sparse-block Gauss-Seidel (:enumerator:`SICONOS_LCP_NSGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.

driver: :func:`lcp_nsgs_SBM()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

internal solver : :enumerator:`SICONOS_LCP_PSOR`

Equation-based solvers
----------------------

Nonsmooth Newton, min formulation (:enumerator:`SICONOS_LCP_NEWTONMIN`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton method based on the min formulation of the LCP.

.. math::

   0 \le z \perp w \ge 0 \Longrightarrow \min(w,\rho z)=0 \Longrightarrow w = \max(0,w - \rho z)

   H(z) = H(\left[ \begin{array}{c} z \\ w \end{array}\right])= \left[ \begin{array}{c} w-Mz-q \\ min(w,\rho z) \end{array}\right] =0

References: Alart & Curnier 1990, Pang 1990


driver: :func:`lcp_newton_min()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6


Nonsmooth Newton, Fisher-Burmeister (:enumerator:`SICONOS_LCP_NEWTON_FB_FBLSA`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton method based on the Fischer-Bursmeister NCP function.

.. math::

   0 \le z \perp w \ge 0 \Longrightarrow \phi(z,w)=\sqrt{z^2+w^2}-(z+w)=0

   \Phi(z) = \left[ \begin{array}{c}  \phi(z_1,w_1) \\ \phi(z_1,w_1) \\ \vdots \\  \phi(z_n,w_n)  \end{array}\right] =0

References: Alart & Curnier 1990, Pang 1990

driver: :func:`lcp_newton_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-10
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;

Nonsmooth Newton, Fisher-Burmeister (:enumerator:`SICONOS_LCP_NEWTON_MIN_FBLSA`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton method combining the min and FB functions.

driver: :func:`lcp_newton_minFB()`

parameters: same as :enumerator:`SICONOS_LCP_NEWTON_FB_FBLSA`.

GAMS solver (:enumerator:`SICONOS_LCP_GAMS`)
""""""""""""""""""""""""""""""""""""""""""""

Optimization solvers from `GAMS <https://www.gams.com/optimization-solvers/>`_.

*Works only if Siconos has been built with GAMS support (see :ref:`siconos_install_guide`)**

driver: :func:`lcp_gams()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

QP-reformulation
----------------

quadratic programm formulation (:enumerator:`SICONOS_LCP_QP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :func:`lcp_qp()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6

quadratic programm formulation (:enumerator:`SICONOS_LCP_NSQP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Quadratic programm formulation to solve a non symmetric LCP.

driver: :func:`lcp_nsqp()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6

AVI reformulation
-----------------

AVI with Cao/Ferris solver (:enumerator:`SICONOS_AVI_CAOFERRIS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reformulates the LCP as an :ref:`avi_problem`, then uses the solver by Cao and
Ferris.

driver: :func:`lcp_avi_caoferris()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
