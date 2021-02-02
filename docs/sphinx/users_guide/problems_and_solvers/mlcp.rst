.. index::
   single: Mixed Linear Complementarity Problems (MLCP)
   
.. contents::

.. _mlcp_problem:

Mixed Linear Complementarity Problems (MLCP)
********************************************

Problem statement
=================

Find :math:`(z,w)` such that:

.. math::

   \left\{ \begin{array}{l}
   M \ z + q = w \\ w_1=0 \\
   0 \le w_{2} \perp v \ge 0
   \end{array} \right. \text{ with } z= \left[ \begin{array}{c} u\\ v\\ \end{array} \right] \text{ and } w= \left[ \begin{array}{c} w_{1}\\ w_{2}\\ \end{array} \right]
   

:math:`u, w_{1}` are vectors of size n.

:math:`v, w_{2}` are vectors of size m.

Another storage is also possible for the MLCP problem:

Try :math:`(u,v,w)` such that:

:math:`\left\lbrace \begin{array}{l} A u + Cv +a =0\\ D u + Bv +b = w \\ 0 \le v \perp w \ge 0\\ \end{array} \right.`

where A is an ( :math:`n \times n` ) matrix, B is an ( :math:`m \times m` ) matrix, C is an ( :math:`n \times m` ) matrix,

D is an ( :math:`m \times n` ) matrix, a and u is an ( :math:`n` ) vectors b,v and w is an ( :math:`m` ) vectors.


Implementation in numerics
==========================

Structure to define the problem: :class:`MixedLinearComplementarityProblem`.

The generic driver for all MLCPs is :func:`mlcp_driver()`.

Solvers list  :enum:`MLCP_SOLVER`


MLCP solvers must be initialized before any call of the driver. Here is the standard sequence of calls:

#. Initialize the solver with :func:`mlcp_driver_init()`.

#. Solve the problemt with :func:`mlcp_driver()`.

#. Reset the solver with :func:`mlcp_driver_reset()`.


.. _mlcp_error:

Error computation
=================

The criterion is based on :

.. math::

   error = \frac{1}{\|q\| }\sum_{i} [ (z[i]*(Mz+q)[i])_{+} + (z[i])_{-} + (Mz+q)[i])_{-} ]
   
with :math:`x_{+} = max(0,x)` and :math:`x_{-} = max(0,-x)`.

* :func:`mlcp_compute_error` returns 0 if :math:`error \leq tolerance`, else 1.
* A call to this function updates the content of the input vector w with :math:`Mz + q`.

 
.. _mlcp_solvers:

MLCP available solvers
======================

Projected Gauss-Seidel (:enumerator:`SICONOS_MLCP_PGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""
Projected Gauss-Seidel,  a basic Projected Gauss-Seidel solver for MLCP.

driver: :func:`mlcp_pgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_MLCP_PGS_EXPLICIT] = 0, 1 for implicit.
* dparam[SICONOS_DPARAM_TOL] = 1e-6


Projected Gauss-Seidel, SBM (:enumerator:`SICONOS_MLCP_PGS_SBM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Gauss-Seidel with sparse-block storage.

driver: :func:`mlcp_pgs_sbm()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

internal solver : :enumerate:`SICONOS_LCP_PGS`.
 
out
* iparam[SICONOS_IPARAM_MLCP_PGS_SUM_ITER], sum of local number of iterations (output from local_driver)
* dparam[SICONOS_DPARAM_MLCP_PGS_SUM_ERRORS] sum of local errors (output from local_driver)


Regularized Projected Gauss-Seidel (:enumerator:`SICONOS_MLCP_RPGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Regularized Projected Gauss-Seidel, solver for MLCP, able to handle with matrices with null diagonal terms

driver: :func:`mlcp_rpgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_RHO] = 0.5


PSOR (:enumerator:`SICONOS_MLCP_PSOR`)
""""""""""""""""""""""""""""""""""""""

Projected Succesive over relaxation solver.

driver: :func:`mlcp_psor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_OMEGA] = 2


RPSOR (:enumerator:`SICONOS_MLCP_RPSOR`)
""""""""""""""""""""""""""""""""""""""""

Regularized projected successive overrelaxation method.

driver: :func:`mlcp_rpsor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_OMEGA] = 2
* dparam[SICONOS_DPARAM_MLCP_RHO] = 0.5


PATH (Ferris) solver (:enumerator:`SICONOS_MLCP_PATH`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Path (Ferris) Solver.

*Works only if Siconos has been built with path support (if PathFerris or PathVI has been found, see :ref:`siconos_install_guide`)*

driver: :func:`mlcp_path()`

parameters
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Enumerative solver (:enumerator:`SICONOS_MLCP_ENUM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""

Brute-force method which tries every possible solution.

driver: :func:`mlcp_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0 (0 : use dgesv, 1: use dgels)
* dparam[SICONOS_DPARAM_TOL] = 1e-12


PATH + enum solver (:enumerator:`SICONOS_MLCP_PATH_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First try with Path (Ferris) Solver then use enum if the solver failed.

*Works only if Siconos has been built with path support (if PathFerris or PathVI has been found, see :ref:`siconos_install_guide`)*

driver: :func:`mlcp_path_enum()`

parameters : same as :enumerator:`SICONOS_MLCP_ENUM`.

  
Direct + enum solver (:enumerator:`SICONOS_MLCP_DIRECT_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First try direct method and then use enum if the solver failed.
driver: :func:`mlcp_direct_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: A positive value, tolerance to consider that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: A positive value, tolerance to consider that a var is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3 : Number of registered configurations.
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): Number of case the direct solved failed.


  
Simplex solver (:enumerator:`SICONOS_MLCP_SIMPLEX`)
"""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :func:`mlcp_simplex()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Direct/Simplex solver (:enumerator:`SICONOS_MLCP_DIRECT_SIMPLEX`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Try direct method and switch to simplex if it fails.

driver: :func:`mlcp_direct_simplex()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: A positive value, tolerance to consider that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: A positive value, tolerance to consider that a var is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): Number of case the direct solved failed.


Direct/Path solver (:enumerator:`SICONOS_MLCP_DIRECT_PATH`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Try direct method and switch to Path if it fails.

driver: :func:`mlcp_direct_path()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: A positive value, tolerance to consider that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: A positive value, tolerance to consider that a var is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): Number of case the direct solved failed.


Direct/Path/enum solver (:enumerator:`SICONOS_MLCP_DIRECT_PATH_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Try direct then switch to PATH and finish with enum.

driver: :func:`mlcp_direct_path_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: A positive value, tolerance to consider that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: A positive value, tolerance to consider that a var is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;
* iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0 (0 : use dgesv, 1: use dgels)

* iparam[7] (out): Number of case the direct solved failed.


Nonsmooth Newton solver, Fisher-Burmeister (:enumerator:`SICONOS_MLCP_FB`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :func:`mlcp_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Direct + Nonsmooth Newton solver, Fisher-Burmeister (:enumerator:`SICONOS_MLCP_DIRECT_FB`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Try direct solver then switch to Fisher-Burmeister.

driver: :func:`mlcp_direct_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS]
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG];
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS];
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED]

return iparam[7] for iters

