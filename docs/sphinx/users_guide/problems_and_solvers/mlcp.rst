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

The MLCP can also be stored in another form:

Find :math:`(u,v,w)` such that:

:math:`\left\lbrace \begin{array}{l} A u + Cv +a =0\\ D u + Bv +b = w \\ 0 \le v \perp w \ge 0\\ \end{array} \right.`

where A is an ( :math:`n \times n` ) matrix, B is an ( :math:`m \times m` ) matrix, C is an ( :math:`n \times m` ) matrix,

D is an ( :math:`m \times n` ) matrix, a and u are vectors of size :math:`n`, and b, v and w are vectors of size :math:`m`.


Implementation in numerics
==========================

Structure to define the problem: :cpp:class:`MixedLinearComplementarityProblem`.

The generic driver for all MLCPs is :cpp:func:`mlcp_driver()`.

List of solvers: :cpp:enum:`MLCP_SOLVER`


MLCP solvers must be initialized before any call to the driver. The standard sequence of calls is:

#. Initialize the solver with :cpp:func:`mlcp_driver_init()`.

#. Solve the problem with :cpp:func:`mlcp_driver()`.

#. Reset the solver with :cpp:func:`mlcp_driver_reset()`.


.. _mlcp_error:

Error computation
=================

The criterion is based on:

.. math::

   error = \frac{1}{\|q\| }\sum_{i} [ (z[i]*(Mz+q)[i])_{+} + (z[i])_{-} + (Mz+q)[i])_{-} ]

with :math:`x_{+} = max(0,x)` and :math:`x_{-} = max(0,-x)`.

* :cpp:func:`mlcp_compute_error` returns 0 if :math:`error \leq tolerance`, else 1.
* A call to this function updates the content of the input vector w with :math:`Mz + q`.


.. _mlcp_solvers:

MLCP available solvers
======================

Projected Gauss-Seidel (:cpp:enumerator:`SICONOS_MLCP_PGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A basic projected Gauss-Seidel solver for MLCP.

driver: :cpp:func:`mlcp_pgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_MLCP_PGS_EXPLICIT] = 0 (1 for implicit).
* dparam[SICONOS_DPARAM_TOL] = 1e-6


Projected Gauss-Seidel, SBM (:cpp:enumerator:`SICONOS_MLCP_PGS_SBM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Gauss-Seidel with sparse-block storage.

driver: :cpp:func:`mlcp_pgs_sbm()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6

internal solver: :cpp:enumerator:`SICONOS_LCP_PGS`.

Output:

* iparam[SICONOS_IPARAM_MLCP_PGS_SUM_ITER]: sum of the local numbers of iterations (output of local_driver)
* dparam[SICONOS_DPARAM_MLCP_PGS_SUM_ERRORS]: sum of the local errors (output of local_driver)


Regularized Projected Gauss-Seidel (:cpp:enumerator:`SICONOS_MLCP_RPGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Regularized projected Gauss-Seidel solver for MLCP, able to handle matrices with zero diagonal terms.

driver: :cpp:func:`mlcp_rpgs()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_RHO] = 0.5


PSOR (:cpp:enumerator:`SICONOS_MLCP_PSOR`)
""""""""""""""""""""""""""""""""""""""""""

Projected successive over-relaxation solver.

driver: :cpp:func:`mlcp_psor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_OMEGA] = 2


RPSOR (:cpp:enumerator:`SICONOS_MLCP_RPSOR`)
""""""""""""""""""""""""""""""""""""""""""""

Regularized projected successive over-relaxation method.

driver: :cpp:func:`mlcp_rpsor()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[SICONOS_DPARAM_MLCP_OMEGA] = 2
* dparam[SICONOS_DPARAM_MLCP_RHO] = 0.5


PATH (Ferris) solver (:cpp:enumerator:`SICONOS_MLCP_PATH`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

PATH (Ferris) solver.

*Works only if Siconos has been built with PATH support (if PathFerris or PathVI has been found, see* :ref:`siconos_install_guide` *)*

driver: :cpp:func:`mlcp_path()`

parameters:

* dparam[SICONOS_DPARAM_TOL] = 1e-12

Enumerative solver (:cpp:enumerator:`SICONOS_MLCP_ENUM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Brute-force method which tries every possible solution.

driver: :cpp:func:`mlcp_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0 (0: use dgesv, 1: use dgels)
* dparam[SICONOS_DPARAM_TOL] = 1e-12


PATH + enum solver (:cpp:enumerator:`SICONOS_MLCP_PATH_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First tries the PATH (Ferris) solver, then the enumerative solver if PATH fails.

*Works only if Siconos has been built with PATH support (if PathFerris or PathVI has been found, see* :ref:`siconos_install_guide` *)*

driver: :cpp:func:`mlcp_path_enum()`

parameters: same as :cpp:enumerator:`SICONOS_MLCP_ENUM`.


Direct + enum solver (:cpp:enumerator:`SICONOS_MLCP_DIRECT_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

First tries the direct method, then the enumerative solver if it fails.

driver: :cpp:func:`mlcp_direct_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0;
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: a positive value, the tolerance used to decide that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: a positive value, the tolerance used to decide that a variable is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3: number of registered configurations.
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): number of cases where the direct solver failed.



Simplex solver (:cpp:enumerator:`SICONOS_MLCP_SIMPLEX`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :cpp:func:`mlcp_simplex()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Direct/Simplex solver (:cpp:enumerator:`SICONOS_MLCP_DIRECT_SIMPLEX`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Tries the direct method and switches to the simplex method if it fails.

driver: :cpp:func:`mlcp_direct_simplex()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: a positive value, the tolerance used to decide that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: a positive value, the tolerance used to decide that a variable is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): number of cases where the direct solver failed.


Direct/Path solver (:cpp:enumerator:`SICONOS_MLCP_DIRECT_PATH`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Tries the direct method and switches to PATH if it fails.

driver: :cpp:func:`mlcp_direct_path()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: a positive value, the tolerance used to decide that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: a positive value, the tolerance used to decide that a variable is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;

* iparam[7] (out): number of cases where the direct solver failed.


Direct/Path/enum solver (:cpp:enumerator:`SICONOS_MLCP_DIRECT_PATH_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Tries the direct method, then PATH, and finally the enumerative solver.

driver: :cpp:func:`mlcp_direct_path_enum()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS] = 1e-12: a positive value, the tolerance used to decide that complementarity holds.
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG] = 1e-12: a positive value, the tolerance used to decide that a variable is negative.
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS] = 3;
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED] = 0;
* iparam[SICONOS_IPARAM_MLCP_ENUM_USE_DGELS] = 0 (0: use dgesv, 1: use dgels)

* iparam[7] (out): number of cases where the direct solver failed.


Nonsmooth Newton solver, Fischer-Burmeister (:cpp:enumerator:`SICONOS_MLCP_FB`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

driver: :cpp:func:`mlcp_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12

Direct + Nonsmooth Newton solver, Fischer-Burmeister (:cpp:enumerator:`SICONOS_MLCP_DIRECT_FB`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Tries the direct solver, then switches to Fischer-Burmeister.

driver: :cpp:func:`mlcp_direct_FB()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* iparam[SICONOS_IPARAM_MLCP_NUMBER_OF_CONFIGURATIONS]
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_NEG];
* dparam[SICONOS_DPARAM_MLCP_SIGN_TOL_POS];
* iparam[SICONOS_IPARAM_MLCP_UPDATE_REQUIRED]

iparam[7] returns the number of iterations.
