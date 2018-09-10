.. index:: single: Mixed Linear Complementarity problems (MLCP)
.. _doxid-_m_l_c_problem:

Mixed Linear Complementarity Problems (MLCP)
============================================

.. _doxid-_m_l_c_problem_1mlcpIntro:
.. rubric:: The problem:

Find :math:`(z,w)` such that:

:math:` \left\{ \begin{array}{l} M \ z + q = w \\ w_1=0 \\ 0 \le w_{2} \perp v \ge 0 \end{array} \right. \text{ with } z= \left[ \begin{array}{c} u\\ v\\ \end{array} \right] \text{ and } w= \left[ \begin{array}{c} w_{1}\\ w_{2}\\ \end{array} \right] `

:math:` u, w_{1}` are vectors of size n.

:math:` v, w_{2}` are vectors of size m.

Another storage is also possible for the MLCP problem:

Try :math:`(u,v,w)` such that:

:math:` \left\lbrace \begin{array}{l} A u + Cv +a =0\\ D u + Bv +b = w \\ 0 \le v \perp w \ge 0\\ \end{array} \right. `

where A is an ( :math:` n \times n` ) matrix, B is an ( :math:` m \times m` ) matrix, C is an ( :math:` n \times m` ) matrix,

D is an ( :math:` m \times n` ) matrix, a and u is an ( :math:` n ` ) vectors b,v and w is an ( :math:` m ` ) vectors.

.. _doxid-_m_l_c_problem_1mlcpSolversList:
.. rubric:: Available solvers:

The solvers and their parameters are described in :ref:`Mixed Linear Complementarity Problems Solvers`.

Use the generic function mlcp_solver(), to call one the the specific solvers listed below:

* :func:`mlcp_pgs()` , (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver

* :func:`mlcp_rpgs()` , (Projected Gauss-Seidel) is a basic Projected Gauss-Seidel solver

* :func:`mlcp_psor()` , projected successive overrelaxation method

* :func:`mlcp_rpsor()` , regularized projected successive overrelaxation method

* :func:`mlcp_path()` , path solver

* :func:`mlcp_enum()` , enumeratif solver

* :func:`mlcp_simplex()` , solver based on the simplex algorithm

* :func:`mlcp_direct_path()` , use the last solution to find the solution. If it failed the mlcp_path is called.

* :func:`mlcp_direct_enum()` , use the last solution to find the solution. If it failed the mlcp_enum is called.

* :func:`mlcp_direct_simplex()` , use the last solution to find the solution. If it failed the mlcp_simplex is called.

Note that all the algorithms are not available for the two options of storage M,q or A,B,C,D,a,b

(see the functions/solvers list in ``MLCP_solvers.h`` )

.. index:: single: Mixed Linear Complementary Problems Solvers
.. _doxid-_m_l_c_p_solvers:

Mixed Linear Complementarity Problems Solvers
=============================================

This page gives an overview of the available solvers for MLCP and their required parameters.

This page gives an overview of the available solvers for MLCP and their required parameters.

For each solver, the input argument are:

* a :class:`MixedLinearComplementarityProblem`

* the unknowns (z,w)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

.. _doxid-_m_l_c_p_solvers_1mlcpPGS:
.. rubric:: PGS Solver:

Projected Gauss-Seidel solver

function: :func:`mlcp_pgs()`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

* iparam[2] (in): 0 for implicit, 1 for explicit

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

.. _doxid-_m_l_c_p_solvers_1mlcpPGS_SBM:
.. rubric:: PGS Solver for SBM storage:

Projected Gauss-Seidel solver

* iparam[0] (in): maximum number of iterations allowed for GS process

* iparam[1] (out): number of GS iterations processed

* iparam[2] (out): sum of all local number of iterations (if it has sense for the local solver)

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

* dparam[2] (in): sum of all local error values

.. _doxid-_m_l_c_p_solvers_1mlcpRPGS:
.. rubric:: RPGS Solver:

Regularized Projected Gauss-Seidel, solver for MLCP, able to handle with matrices with null diagonal terms

function: :func:`mlcp_rpgs()`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

* dparam[2] (in): rho

.. _doxid-_m_l_c_p_solvers_1mlcpPSOR:
.. rubric:: PSOR Solver:

Projected Succesive over relaxation solver for MLCP. See cottle, Pang Stone Chap 5 function: :func:`mlcp_psor()`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

* dparam[2] (in): omega

.. _doxid-_m_l_c_p_solvers_1mlcpRPSOR:
.. rubric:: RPSOR Solver:

Regularized Projected Succesive over relaxation solver for MLCP function: :func:`mlcp_rpsor()`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

* dparam[0] (in): tolerance

* dparam[1] (out): resulting error

* dparam[2] (in): omega

* dparam[3] (in): rho

.. _doxid-_m_l_c_p_solvers_1mlcpPath:
.. rubric:: Path (Ferris) Solver:

The path solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_path()`

parameters:

* dparam[0] (in): tolerance

.. _doxid-_m_l_c_p_solvers_1mlcpENUM:
.. rubric:: ENUM Solver:

The enumeratif solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_enum()`

The enumeratif solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

parameters:

* dparam[0] (in): a positive value, tolerane about the sign.

* iparam[4] (in) : use DGELS (1) or DGESV (0).

* dWork : working float zone size : The number of doubles is retruned by the function :func:`mlcp_driver_get_dwork()`. MUST BE ALLOCATED BY THE USER.

* iWork : working int zone size : . The number of double is retruned by the function :func:`mlcp_driver_get_iwork()`. MUST BE ALLOCATED BY THE USER.

.. _doxid-_m_l_c_p_solvers_1mlcpSIMPLEX:
.. rubric:: SIMPLEX Solver:

The simplex solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_simplex()`

The simplex solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

parameters:

* iparam[0] (in): Max number of iteration (example: 1000000).

* dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).

* dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).

* dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).

* dWork : working float zone size : The number of doubles is retruned by the function :func:`mlcp_driver_get_dwork()`. MUST BE ALLOCATED BY THE USER.

* iWork : working int zone size : . The number of double is retruned by the function :func:`mlcp_driver_get_iwork()`. MUST BE ALLOCATED BY THE USER.

.. _doxid-_m_l_c_p_solvers_1mlcpDIRECT_ENUM:
.. rubric:: DIRECT_ENUM Solver:

The direct and enumeratif solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_direct_enum()`.

parameters:

* iparam[5] (in): Number of registered configurations.

* iparam[7] (out): Number of case the direct solved failed.

* dparam[0] (in): A positive value, tolerane about the sign.

* dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).

* dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).

* dWork : working float zone size : The number of doubles is retruned by the function :func:`mlcp_driver_get_dwork()`. MUST BE ALLOCATED BY THE USER.

* iWork : working int zone size : . The number of double is retruned by the function :func:`mlcp_driver_get_iwork()`. MUST BE ALLOCATED BY THE USER.

.. _doxid-_m_l_c_p_solvers_1mlcpDIRECT_PATH:
.. rubric:: DIRECT_PATH Solver:

The path solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_direct_path()`

parameters:

* iparam[5] (in): Number of registered configurations.

* iparam[7] (out): Number of case the direct solved failed.

* dparam[0] (in): Tolerance.

* dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).

* dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).

* dWork : working float zone size : The number of doubles is retruned by the function :func:`mlcp_driver_get_dwork()`. MUST BE ALLOCATED BY THE USER.

* iWork : working int zone size : . The number of double is retruned by the function :func:`mlcp_driver_get_iwork()`. MUST BE ALLOCATED BY THE USER.

.. _doxid-_m_l_c_p_solvers_1mlcpDIRECT_SIMPLEX:
.. rubric:: DIRECT_SIMPLEX Solver:

The direct and simplex solver must be initialize:

1) Initialize the solver with :func:`mlcp_driver_init()`.

2) Use a lot with :func:`mlcp_driver()`.

3) Reset the solver with :func:`mlcp_driver_reset()`.

function: :func:`mlcp_direct_simplex()`.

parameters:

* iparam[0] (in): Max number of iteration (example: 1000000).

* iparam[5] (in): Number of registered configurations.

* iparam[7] (out): Number of case the direct solved failed.

* dparam[0] (in): A positive value, tolerance to consider that a var is null(ex: 10e-12).

* dparam[1] (in): A positive value, tolerance to consider that complementarity holds(ex: 10e-12).

* dparam[2] (in): A positive value, tolerance to consider that a var is negative(ex: 10e-9).

* dparam[5] (in): A tolerance for the direct solver to consider that a var is negative(ex: 1e-12).

* dparam[6] (in): A tolerance for the direct solver to consider that a var is positive(ex: 1e-12).

* dWork : working float zone size : The number of doubles is retruned by the function :func:`mlcp_driver_get_dwork()`. MUST BE ALLOCATED BY THE USER.

* iWork : working int zone size : . The number of double is retruned by the function :func:`mlcp_driver_get_iwork()`. MUST BE ALLOCATED BY THE USER.

