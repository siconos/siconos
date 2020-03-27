.. index::
   single: Friction-contact problems (2 or 3D)

.. contents::

.. _fc_problem:

Friction-contact problems (2 or 3D)
***********************************

Problem statement
=================


Given

* a symmetric positive semi definite matrix :math:`{M} \in {{\mathrm{I\!R}}}^{n \times n}`

* a vector :math:`{q} \in {{\mathrm{I\!R}}}^n`

* a vector of coefficients of friction :math:`\mu \in{{\mathrm{I\!R}}}^{n_c}`

the (reduced or dual) frictional contact problem is to find two vectors :math:`u\in{{\mathrm{I\!R}}}^n` , the relative local velocity and :math:`r\in {{\mathrm{I\!R}}}^n` , the contact forces denoted by :math:`\mathrm{FC}(M,q,\mu)` such that

.. math::

    \begin{eqnarray*} \begin{cases}
    u = M r + q \\
    \hat u = u +\left[ \left[\begin{array}{c} \mu^\alpha \|u^\alpha_{T}\|\\ 0 \\ 0 \end{array}\right]^T, \alpha = 1 \ldots n_c \right]^T \\ \ \ C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu}
    \end{cases} \end{eqnarray*}

where the Coulomb friction cone is defined by :math:`C_{\mu} = \prod\limits_{\alpha=1\ldots n_c} C^{\alpha}_{\mu^\alpha}`

with :math:`C^{\alpha}_{\mu^\alpha} =\{ r^\alpha, \|r_{t}\| \leq \mu_{\alpha} |r^\alpha_{n}|\}` , and the set :math:`C^{\alpha,\star}_{\mu^\alpha}` its dual.

The modified local velocity :math:`\widehat u ` is not considered as an unknown since it can be obtained uniquely from the local velocity :math:`u` . Coulomb's friction law with Signorini's condition for the unilateral contact written in terms of second order complementarity condition

.. math::

    \begin{eqnarray} C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu} \end{eqnarray}

can be interpreted in a more usual form

.. math::

    \begin{eqnarray} \begin{cases} 0 \leq u_{N} \perp r_N \geq 0 \quad\quad\text{ Signorini condition}\\ u_T = 0 \Rightarrow \|r_T\| \leq \mu |r_n| \quad\quad\text{ Sticking mode} \\ u_T \neq 0 \Rightarrow r_T = - \mu |r_n| \frac{u_T }{\|u_T\|} \quad\quad\text{ Sliding mode} \end{cases} \end{eqnarray}

This problem models any instance of discretized frictional contact problem obtained from

* the time-discretization of dynamics contact problems with event-capturing of event-tracking schemes,

* the time-discretization of quasi-static contact problems,

* the modeling of static contact problems. In this last case, :math:`u` plays the role of the relative displacement at contact

Implementation in numerics
==========================

Structure to define the problem: :class:`FrictionContactProblem`.

Solvers list  :enum:`FRICTION_SOLVER`

The generic drivers for friction-contact problems are:

* :func:`fc2d_driver` (id contains FRICTION_2D)
* :func:`fc3d_driver` (id contains FRICTION_3D)
* :func:`gfc3d_driver` (id contains GLOBAL_FRICTION)
* :func:`rolling_fc3d_driver` (id contains ROLLING_FRICTION_3D)


For details regarding global formulation and rolling-friction problems, see :ref:`gfc_problem` or :ref:`rfc_problem`.
  
.. _fc_error:

Error strategy
==============

To set internal solver tolerance (when it makes sense!) use one of the following functions :

:func:`fc3d_set_internalsolver_tolerance`, :func:`gfc3d_set_internalsolver_tolerance`, :func:`rolling_fc3d_set_internalsolver_tolerance`

The computation of the tolerance depends on the value of iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY].

It can be:

* SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE

  internal solver tolerance = error/dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO]
  
* SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE_N_CONTACT

  internal solver tolerance = error/dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] * number of contacts
  
* SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE
    
  internal solver tolerance = value provided during initialisation of the local solver.

Warning : iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] and dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] must be set properly for all solvers that are using Xfc3d_set_internal_tolerance function.
  
  
.. _fc2d_solvers:

Friction 2D available solvers
=============================

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_2D_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

**Driver:** :func:`fc2d_nsgs`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_NSGS_SHUFFLE] = 0
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] = SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL;
* dparam[SICONOS_DPARAM_TOL] = 1e-4

Conjugated projected gradient (:enumerator:`SICONOS_FRICTION_2D_CPG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc2d_cpg`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] = 1e-4

Lemke solver (:enumerator:`SICONOS_FRICTION_2D_LEMKE`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Convert friction problem into a LCP and solve it using Lemke solver.

**Driver:** :func:`fc2d_lemke`

**Parameters:** same as :enumerator:`SICONOS_LCP_LEMKE`, see :ref:`lcp_solvers`.


Enumerative solver (:enumerator:`SICONOS_FRICTION_2D_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Convert friction problem into a LCP and solve it using enumerative solver.

**Driver:** :func:`fc2d_enum`

**Parameters:** same as :enumerator:`SICONOS_LCP_ENUM`, see :ref:`lcp_solvers`.

.. _fc3d_solvers:

Friction 3D available solvers
=============================

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_3D_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
**Driver:** :func:`fc3d_nsgs`

**Parameters:**


* iparam[SICONOS_IPARAM_MAX_ITER] = 1000 : Maximum iteration number
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION] : error computation method,
  
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_FULL : Full error computation with velocity computation
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL (DEFAULT): Light error computation with incremental values on reaction verification of absolute error at the end
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_LIGHT : only light error computation (velocity not computed)
  * SICONOS_FRICTION_3D_NSGS_ERROR_EVALUATION_ADAPTIVE :  we adapt the frequency of the full erro evaluation.

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 0,  error computation frequency

* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE] : shuffle the contact indices in the loop
  
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_FALSE : no shuffle
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE : shuffle only at the beginning
  * SICONOS_FRICTION_3D_NSGS_SHUFFLE_TRUE_EACH_LOOP : shuffle in each iteration

* iparam[SICONOS_FRICTION_3D_NSGS_SHUFFLE_SEED] = 0 : seed for the random generator in shuffling  contacts

* iparam[SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION] : filter local solution if the local error is greater than 1.0

  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_FALSE (default) the filter is not applied
  * SICONOS_FRICTION_3D_NSGS_FILTER_LOCAL_SOLUTION_TRUE  the filter is applied

* iparam[SICONOS_FRICTION_3D_NSGS_RELAXATION] : method uses overrelaxation

  * SICONOS_FRICTION_3D_NSGS_RELAXATION_FALSE (default) relaxation is not used,
  * SICONOS_FRICTION_3D_NSGS_RELAXATION_TRUE  relaxation is used with parameter dparam[8],

  
* dparam[SICONOS_DPARAM_TOL] = 1e-4, user tolerance on the loop
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0
* dparam[SICONOS_FRICTION_3D_NSGS_RELAXATION_VALUE]  the relaxation parameter omega

out

*  iparam[SICONOS_IPARAM_ITER_DONE] = iter number of performed iterations
* dparam[SICONOS_DPARAM_RESIDU]  reached error

Default internal solver : :enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID`.
      


Nonsmooth Gauss-Seidel, velocity version (:enumerator:`SICONOS_FRICTION_3D_NSGSV`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_nsgs_velocity`

**Parameters:**


* iparam[SICONOS_IPARAM_MAX_ITER] = 1000 : Maximum iteration number
* dparam[SICONOS_DPARAM_TOL] = 1e-4, user tolerance on the loop
out

*  iparam[7] as number of performed iterations

 dparam[SICONOS_DPARAM_RESIDU(1)]  reached error

Default internal solver : :enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN`.
      

Proximal point solver (:enumerator:`SICONOS_FRICTION_3D_PROX`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_proximal`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000 : Maximum iteration number
* iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_STRATEGY]

  * SICONOS_FRICTION_3D_PROXIMAL_REGULARIZATION) 
  * SICONOS_FRICTION_3D_PROXIMAL_PROX (default)
* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0
    
* dparam[SICONOS_DPARAM_TOL] = 1e-4, user tolerance on the loop
* dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_ALPHA] = 1e4
* dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA] = 5.
* dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_NU] = 1.

out

iparam[SICONOS_FRICTION_3D_PROXIMAL_IPARAM_CUMULATIVE_ITER_DONE]

Default internal solver : :enumerator:`SICONOS_FRICTION_3D_NSN_AC`.

Fixed-point (Tresca) (:enumerator:`SICONOS_FRICTION_3D_TFP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Fixed point solver for friction-contact 3D problem based on the Tresca
  problem with fixed friction threshold

**Driver:** :func:`fc3d_TrescaFixedPoint`

**Parameters:**


* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;
* dparam[SICONOS_DPARAM_TOL] = 1e-14
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;


Default internal solver : :enumerator:`SICONOS_FRICTION_3D_NSGS` with
:enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinderWithLocalIteration`
as internal solver.

Nonsmooth Newton/ Alart-Curnier (:enumerator:`SICONOS_FRICTION_3D_NSN_AC`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_nonsmooth_Newton_AlartCurnier`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 200;

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD 
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED, (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL

* iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY]

  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO (default)
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP : Loop PLI-NSN strategy 
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP : NSN and after Loop PLI-NSN strategy for the hybrid solver 
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN : VI_EG preconditionning to NSN

* iparam[3] = 100000; /* nzmax*/
* iparam[5] = 1;

* iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT : uses constant value (dparam[SICONOS_FRICTION_3D_NSN_RHO]) for rho
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE

* iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1,  mpi com fortran 

* iparam[SICONOS_FRICTION_3D_NSN_LINEAR_SOLVER] Linear solver used at each Newton iteration
  * SICONOS_FRICTION_3D_NSN_USE_CSLUSOL
  * SICONOS_FRICTION_3D_NSN_USE_MUMPS

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1; (must be > 0 !)

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]
  
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE (default)
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_NO

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100  maximum number of iterations allowed for the line search.
 
* dparam[SICONOS_DPARAM_TOL] = 1e-3
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1


Nonsmooth Newton/ Alart-Curnier (test) (:enumerator:`SICONOS_FRICTION_3D_NSN_AC_TEST`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_nonsmooth_Newton_AlartCurnier2`

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED,
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL

* iparam[SICONOS_IPARAM_LSA_SEARCH_CRITERION] = SICONOS_LSA_GOLDSTEIN;
* iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY]
* optionsiparam[SICONOS_IPARAM_MAX_ITER] = 1000;
* optionsdparam[SICONOS_DPARAM_TOL] = 1e-10;

* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0;
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0;
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16;
* dparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_RESIDU;

   
Fixed-Point (De Saxce formulation) (:enumerator:`SICONOS_FRICTION_3D_DSFP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_DeSaxceFixedPoint`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;

  
Fixed-Point projection (VI reformulation) (:enumerator:`SICONOS_FRICTION_3D_VI_FPP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_VI_FixedPointProjection`

**Parameters:** same as :enumerator:`SICONOS_VI_FPP`, see :ref:`vi_solvers`.

Fixed-Point projection on cylinder (VI reformulation) (:enumerator:`SICONOS_FRICTION_3D_VI_FPP_Cylinder`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_VI_FixedPointProjection_Cylinder`

**Parameters:** same as :enumerator:`SICONOS_VI_FPP`, see :ref:`vi_solvers`.

Extra Gradient (VI reformulation) (:enumerator:`SICONOS_FRICTION_3D_VI_EG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_VI_ExtraGradient`

**Parameters:** same as :enumerator:`SICONOS_VI_EG`, see :ref:`vi_solvers`.


Hyperplane Projection (:enumerator:`SICONOS_FRICTION_3D_HP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_HyperplaneProjection`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 50.;

* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_PROXIMAL_DPARAM_SIGMA] = 0.99
  
Fixed-Point projection (:enumerator:`SICONOS_FRICTION_3D_FPP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_fixedPointProjection`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;


Extra Gradient (:enumerator:`SICONOS_FRICTION_3D_EG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_ExtraGradient`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;
* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = -1.;


Nonsmooth Newton (Fischer-Burmeister formulation) (:enumerator:`SICONOS_FRICTION_3D_NSN_FB`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_nonsmooth_Newton_FischerBurmeister`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 200;

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED,
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL
  
* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]
  
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE (default)
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_NO
    
* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100;

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1; (must be > 0 !)

* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;


PATH (via GAMS) + AVI reformulation (:enumerator:`SICONOS_FRICTION_3D_GAMS_PATH`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
solver using PATH (via GAMS) for friction-contact 3D problem based on an AVI reformulation

**Driver:** :func:`fc3d_AVI_gams_path`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
* dparam[SICONOS_DPARAM_TOL] = 1e-9;

out

* dparam[TOTAL_TIME_USED]
* iparam[TOTAL_ITER]
* iparam[LAST_MODEL_STATUS]
* iparam[LAST_SOLVE_STATUS]

PATHVI (via GAMS) + AVI reformulation :enumerator:`SICONOS_FRICTION_3D_GAMS_PATHVI`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
solver using PATHVI (via GAMS) for friction-contact 3D problem based on an AVI reformulation

**Driver:** :func:`fc3d_AVI_gams_pathvi`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
* dparam[SICONOS_DPARAM_TOL] = 1e-9;

out

* dparam[TOTAL_TIME_USED]
* iparam[TOTAL_ITER]
* iparam[LAST_MODEL_STATUS]
* iparam[LAST_SOLVE_STATUS]

 ACLM Fixed-Point (:enumerator:`SICONOS_FRICTION_3D_ACLMFP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_ACLMFixedPoint`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_GIVEN_VALUE
* dparam[SICONOS_DPARAM_TOL] = 1e-4;
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 10.0

Internal solver: :enumerator:`SICONOS_SOCLCP_NSGS`, see :ref:`soclcp_solvers`.

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_3D_SOCLCP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
**Driver:** :func:`fc3d_SOCLCP`

**Parameters:** same as :enumerator:`SICONOS_SOCLCP_NSGS`, see : ref:`soclcp_solvers`.


PATH (via GAMS) + LCP reformulation (:enumerator:`SICONOS_FRICTION_3D_GAMS_LCP_PATH`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
solver using PATH (via GAMS) for friction-contact 3D problem based on an LCP reformulation

**Driver:** :func:`fc3d_lcp_gams_path`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
* dparam[SICONOS_DPARAM_TOL] = 1e-9;

PATHVI (via GAMS) + LCP reformulation :enumerator:`SICONOS_FRICTION_3D_GAMS_LCP_PATHVI`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
solver using PATHVI (via GAMS) for friction-contact 3D problem based on an LCP reformulation

**Driver:** :func:`fc3d_lcp_gams_pathvi`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000;
* dparam[SICONOS_DPARAM_TOL] = 1e-9;

Nonsmooth Newton, Natural Map (:enumerator:`SICONOS_FRICTION_3D_NSN_NM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Nonsmooth Newton solver based on the Natural--Map function for
the local (reduced) frictional contact problem in the dense form.

**Driver:** :func:`fc3d_nonsmooth_Newton_NaturalMap`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 200;

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED,
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]
  
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE (default)
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_NO
    
* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100;
* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1;

* dparam[SICONOS_DPARAM_TOL] = 1e-3;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.;


Fixed point, Panagiotopoulos (:enumerator:`SICONOS_FRICTION_3D_PFP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Fixed point solver for friction-contact 3D problem based on the Panagiotopoulos
method based on an alternative technique between the normal problem and the tangential one.

**Driver:** :func:`fc3d_Panagiotopoulos_FixedPoint`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 200;
* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] =  SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE;

* dparam[SICONOS_DPARAM_TOL] = 1e-4;
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] =10.0;

Two internal solvers: :enumerator:`SICONOS_LCP_PGS` and :enumerator:`SICONOS_CONVEXQP_VI_FPP`.

ADMM (:enumerator:`SICONOS_FRICTION_3D_ADMM`)
"""""""""""""""""""""""""""""""""""""""""""""

Solver based on `ADMM method <https://stanford.edu/~boyd/admm.html>`_.

**Driver:** :func:`fc3d_admm`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SYMMETRY]
  
  * SICONOS_FRICTION_3D_ADMM_FORCED_SYMMETRY (default)
  * SICONOS_FRICTION_3D_ADMM_FORCED_ASYMMETRY
  * SICONOS_FRICTION_3D_ADMM_CHECK_SYMMETRY

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]

  * SICONOS_FRICTION_3D_ADMM_ACCELERATION
  * SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART (default)
  * SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE]

  * SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE
  * SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE (default)

* dparam[SICONOS_FRICTION_3D_ADMM_INITIAL_RHO] =

  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN (default)
  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF
  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]

  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING
  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT (default)

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO]

  * SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO (default)
  * SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES
    
* iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]

  * SICONOS_FRICTION_3D_RESCALING_NO (default)
  * SICONOS_FRICTION_3D_RESCALING_SCALAR,
  * SICONOS_FRICTION_3D_RESCALING_BALANCING_M,
  * SICONOS_FRICTION_3D_RESCALING_BALANCING_MH

* dparam[SICONOS_DPARAM_TOL] = 1e-6;
* dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 1.
* dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA] = 0.999;
* dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU] = 2.
* dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI] = 2.;


"One contact" solvers
^^^^^^^^^^^^^^^^^^^^^

Newton(:enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN`, ...)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`fc3d_onecontact_nonsmooth_Newton_solvers_solve`

which switches to one of the local drivers below:

.. csv-table:: Projection on cone solvers
   :header: "Solver id", "Driver"
   :widths: 15, 30

   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN`",":func:`fc3d_onecontact_nonsmooth_Newton_solvers_solve_direc`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN_GP`",":func:`fc3d_onecontact_nonsmooth_Newton_solvers_solve_dampe`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_NSN_GP_HYBRID`",":func:`fc3d_onecontact_nonsmooth_Newton_solvers_solve_hybrid`"

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 10
* iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED,
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL

* iparam[SICONOS_FRICTION_3D_NSN_RHO_STRATEGY]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_CONSTANT : uses constant value (dparam[SICONOS_FRICTION_3D_NSN_RHO]) for rho
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPECTRAL_NORM
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM_COND (default for NSN)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_SPLIT_SPECTRAL_NORM (default for NSN_GP and NSN_GP_HYBRID)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_RHO_STRATEGY_ADAPTIVE

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]
  
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE (default for NSN_GP and NSN_GP_HYBRID)
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_NO (default for NSN)

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 10;

* iparam[SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY]

  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NO
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_PLI_NSN_LOOP (default for NSN)
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_NSN_AND_PLI_NSN_LOOP (default for NSN_GP and NSP_GP_HYBRID)
  * SICONOS_FRICTION_3D_NSN_HYBRID_STRATEGY_VI_EG_NSN

* iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_LOOP] = 1;
* iparam[SICONOS_FRICTION_3D_NSN_HYBRID_MAX_ITER] = 10 (for NSN), 100 (for NSN_GP and NSN_GP_HYBRID);

    
* dparam[SICONOS_DPARAM_TOL] =1e-14;
* dparam[SICONOS_FRICTION_3D_NSN_RHO] =1.0;
 

Projection on cone or cylinder (:enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone`, ...)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. csv-table:: Projection on cone solvers
   :header: "Solver id", "Driver"
   :widths: 15, 30

   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone`",":func:`fc3d_projectionOnCone_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithRegularization`",":func:`fc3d_projectionOnCone_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithLocalIteration`",":func:`fc3d_projectionOnConeWithLocalIteration_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnConeWithDiagonalization`",":func:`fc3d_projectionOnConeWithDiagonalization_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCone_velocity`",":func:`fc3d_projectionOnCone_velocity_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnCylinder`",":func:`fc3d_projectionOnCylinder_solve`"
   ":enumerator:`SICONOS_FRICTION_3D_ONECONTACT_ProjectionOnProjectionOnCylinderWithLocalIteration`",":func:`fc3d_projectionOnCylinderWithLocalIteration_solve`"
  
**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_FRICTION_3D_CURRENT_CONTACT_NUMBER]
* dparam[SICONOS_DPARAM_TOL] =1e-14
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1., used only in the 'with regularization' case


NCP Fixed Point solver (:enumerator:`SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint`, ...)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


.. csv-table:: NCP Fixed-point solvers
   :header: "Solver id", "Driver", "Update"
   :widths: 15, 30, 30

   ":enumerator:`SICONOS_FRICTION_3D_NCPGlockerFBFixedPoint`",":func:`fc3d_FixedP_solve`",                                ":func:`NCPGlocker_update`"
   ":enumerator:`SICONOS_FRICTION_3D_NCPGlockerFBNewton`",    ":func:`fc3d_onecontact_nonsmooth_Newton_solvers_solve`",   ":func:`NCPGlocker_update`"
   ":enumerator:`SICONOS_FRICTION_3D_NCPGlockerFBPATH`",      ":func:`fc3d_Path_solve`",                                  ":func:`NCPGlocker_update`"

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* dparam[SICONOS_DPARAM_TOL] =1e-12


Quartic (:enumerator:`SICONOS_FRICTION_3D_ONECONTACT_QUARTIC`, ...)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

ids : :enumerator:`SICONOS_FRICTION_3D_ONECONTACT_QUARTIC`, :enumerator:`SICONOS_FRICTION_3D_ONECONTACT_QUARTIC_NU`

**Driver:** :func:`fc3d_unitary_enumerative`

**Parameters:**

* dparam[SICONOS_DPARAM_TOL] =1e-12





As Convex QP (:enumerator:`SICONOS_FRICTION_3D_CONVEXQP_CYLINDER`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Reformulate the problem as a convex QP and solve using :enumerator:`SICONOS_CONVEXQP_PG`.


**Driver:** :func:`fc3d_ConvexQP_ProjectedGradient_Cylinder`

**Parameters:** same as :enumerator:`SICONOS_CONVEXQP_PG, see :ref:`convex_qp_solvers`.

