.. index::
   single: Global friction-contact problems (2 or 3D)

.. contents::

.. _gfc_problem:

Global-Friction-contact problems (2D or 3D)
*******************************************

Problem statement
=================

Given

* a symmetric positive semi definite matrix :math:`{M} \in {{\mathrm{I\!R}}}^{n \times n}`

* a matrix :math:`{H} \in {{\mathrm{I\!R}}}^{n \times {d\, n_c}}`

* a vector :math:`{q} \in {{\mathrm{I\!R}}}^n`

* a vector :math:`{b} \in {{\mathrm{I\!R}}}^{d\, n_c}`

* a vector of coefficients of friction :math:`\mu \in{{\mathrm{I\!R}}}^{n_c}`

the global frictional contact problem denoted by :math:`\mathrm{PFC}(M,H,q,b,\mu)` consists in finding three vectors,

* the global velocity :math:`v\in{{\mathrm{I\!R}}}^n`,

* the relative local velocity :math:`u\in{{\mathrm{I\!R}}}^{d\,n_c}`,

* the contact forces :math:`r\in {{\mathrm{I\!R}}}^{d,n_c}`,

such that :

.. math::

    \begin{eqnarray*} \begin{cases} M v = q + H r \\ u = H^\top v + b \\ \hat u = u +\left[ \left[\begin{array}{c} \mu^\alpha \|u^\alpha_{T}\|\\ 0 \\ 0 \end{array}\right]^T, \alpha = 1 \ldots n_c \right]^T \\ \ \ C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu} \end{cases} \end{eqnarray*}

where the Coulomb friction cone is defined by :math:`C_{\mu} = \prod\limits_{\alpha=1\ldots n_c} C^{\alpha}_{\mu^\alpha}`

with :math:`C^{\alpha}_{\mu^\alpha} =\{ r^\alpha, \|r_{t}\| \leq \mu_{\alpha} |r^\alpha_{n}|\}` , and the set :math:`C^{\alpha,\star}_{\mu^\alpha}` its dual.

The modified local velocity :math:`\widehat u` is not considered as an unknown since it can be obtained uniquely from the local velocity :math:`u`.

Coulomb's friction law with Signorini's condition for the unilateral contact written in terms of second order complementarity condition

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

Structure to define the problem: :class:`GlobalFrictionContactProblem`.

Solvers list : :enum:`FRICTION_SOLVER`, id containing 'GLOBAL_FRICTION_3D'

The generic drivers for global friction-contact problems is :func:`gfc3d_driver`.

.. _gfc_error:

Error strategy
==============

To set internal solver tolerance (when it makes sense!) use :func:`gfc3d_set_internalsolver_tolerance`.

Check details in :ref:`fc_error`


.. _gfc3d_solvers:

Global Friction 3D available solvers
====================================

NSGS (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_NSGS`)
""""""""""""""""""""""""""""""""""""""""""""""""""""

Non-Smooth Gauss Seidel solver with reformulation.

**Driver:** :func:`gfc3d_nsgs`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_NSGS`.

Warning : default iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] is 0, which may lead to
very expensive computation for error checking. Increase this value to improve performances.

Nonsmooth Newton, Alart-Curnier, (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_NSN_AC`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_nonsmooth_Newton_AlartCurnier_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_NSN_AC`.

* iparam[SICONOS_IPARAM_MAX_ITER] = 200;

* iparam[SICONOS_FRICTION_3D_NSN_FORMULATION]

  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_STD (default)
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_STD
  * SICONOS_FRICTION_3D_NSN_FORMULATION_ALARTCURNIER_GENERATED,
  * SICONOS_FRICTION_3D_NSN_FORMULATION_JEANMOREAU_GENERATED
  * SICONOS_FRICTION_3D_NSN_FORMULATION_NULL

* iparam[SICONOS_FRICTION_3D_NSN_MEMORY_ALLOCATED] = 0, 0 if memory for internal work arrays must be allocated, else 1.

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH]
  
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_GOLDSTEINPRICE (default)
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_ARMIJO
  * SICONOS_FRICTION_3D_NSN_LINESEARCH_NO

* iparam[SICONOS_FRICTION_3D_NSN_LINESEARCH_MAX_ITER] = 100  maximum number of iterations allowed for the line search.
* iparam[SICONOS_FRICTION_3D_NSN_MPI_COM] = -1
    
* dparam[SICONOS_DPARAM_TOL] = 1e-10
  
* dparam[SICONOS_FRICTION_3D_NSN_RHO] = 1.

* iparam[SICONOS_FRICTION_3D_IPARAM_ERROR_EVALUATION_FREQUENCY] = 1

PATH (GAMS) (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_GAMS_PATH`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_AVI_gams_path`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_GAMS_PATH`.

PATHVI (GAMS) (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_GAMS_PATHVI`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_AVI_gams_pathvi`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_GAMS_PATHVI`.

Fixed-Point projection (VI reformulation) (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_VI_FPP`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_VI_FixedPointProjection`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_VI_FPP`.

Extra-Gradient (VI reformulation) (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_VI_EG`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_VI_ExtraGradient`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_VI_EG`.

ACLM Fixed point (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_ACLMFP`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_ACLMFixedPoint`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
* iparam[SICONOS_FRICTION_3D_IPARAM_INTERNAL_ERROR_STRATEGY] = SICONOS_FRICTION_3D_INTERNAL_ERROR_STRATEGY_ADAPTIVE
* dparam[SICONOS_DPARAM_TOL] = 1e-4;
* dparam[SICONOS_FRICTION_3D_DPARAM_INTERNAL_ERROR_RATIO] = 2.0

Internal solver: :enumerator:`SICONOS_CONVEXQP_ADMM`, see :ref:`convex_qp_solvers`.

ADMM (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_ADMM`)
""""""""""""""""""""""""""""""""""""""""""""""""""""

Solver based on `ADMM method <https://stanford.edu/~boyd/admm.html>`_.

**Driver:** :func:`gfc3d_ADMM`

**Parameters:**

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000;

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_ACCELERATION]

  * SICONOS_FRICTION_3D_ADMM_ACCELERATION
  * SICONOS_FRICTION_3D_ADMM_ACCELERATION_AND_RESTART (default)
  * SICONOS_FRICTION_3D_ADMM_NO_ACCELERATION

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_SPARSE_STORAGE]

  * SICONOS_FRICTION_3D_ADMM_FORCED_SPARSE_STORAGE
  * SICONOS_FRICTION_3D_ADMM_KEEP_STORAGE (default)

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_INITIAL_RHO] =

  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_GIVEN (default)
  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_NORM_INF
  * SICONOS_FRICTION_3D_ADMM_INITIAL_RHO_EIGENVALUES;
    
* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_RHO_STRATEGY]

  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_RESIDUAL_BALANCING
  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_SCALED_RESIDUAL_BALANCING
  * SICONOS_FRICTION_3D_ADMM_RHO_STRATEGY_CONSTANT (default)

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_GET_PROBLEM_INFO]

  * SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_NO (default)
  * SICONOS_FRICTION_3D_ADMM_GET_PROBLEM_INFO_YES

* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_FULL_H] = SICONOS_FRICTION_3D_ADMM_FULL_H_NO;
* iparam[SICONOS_FRICTION_3D_ADMM_IPARAM_UPDATE_S] = SICONOS_FRICTION_3D_ADMM_UPDATE_S_YES;

* iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING]

  * SICONOS_FRICTION_3D_RESCALING_NO (default)
  * SICONOS_FRICTION_3D_RESCALING_SCALAR,
  * SICONOS_FRICTION_3D_RESCALING_BALANCING_M,
  * SICONOS_FRICTION_3D_RESCALING_BALANCING_MH

* iparam[SICONOS_FRICTION_3D_IPARAM_RESCALING_CONE]=SICONOS_FRICTION_3D_RESCALING_CONE_NO;
* dparam[SICONOS_DPARAM_TOL] = 1e-6;
* dparam[SICONOS_FRICTION_3D_ADMM_RHO] = 0.1;
* dparam[SICONOS_FRICTION_3D_ADMM_RESTART_ETA] = 0.999;
* dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_TAU] = 2.
* dparam[SICONOS_FRICTION_3D_ADMM_BALANCING_RESIDUAL_PHI] = 10.;


Solvers with reformulation
--------------------------

All solvers with id ending with '_WR' (which stands for With Reformulation)
starts with a reformulation of the global problem into a local one, which is solved
with one of the fc3d solvers.


NSGS, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_NSGS_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Non-Smooth Gauss Seidel solver with reformulation.

**Driver:** :func:`gfc3d_nsgs_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_NSGS`.

NSGS, velocity, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_NSGS_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Non-Smooth Gauss Seidel solver with reformulation.

**Driver:** :func:`gfc3d_nsgs_velocity_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_NSGSV`.

Proximal point, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_PROX_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_proximal_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_PROX`.

DeSaxce FixedPoint, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_DSFP_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_DeSaxceFixedPoint_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_DSFP`.

Tresca FixedPoint, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_TFP_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_TrescaFixedPoint_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_TFP`.

Nonsmooth Newton, Alart-Curnier, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_NSN_AC_WR`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_nonsmooth_Newton_AlartCurnier_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_NSN_AC`.

ADMM, with reformulation (:enumerator:`SICONOS_GLOBAL_FRICTION_3D_ADMM_WR`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

**Driver:** :func:`gfc3d_admm_wr`

**Parameters:** same as :enumerate:`SICONOS_FRICTION_3D_ADMM`.

