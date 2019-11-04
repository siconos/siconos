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
* :func:`gfc3d_driver` (id contains GLOBAL_FRICTION_2D)
* :func:`rolling_fc3d_driver` (id contains ROLLING_FRICTION_3D)

  
.. _fc2d_solvers:

Friction 2D available solvers
=============================

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_2D_NSGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

direct solver for LCP based on pivoting method principle for degenerate problem: the choice of pivot variable is performed via lexicographic ordering.

driver :func:`lcp_lexicolemke`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
* iparam[SICONOS_LCP_IPARAM_PIVOTING_METHOD_TYPE] = 0
* dparam[SICONOS_DPARAM_TOL] = 1e-6
* dparam[2] = 0.0 
* dparam[3] = 0.0

Nonsmooth Gauss-Seidel (:enumerator:`SICONOS_FRICTION_2D_PGS`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Conjugated projected gradient (:enumerator:`SICONOS_FRICTION_2D_CPG`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Latin method (:enumerator:`SICONOS_FRICTION_2D_LATIN`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""
Lemke solver (:enumerator:`SICONOS_FRICTION_2D_LEMKE`)
""""""""""""""""""""""""""""""""""""""""""""""""""""""
Enumerative solver (:enumerator:`SICONOS_FRICTION_2D_ENUM`)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

  
.. _fc3d_solvers:

Friction 3D available solvers
=============================

.. _gfc3d_solvers:

Global Friction 3D available solvers
====================================

.. _rfc3d_solvers:

Rolling-Friction 3D available solvers
=====================================





* :ref:`fc3d_nsgs() <doxid-fc3d___solvers_8h_1ad5fdf1a37ff65852645b460ada37cb71>` : non-smooth Gauss-Seidel solver. SolverId : SICONOS_FRICTION_3D_NSGS =500,

* :ref:`fc3d_nsgs_velocity() <doxid-fc3d___solvers_8h_1a84434da5b88d4b8760689d1bfa962342>` : non-smooth Gauss-Seidel solver based on velocity updates SolverId : SICONOS_FRICTION_3D_NSGSV =501,

* :ref:`fc3d_proximal() <doxid-fc3d___solvers_8h_1afdb77f4dd3a5e466b250ad25a0cff36d>` : Proximal point solver for friction-contact 3D problem SolverId : SICONOS_FRICTION_3D_PROX =502,

* :ref:`fc3d_TrescaFixedPoint() <doxid-fc3d___solvers_8h_1a39c8157d2a05ad422dfb327fcc2e209e>` : Fixed point solver for friction-contact 3D problem based on the Tresca problem with fixed friction threshold SolverId : SICONOS_FRICTION_3D_TFP =503,

* fc3d_globalAlartCurnier() : Global Alart Curnier solver SolverId : SICONOS_FRICTION_3D_NSN_AC =504,

* :ref:`fc3d_DeSaxceFixedPoint() <doxid-fc3d___solvers_8h_1aadc2e7d0e90773b4eb10d93085d72bea>` : Fixed Point solver for friction-contact 3D problem based on the De Saxce Formulation SolverId : SICONOS_FRICTION_3D_DSFP=505,

* :ref:`fc3d_ExtraGradient() <doxid-fc3d___solvers_8h_1a9d21d7d3c9beea711644e6b27648187b>` : Extra Gradient solver for friction-contact 3D problem based on the De Saxce Formulation SolverId : SICONOS_FRICTION_3D_EG=506,

* :ref:`fc3d_HyperplaneProjection() <doxid-fc3d___solvers_8h_1ab6575d95a1bf15da12ca0ad9abe4d4bf>` : Hyperplane Projection solver for friction-contact 3D problem based on the De Saxce Formulation SolverId : SICONOS_FRICTION_3D_HP=507,

(see the functions/solvers list in ``fc3d_Solvers.h`` )

.. _doxid-fc_problem_1fc3DParam:
.. rubric:: Required and optional parameters:

fc3d problems needs some specific parameters, given to the :ref:`fc3d_driver() <doxid-_non_smooth_drivers_8h_1ae89ab73414684ab73c0974d1a146bdc8>` function thanks to a SolverOptions structure.

.. _doxid-fc_problem_1fc2DSolversList:
.. rubric:: Available solvers for Friction Contact 2D:

* :ref:`fc2d_nsgs() <doxid-fc2d___solvers_8h_1a5f338a862ee4b2105b923d4eea9a8768>` , Non Linear Gauss Seidel solver. SolverId SICONOS_FRICTION_2D_NSGS =400,

* :ref:`fc2d_cpg() <doxid-fc2d___solvers_8h_1a5e61270d2465dd97040a52f19e679871>` , conjugate projected gradient SolverId SICONOS_FRICTION_2D_CPG =401,

* fc2d_pgs(), projected Gauss Seidel solver. SolverId SICONOS_FRICTION_2D_PGS =402,

* :ref:`fc2d_latin() <doxid-fc2d___solvers_8h_1a1ef4633f903533150224f306435060ef>` , latin solver. SolverId SICONOS_FRICTION_2D_LATIN =403,

* :ref:`fc2d_lexicolemke() <doxid-fc2d___solvers_8h_1a78b1831d4dbd5c364240885cf0fff6e6>` , lemke solver. SolverId SICONOS_FRICTION_2D_LMEKE =404,

