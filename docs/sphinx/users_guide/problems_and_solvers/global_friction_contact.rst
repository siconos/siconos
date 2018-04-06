.. index:: single: Global-Friction-contact problems (2 or 3-dimensional)
.. _doxid-global_fc_problem:

Global-Friction-contact problems (2D or 3D)
===========================================

.. _doxid-global_fc_problem_1pfcIntro:
.. rubric:: Problem statement.:

Given

* a symmetric positive semi definite matrix :math:`{M} \in {{\mathrm{I\!R}}}^{n \times n}`

* a matrix :math:`{H} \in {{\mathrm{I\!R}}}^{n \times {d\, n_c}}`

* a vector :math:`{q} \in {{\mathrm{I\!R}}}^n`

* a vector :math:`{b} \in {{\mathrm{I\!R}}}^{d\, n_c}`

* a vector of coefficients of friction :math:`\mu \in{{\mathrm{I\!R}}}^{n_c}`

the (global or global) frictional contact problem is to find three vectors :math:`v\in{{\mathrm{I\!R}}}^n` , the (global) velocity :math:`u\in{{\mathrm{I\!R}}}^{d\,n_c}` , the relative local velocity and :math:`r\in {{\mathrm{I\!R}}}^{d,n_c}` , the contact forces denoted by :math:`\mathrm{PFC}(M,H,q,b,\mu)` such that

.. math::

    \begin{eqnarray*} \begin{cases} M v = q + H r \\ u = H^\top v + b \\ \hat u = u +\left[ \left[\begin{array}{c} \mu^\alpha \|u^\alpha_{T}\|\\ 0 \\ 0 \end{array}\right]^T, \alpha = 1 \ldots n_c \right]^T \\ \ \ C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu} \end{cases} \end{eqnarray*}

where the Coulomb friction cone is defined by :math:`C_{\mu} = \prod\limits_{\alpha=1\ldots n_c} C^{\alpha}_{\mu^\alpha}`

with :math:`C^{\alpha}_{\mu^\alpha} =\{ r^\alpha, \|r_{t}\| \leq \mu_{\alpha} |r^\alpha_{n}|\}` , and the set :math:`C^{\alpha,\star}_{\mu^\alpha}` its dual.

The modified local velocity :math:`\widehat u ` is not considered as an unknown since it can obtained uniquely from the local velocity :math:`u` . Coulomb's friction law with Signorini's condition for the unilateral contact written in terms of second order complementarity condition

.. math::

    \begin{eqnarray} C^\star_{\mu} \ni {\hat u} \perp r \in C_{\mu} \end{eqnarray}

can be interpreted in a more usual form

.. math::

    \begin{eqnarray} \begin{cases} 0 \leq u_{N} \perp r_N \geq 0 \quad\quad\text{ Signorini condition}\\ u_T = 0 \Rightarrow \|r_T\| \leq \mu |r_n| \quad\quad\text{ Sticking mode} \\ u_T \neq 0 \Rightarrow r_T = - \mu |r_n| \frac{u_T }{\|u_T\|} \quad\quad\text{ Sliding mode} \end{cases} \end{eqnarray}

This problem models any instance of discretized frictional contact problem obtained from

* the time-discretization of dynamics contact problems with event-capturing of event-tracking schemes,

* the time-discretization of quasi-static contact problems,

* the modeling of static contact problems. In this last case, :math:`u` plays the role of the relative displacement at contact

The problem is stored and given to the solver in Siconos/Numerics thanks to a C structure :ref:`GlobalFrictionContactProblem <doxid-struct_global_friction_contact_problem>` .

.. _doxid-global_fc_problem_1pfc3DSolversList:
.. rubric:: Available solvers for Friction Contact 3D:

Use the generic function :ref:`gfc3d_driver() <doxid-_non_smooth_drivers_8h_1a44096d1e6519a0183219f779b1627424>` to call one the the specific solvers listed below:

* :ref:`gfc3d_nsgs() <doxid-gfc3d___solvers_8h_1a8744cb3657964afe57ebe1f68e15305a>` : non-smooth Gauss-Seidel solver (see the functions/solvers list in ``gfc3d_Solvers.h`` )

.. _doxid-global_fc_problem_1pfc3DParam:
.. rubric:: Required and optional parameters:

gfc3d problems needs some specific parameters, given to the :ref:`gfc3d_driver() <doxid-_non_smooth_drivers_8h_1a44096d1e6519a0183219f779b1627424>` function thanks to a SolverOptions structure.

.. index:: single: Global Friction-Contact 3D problems Solvers
.. _doxid-_global_f_c3_d_solvers:

.. rubric:: Global Friction-Contact 3D problems Solvers:

This page gives an overview of the available solvers for friction-contact (3D) problems and their required parameters.

For each solver, the input argument are:

* a :ref:`FrictionContactProblem <doxid-struct_friction_contact_problem>`

* the unknowns (reaction,velocity)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

.. _doxid-_global_f_c3_d_solvers_1pfc3Dnsgs:
.. rubric:: Non-Smooth Gauss Seidel Solver:

function: :ref:`fc3d_nsgs() <doxid-fc3d___solvers_8h_1ad5fdf1a37ff65852645b460ada37cb71>` parameters:

