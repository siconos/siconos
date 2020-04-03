.. index::
   single: Variational Inequality (VI)
   
.. contents::

.. _vi_problem:

Variational Inequality (VI)
***************************

Problem statement
=================

Given

* an integer :math:`n` , the dimension of the ambient space,

* a mapping :math:`F\colon \mathrm{I\!R}^n \rightarrow \mathrm{I\!R}^n`

* a set :math:`{X} \in {{\mathrm{I\!R}}}^n`

the variational inequality problem consists in finding a vector :math:`z\in{{\mathrm{I\!R}}}^n` , such that

.. math::

    \begin{equation*} F(z)^T(y-z) \geq 0,\quad \text{ for all } y \in X \end{equation*}

or equivalently,

.. math::

    \begin{equation*} - F(z) \in \mathcal{N}_X(z) \end{equation*}

where :math:`\mathcal{N}_X` is the normal cone to :math:`X` at :math:`z` .

*References* : :cite:`Facchinei.2003`, :cite:`Acary.Brogliato2008`.

Implementation in numerics
==========================

Structure to define the problem: :struct:`VariationalInequality`.

The generic driver for all VI is :func:`variationalInequality_driver()`.

solvers list  :enum:`VI_SOLVER`

.. _vi_solvers:

VI Available solvers
====================

Extra gradient (:enumerator:`SICONOS_VI_EG`)
--------------------------------------------

Extra Gradient solver forvariational inequality problem based on the De Saxce Formulation

driver :func:`variationalInequality_ExtraGradient()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000
* iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION] = SICONOS_VI_ERROR_EVALUATION_LIGHT_WITH_FULL_FINAL
* iparam[SICONOS_VI_IPARAM_ERROR_EVALUATION_FREQUENCY] (set but not used)
* iparam[SICONOS_VI_IPARAM_LINESEARCH_METHOD] = SICONOS_VI_LS_ARMIJO
  
  allowed values :

  * SICONOS_VI_LS_ARMIJO : Armijo rule with Khotbotov ratio (default)
  * SICONOS_VI_LS_SOLODOV : Armijo rule with Solodov.Tseng ratio
  * SICONOS_VI_LS_HANSUN : Armijo rule with Han.Sun ratio

* iparam[SICONOS_VI_IPARAM_ACTIVATE_UPDATE] = 0;
* iparam[SICONOS_VI_IPARAM_DECREASE_RHO] = 0;

* dparam[SICONOS_DPARAM_TOL] = 1e-3, in-out parameter
* dparam[SICONOS_VI_DPARAM_RHO] = -1., in-out parameter
* dparam[SICONOS_VI_DPARAM_LS_TAU] = 2/3
* dparam[SICONOS_VI_DPARAM_LS_TAUINV] = 3/2
* dparam[SICONOS_VI_DPARAM_LS_L] = 0.9
* dparam[SICONOS_VI_DPARAM_LS_LMIN] = 0.3


Fixed-point  projection (:enumerator:`SICONOS_VI_FPP`)
------------------------------------------------------

Fixed Point Projection solver for variational inequality problem based on the De Saxce Formulation.

driver: :func:`variationalInequality_FixedPointProjection()`

parameters: same as :enumerator:`SICONOS_VI_EG.`

Hyperplane  projection (:enumerator:`SICONOS_VI_HP`)
----------------------------------------------------

driver: :func:`variationalInequality_HyperplaneProjection()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 20000
  
* iparam[SICONOS_VI_IPARAM_LS_MAX_ITER] = 100
  
* dparam[SICONOS_DPARAM_TOL] = 1e-3
  
* dparam[SICONOS_VI_DPARAM_LS_TAU] = 1.0, tau
  
* dparam[SICONOS_VI_DPARAM_SIGMA] = 0.8, sigma
  
out :
  
* iparam[SICONOS_IPARAM_ITER_DONE] : number of iterations
    

SICONOS_VI_BOX_QI (:enumerator:`SICONOS_VI_BOX_QI`)
---------------------------------------------------

Solver using the merit function proposed by Qi for box-constrained Newton QI LSA

id: 

driver : :func:`variationalInequality_box_newton_QiLSA()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 1000
* iparam[SICONOS_IPARAM_PREALLOC] = 0  
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
  
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0
  
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0 (set but not used)
  
* iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH] = 1

* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16 
  
* dparam[SICONOS_DPARAM_TOL] = 1e-10
  
SICONOS_VI_BOX_AVI_LSA (:enumerator:`SICONOS_VI_BOX_AVI_LSA`)
-------------------------------------------------------------

driver : :func:`vi_box_AVI_LSA()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 100
* iparam[SICONOS_IPARAM_LSA_FORCE_ARCSEARCH] = 1
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS] = 0
  
* iparam[SICONOS_IPARAM_LSA_NONMONOTONE_LS_M] = 0 (set but not used)
* iparam[SICONOS_IPARAM_STOPPING_CRITERION] = SICONOS_STOPPING_CRITERION_USER_ROUTINE;
* dparam[SICONOS_DPARAM_TOL] = 1e-12
* dparam[SICONOS_DPARAM_LSA_ALPHA_MIN] = 1e-16 
  
internal solver : :enumerator:`SICONOS_RELAY_AVI_CAOFERRIS`

SICONOS_VI_BOX_PATH (:enumerator:`SICONOS_VI_BOX_PATH`)
-------------------------------------------------------

driver : :func:`vi_box_path()`

parameters:

* iparam[SICONOS_IPARAM_MAX_ITER] = 10000
  
* dparam[SICONOS_DPARAM_TOL] = 1e-12

