.. index:: single: Relay or box-constrained AVI problems
.. _doxid-_relay_problem:

Relay or box-constrained AVI problems
=====================================

.. _doxid-_relay_problem_1relayIntro:
.. rubric:: The problem:

Find :math:`(z,w)` such that:

.. math::

    \begin{equation*} \left\lbrace \begin{array}{l} w = M z + q\\ -w \in \mathcal{N}_{K}(z)\\ \end{array}, \right. \end{equation*}

where M is an ( :math:` n \times n ` )-matrix, q, z and w are n-dimensional vectors, K is the box defined by :math:`K=\{x\in\mathbb{R}^n \mid lb_i \leq x_i \leq ub_i, i = 1, ..., n \}` and :math:`\mathcal{N}_K(z)` is the normal cone to :math:`K` at :math:`z` .

The solvers and their parameters are described in :ref:`Relay Problems Solvers <doxid-_relay_solvers>` .

.. _doxid-_relay_problem_1relaySolversList:
.. rubric:: Available solvers:

The "direct" solvers are

* :ref:`relay_avi_caoferris() <doxid-_relay___solvers_8h_1a02bb6df7e1147154ba01c36d79cb3a11>` based on an algorithm by Cao and Ferris for AVI with a polytopic set :math:`K` .

* :ref:`relay_path() <doxid-_relay___solvers_8h_1a343139518641407ce13174b1ad360a32>` using the PATH solver

Using an LCP reformulation (splitting z in positive and negative part), we have the following available solvers:

* :ref:`relay_enum() <doxid-_relay___solvers_8h_1acb9e15dbc9ae90c1c02ac1f5a7937790>` which solves the LCP using the enumerative method

* :ref:`relay_lexicolemke() <doxid-_relay___solvers_8h_1af0ac56c5358a7148ce5c56eb68eabbc2>` which solves the LCP using Lemke's algorithm

(see the functions/solvers list in ``Relay_Solvers.h`` )

.. index:: single: Relay Problems Solvers
.. _doxid-_relay_solvers:

Relay Problems Solvers
======================

This page gives an overview of the available solvers for relay problems and their required parameters.

This page gives an overview of the available solvers for relay problems and their required parameters.

For each solver, the input argument are:

* a :ref:`RelayProblem <doxid-struct_relay_problem>`

* the unknowns (z,w)

* info, the termination value (0: convergence, >0 problem which depends on the solver)

* a SolverOptions structure, which handles iparam and dparam

.. _doxid-_relay_solvers_1relayENUM:
.. rubric:: Enumerative solver:

The relay problem is reformulated as a LCP and solved with the enumerative solver

function: :ref:`relay_enum() <doxid-_relay___solvers_8h_1acb9e15dbc9ae90c1c02ac1f5a7937790>`

parameters:

* dparam[0] (in): tolerance

* iparam[0] (in) : search for multiple solutions if 1

* iparam[1] (out) : key of the solution

* iparam[3] (in) : starting key values (seed)

* iparam[4] (in) : use DGELS (1) or DGESV (0).

.. _doxid-_relay_solvers_1relayPATH:
.. rubric:: PATH solver:

The relay problem is reformulated as a LCP and solved with the PATH solver

function: :ref:`relay_path() <doxid-_relay___solvers_8h_1a343139518641407ce13174b1ad360a32>`



* dparam[0] (in): tolerance

.. _doxid-_relay_solvers_1relayLEMKE:
.. rubric:: Lemke solver:

The relay problem is reformulated as a LCP and solved with Lemke's method

function: :ref:`relay_lexicolemke() <doxid-_relay___solvers_8h_1af0ac56c5358a7148ce5c56eb68eabbc2>`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

.. _doxid-_relay_solvers_1relayAVI_CaoFerris:
.. rubric:: CaoFerris solver:

The relay problem is reformulated as an AVI and solved with the solver proposed by Cao and Ferris

function: :ref:`relay_avi_caoferris() <doxid-_relay___solvers_8h_1a02bb6df7e1147154ba01c36d79cb3a11>`

parameters:

* iparam[0] (in): maximum number of iterations allowed

* iparam[1] (out): number of iterations processed

