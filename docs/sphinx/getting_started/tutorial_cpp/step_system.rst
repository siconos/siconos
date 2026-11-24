.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _step_system_example:

Gene regulatory networks
========================

.. highlight:: c++

This example describes the simulation of a gene regulatory network in Siconos.
For a complete description of the problem, see :cite:`Acary.DeJong.Brogliatio2014`.

The gene regulatory network is described by a piecewise linear model, in which the variables are concentrations of gene products.

We consider a two-gene system with the following dynamics:

.. math::

   \dot x_0 &=  -\gamma_0 x_0 + \kappa_0 S^+(x_1, \theta_1^0)S^-(x_0, \theta_0^1) \\
   \dot x_1 &=  -\gamma_1 x_1 + \kappa_1 S^+(x_0, \theta_0^0)S^-(x_1, \theta_1^1) \\

with the step functions:

.. math::
   :nowrap:

   \begin{equation}
   S^+(x_j, \theta_j^k) = \left\{\begin{array}{cc}
      1 & x_j > \theta_j^k \\ \ [0,1] & x_j = \theta_j^k \\
      0 & x_j < \theta_j^k\end{array}\right. \ \
   S^-(x_j, \theta_j^k) = \left\{\begin{array}{cc}
      0 & x_j > \theta_j^k \\ \ [0,1] & x_j = \theta_j^k \\
      1 & x_j < \theta_j^k\end{array}\right.
   \end{equation}

The :math:`\theta_j^k` are thresholds that control the inhibition or activation of the expression of a gene.
