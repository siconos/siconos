.. SPDX-FileCopyrightText: 2026 INRIA
..
.. SPDX-License-Identifier: Apache-2.0

.. _interactions:

Interactions between dynamical systems
======================================

An Interaction is an "object" that defines how some Dynamical Systems are linked and how they behave together. For example, for a set of rigid bodies, Interactions define what happens when bodies come into contact.

First, we introduce the set of all possible interactions, i.e. all the interactions declared by the user,

.. math::

   \mathcal{I}_0 = \{ \alpha \ / \ I_{\alpha} \in NSDS \}.

NSDS is the nonsmooth dynamical system and :math:`I_{\alpha}` a single interaction.

An Interaction applies to a set of Dynamical Systems; the set :math:`\mathcal{DS}^{\alpha}=\{ d \ / \ ds_{d} \in I_{\alpha} \}`
is the set of all dynamical systems involved in :math:`I_{\alpha}`.
Finally, we denote by :math:`\mathcal{DS}_{\alpha,\beta}` the set of dynamical systems involved in both interactions :math:`\alpha` and :math:`\beta`:

.. math::

   \mathcal{DS}_{\alpha,\beta} =  \mathcal{DS}_{\alpha}\cap \mathcal{DS}_{\beta}

For an interaction :math:`I_{\alpha}`, upper-case letters denote the concatenation of the variables of the dynamical systems of :math:`\mathcal{DS}^{\alpha}`.
Note that the index :math:`\alpha` of variables or operators specific to :math:`I_{\alpha}` is omitted whenever possible to lighten the notation.

This gives the following vectors of global coordinates, X (or Q in the Lagrangian case):

.. math::

   X=\left[\begin{array}{c}
   x_0 \\
   x_1 \\
   ...
   \end{array}\right], \ or \ Q=
   \left[\begin{array}{c}
   q_0\\
   q_1 \\
   ...
   \end{array}\right]

:math:`x_i \ (or \ q_i)` being the vectors of global coordinates of the Dynamical Systems involved in the Interaction.

Remember also that each DynamicalSystem has a variable called :math:`r` (or :math:`p` in the Lagrangian case), an input vector related to the nonsmooth behaviour or law, with:

.. math::

   r_d = \sum_{\beta \in \mathcal{I}^d} r_d^{\beta} \\
   \\
   \mathcal{I}^d = \{ \alpha \in \mathcal{I}_0 \ / \ d \in \mathcal{DS}_{\alpha} \}

We thus define :math:`R` for the Interaction as:

.. math::

   R^{\alpha}=\left[\begin{array}{c}
   r^{\alpha}_0 \\
   r^{\alpha}_1 \\
   ...
   \end{array}\right], \ or \ P^{\alpha}=
   \left[\begin{array}{c}
   p^{\alpha}_0\\
   p^{\alpha}_1 \\
   ...
   \end{array}\right]

*Warning: first- and second-order Dynamical Systems cannot be mixed in a single Interaction.*

An Interaction is characterized by "local" variables, :math:`y` (also called output, :math:`R` being the input) and :math:`\lambda`. Both are "vectors of vectors":

:math:`y[i]` is the :math:`i`-th time derivative of :math:`y`. Each :math:`y[i]` or :math:`\lambda[i]` is a vector of size *interactionSize*.

Note that the number of saved derivatives depends on the problem type.

An Interaction then provides:

* a "Nonsmooth Law" that links y and :math:`\lambda`;
* a "Relation" between the local variables :math:`(y,\lambda)` and the global ones (those of the Dynamical Systems), :math:`(X,R)` (the constraints).

As an example, consider again a ball bouncing on the ground:

* the Interaction includes the two dynamical systems (ball and ground);
* the relation defines y as the distance between the ground and the ball, and :math:`\lambda` as the reaction of the ground at contact
  (a Lagrange multiplier);
* the nonsmooth law "says" that if there is contact, the reaction is nonnegative (y = 0 and :math:`\lambda \geq 0`), and otherwise there is no reaction (y > 0 and :math:`\lambda = 0`); and also that the velocity after the contact is the opposite of the velocity before the contact multiplied by a restitution coefficient (this is a Newton impact law, see below).

The different types of relations and nonsmooth laws are defined and described in the sections below.


.. toctree::
   :maxdepth: 4

   relations
   non_smooth_laws
