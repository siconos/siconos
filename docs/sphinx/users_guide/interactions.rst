.. _interactions:

Interactions between dynamical systems
======================================

An Interaction is an "object" that defines the way some Dynamical Systems are linked, how they behave together. For example if you consider a set of rigid bodies, Interactions will define what will happen when contact between bodies occurs.

First of all, we introduce the set of all possible interactions, indeed all interactions declared by user, 

.. math::

   \mathcal{I}_0 = \{ \alpha \ / \ I_{\alpha} \in NSDS \}.
   
NSDS is the non-smooth dynamical system and :math:`I_{\alpha}` a single interaction.

An Interaction is applied to a set of Dynamical Systems, then the set :math:`\mathcal{DS}^{\alpha}=\{ d \ / \ ds_{d} \in I_{\alpha} \}`
is the set of all dynamical systems involved in :math:`I_{\alpha}`.
Finally, we denote :math:`\mathcal{DS}_{\alpha,\beta}` the set of dynamical systems that are involved in interactions :math:`\alpha` and :math:`\beta`;

.. math::
   
   \mathcal{DS}_{\alpha,\beta} =  \mathcal{DS}_{\alpha}\cap \mathcal{DS}_{\beta}

Considering an interaction :math:`I_{\alpha}`, upper case letters will be used to represent the concatenation of variables from the dynamical systems of :math:`\mathcal{DS}^{\alpha}`.
Note that Index :math:`\alpha` for variables or operators specific to :math:`I_{\alpha}`, will be omitted as soon as possible to lighten notations.

Then one get the following vectors of global coordinates, say X (or Q in the Lagrangian case):

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

Remember also that each DynamicalSystem as a variable called :math:`r` (or :math:`p` in the Lagrangian case), an input vector related to the non-smooth behavior or law, with:

.. math::

   r_d = \sum_{\beta \in \mathcal{I}^d} r_d^{\beta} \\
   \\
   \mathcal{I}^d = \{ \alpha \in \mathcal{I}_0 \ / \ d \in \mathcal{DS}_{\alpha} \}

Thus we define :math:`R` for the Interaction as:

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

*Warning: it is forbidden to mix first and second order Dynamical Systems in a single Interaction.*

An Interaction is characterized by some "local" variables, :math:`y` (also called output, :math:`R` being the input) and :math:`\lambda`. Both of them are "vector of vectors":

:math:`y[i]` is a vector that represents the derivative number :math:`i` of variable :math:`y` according to time. Each :math:`y[i]` or :math:`\lambda[i]` is a vector of size *interactionSize*.

Not that the number of saved derivatives depends on the problem type.

Then an Interaction proposes:
* a "Non Smooth Law" that links y and :math:`\lambda`
* a "Relation" between the local variables :math:`(y,\lambda)` and the global ones (those of the Dynamical Systems), :math:`(X,R)` (the constraints).

As an example consider again the case of a ball bouncing on the ground:
* the Interaction will include the two dynamical systems (ball and ground)
* the relation will consist in defining y as the distance between the ground and the ball and :math:`\lambda` as something like the reaction of the ground at contact
(Lagrangian multipliers indeed).
* the Non-Smooth law will "say": if there is contact, then the reaction is positive (y=0 then :math:`\lambda>0`) and if not, no reaction occurs(if :math:`\lambda=0, y>0`); and also that the velocity after contact is equal to the opposite of the one before contact multiplied by some "damping" coefficient (this corresponds to a Newton Impact Law, see below).

Definition and description of the different types of relations and non-smooth laws are presented in the sections below.


.. toctree::
   :maxdepth: 4

   relations
   non_smooth_laws
