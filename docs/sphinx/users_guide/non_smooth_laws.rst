.. _non_smooth_laws:

Non Smooth Laws
---------------

A non-smooth law is an object used to define the behavior of the systems involved in an Interaction, when a non-smooth event occurs.
For example, in the case of an impact, a Newton impact law will link the pre and post velocities at impact in something like "post-velocity = -e X pre-velocity".

Each non-smooth law is characterized by:

* a type (more or less the name of its class), i.e. what kind of law is required
* the size of vectors involved in the law
* some specific variables depending on its type.

Nonsmooth laws are defined in classes which name ends with "NSL". All of them are derived from an abstract class which defines a generic interface, *NonSmoothLaw*.

*Available classes*: NonSmoothLaw, ComplementarityConditionNSL, EqualityConditionNSL, MixedComplementarityConditionNSL, MultipleImpactNSL,  NewtonImpactNSL, NewtonImpactFrictionNSL, RelayNSL, NormalConeNSL.

.. image:: /figures/classNonSmoothLaw.*
	   
Complementarity Condition
^^^^^^^^^^^^^^^^^^^^^^^^^

nsLawSize: 1.
no specific parameters.

.. math::

   0 \leq y \perp \lambda\geq 0

Newton Impact
^^^^^^^^^^^^^

nsLawSize: 1.

parameter: *e*, the Newton normal coefficient of restitution.

.. math::
   
   if \ y(t)=0,\ 0 \leq \dot y(t^+) +e   \dot y(t^-) \perp \lambda\geq 0

Newton Impact-Friction
^^^^^^^^^^^^^^^^^^^^^^

nsLawSize: 2 or 3 (2 or 3 dimensional friction).

parameters: en, et (Newton impact normal and tangential coefficients) and mu, friction coefficient.

Newton Impact Law plus Coulomb Friction.

In this case, y components are in the following order:

* first relation, normal part
* first relation, tangential part
* ...
* relation n, normal part
* relation n, tangential part

and so on .

Note also that usually only normal part definition is required for y[0]. 

Relay
^^^^^

nsLawSize: 1.

parameters: c and d

.. math::

   \dot y &=0, \  d \leq \lambda \leq c \\
   \dot y &\geq 0, \ \lambda =  c \\
   \dot y &\leq 0, \ \lambda =  d \\
