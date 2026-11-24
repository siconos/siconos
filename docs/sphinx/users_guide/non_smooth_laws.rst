.. _non_smooth_laws:

Nonsmooth laws
--------------

A nonsmooth law is an object that defines the behaviour of the systems involved in an Interaction when a nonsmooth event occurs.
For example, for an impact, a Newton impact law links the pre- and post-impact velocities with a rule like "post-velocity = -e x pre-velocity".

Each non-smooth law is characterized by:

* a type (more or less the name of its class), i.e. what kind of law it is;
* the size of the vectors involved in the law;
* some specific variables depending on its type.

Nonsmooth laws are defined in classes whose names end with "NSL". All of them derive from an abstract class, *NonSmoothLaw*, which defines a generic interface.

*Available classes*: NonSmoothLaw, ComplementarityConditionNSL, EqualityConditionNSL, MixedComplementarityConditionNSL, MultipleImpactNSL,  NewtonImpactNSL, NewtonImpactFrictionNSL, RelayNSL, NormalConeNSL.

.. image:: /figures/classNonSmoothLaw.*

Complementarity Condition
^^^^^^^^^^^^^^^^^^^^^^^^^

nsLawSize: 1.

No specific parameters.

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

nsLawSize: 2 or 3 (2D or 3D friction).

parameters: en, et (Newton impact normal and tangential coefficients) and mu, the friction coefficient.

Newton impact law plus Coulomb friction.

In this case, y components are in the following order:

* first relation, normal part
* first relation, tangential part
* ...
* relation n, normal part
* relation n, tangential part

and so on.

Note also that usually only the normal part needs to be defined for y[0].

Relay
^^^^^

nsLawSize: 1.

parameters: c and d

.. math::

   \dot y &=0, \  d \leq \lambda \leq c \\
   \dot y &\geq 0, \ \lambda =  c \\
   \dot y &\leq 0, \ \lambda =  d \\
