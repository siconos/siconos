.. _osns_problems:

Nonsmooth problems formulation and solve
========================================

When dynamical systems and their interactions have been properly defined inside a model and its nonsmooth dynamical system,
a proper formulation for the latter must be chosen, associated to a nonsmooth solver.

Linear nonsmooth problems
-------------------------

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n}

where :math:`w \in R^{n}, z \in R^{n}` are the unknowns.

* Linear Complementarity Problems (:doxysiconos:`LCP`)

.. math::

   w =  q + M z, M \in R^{n \times n }, q \in R^{n} \\
   w \geq 0, z \geq 0,  z^{T} w =0

* Mixed Linear Complementarity Problems (:doxysiconos:`MLCP`)

.. math::
  
   0 =  Au + Cv + a\\
   z =  Du + Bv + b\\
   v \geq 0, z \geq 0,  z^{T} v =0

where

.. math::
   u \in R^{n},  v \in R^{m}, z \in R^{m} \ the \ unknowns, \\
   a \in R^{n}, b \in R^{m} \\
   A \in R^{n \times n }, B \in R^{m \times m }\\
   C \in R^{n \times m }, D \in R^{m \times n }
   
* :doxysiconos:`Relay`
* :doxysiconos:`Equality`
* :doxysiconos:`AVI`
* 2D or 3D friction contact problem :doxysiconos:`FrictionContact`

.. math::

  velocity =  q + M reaction \\
  velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0

and a Coulomb friction law.
With :math:`velocity \in R^{n}, reaction \in R^{n}` the unknowns,
and :math:`M \in R^{n \times n }, q \in R^{n}`

  
* multiple-impact problem (:doxysiconos:`OSNSMultipleImpact`)
* primal friction contact problems (:doxysiconos:`GlobalFrictionContact`)

.. math::

   M velocity =  q +  reaction \\
   localVelocity = H^T velocity + tildeLocalVelocity\\
   reaction = H localReaction \\

and :math:`localVelocity,localReaction` belongs to the Coulomb friction law with unilateral contact.

With :math:`velocity \in R^{n}, reaction \in R^{n}, localVelocity \in R^{m}, localReaction \in R^{m}` the unknowns,
:math:`M \in R^{n \times n }, q \in R^{n}`. :math:`tildeLocalVelocity \in R^{m}` is the modified local velocity (:math:`e U_{N,k}`), :math:`M \in R^{n \times n }, q \in R^{n}, H \in R^{n \times m }`.
   
* Generic mechanical problem (:doxysiconos:`GenericMechanical`)
  
Complete problem with bilateral equality, complementarity, impact and friction.

