.. _modal_moreau_jean:

Time integration of the dynamics - Exact scheme
===============================================

Time integration for second-order (Lagrangian) systems of the form:

.. math::
	   
   M\ddot q(t) + Kq(t) + C\dot q(t) = p \\
   q(t_0) = q_0 \\
   \dot q(t_0) = v_0


equivalent to

.. math::
   :label: LTIDS

   Mdv + Kqdt + cv^+dt = dr \\
   v(t) = \dot q^+(t) \\
   q(t_0) = q_0, \ \dot q(t_0) = v_0

with

.. math::
   :label: displ

   q(t) = q_0 + \int_{t_0}^t v^+(t)


:math:`v^+(t)` stands for right limit of v in t.

with the vectors :math:`q = [q_k], v = [v_k] k\in[0,ndof-1]`.

      
Let us integrate the dynamics over a time step, :math:`\Delta t = t^{i+1} - t^i`.

We consider two different schemes, the classical (in Siconos) 'Moreau-Jean' and 'Modal-Moreau-Jean' denoted respectively MJ and MMJ in the following.

Notations :

.. math::

   q_k(t^i) = q_k^i \\
   v_k(t^{i}) = v_k^i \\
   
In the following, we will use k for space (bottom) indices and i for time (top) indices.
   
MJ is based on a theta-scheme, for :math:`\theta \in [0,1]`

.. math::

   \int_{t^i}^{t^{i+1}} f(q,v,t)dt & \approx \Delta t\theta f(q^{i+1}, v^{i+1}, t^{i+1}) + \Delta t(1 -\theta)f(q^{i}, v^{i}, t^{i}) \\
                                   & \approx \Delta t\theta f^{i+1} + \Delta t(1 -\theta)f^{i} \\


				   
With MJM we consider diagonal stiffness and damping, 

.. math::
   :nowrap:

   \begin{eqnarray}
      K = diag(\omega_k^2) \\
      C = diag(2\sigma_k)
   \end{eqnarray}


:math:`\omega_k` and :math:`\sigma_k` being respectively the modal pulsation and the damping parameter (values to be taken from 2.2 and 2.3 in JSV paper).

Bilbao exact scheme writes:

.. math::

   \begin{array}{ccc}
   Kq &\approx \Theta Kq^i + \frac{(\mathcal{I}-\Theta)K}{2}(q^{i+1} + q^{i-1}) \\
   C\dot q &\approx \frac{1}{\Delta t}\Sigma^*(q^{i+1} - q^{i-1}) \\
   \end{array}

for :math:`\Theta = diag(\theta_k)` and :math:`\Sigma^* = diag(\sigma_k^*)` some diagonal matrices, with

.. math::
   
   \theta_{k} &= \frac{2}{\omega_k^2\Delta t^2} - \frac{A_k}{1+e_k-A_k}, \\
   \sigma^*_{k} &= \left(\frac{1}{\Delta t} + \frac{\omega_k^2\Delta t}{2} - \theta_k\frac{\omega_k^2\Delta t}{2} \right)\frac{1-e_k}{1+e_k} \\
   A_k &= e^{-\sigma_k\Delta t}\left(e^{\sqrt{\sigma_k^2 - \omega_k^2}\Delta t} + e^{-\sqrt{\sigma_k^2 - \omega_k^2}\Delta t}\right) \\
   e_k &= e^{-2\sigma_k\Delta t} \\

.. math::

   \begin{array}{c|c|c}
   Dynamics       & Moreau-Jean                       &         Modal-Moreau-Jean \\
   \int_{t^i}^{t^{i+1}} Mdv & \approx M(v^{i+1}-v^{i}) & \approx M(v^{i+1}-v^{i}) \\
   \int_{t^i}^{t^{i+1}} Kqdt & \approx \Delta t(\theta Kq^{i+1} + (1 - \theta) Kq^i) & \approx \Delta t\Theta Kq^i + \frac{\Delta t}{2}(\mathcal{I}-\Theta)K(q^{i+1} + q^{i-1}) \\
   \int_{t^i}^{t^{i+1}} Cvdt & \approx \Delta t(\theta Cv^{i+1} + (1 - \theta) Cv^i) & \approx \Sigma^*(q^{i+1} - q^{i-1})\\
   \int_{t^i}^{t^{i+1}} dr & \approx p^{i+1} & \approx p^{i+1} \\
    \end{array}

For MJ, this leads to

.. math::

   M(v^{i+1}-v^{i}) + \Delta t(\theta Kq^{i+1} + (1 - \theta) Kq^i) + \Delta t(\theta Cv^{i+1} + (1 - \theta) Cv^i) &= p^{i+1} \\

using :math:`q^{i+1} = q^i + \Delta t(\theta v^{i+1} + (1 - \theta) v^i)`, we get

.. math::
   
   [M + \Delta t^2\theta^2 K + \Delta t\theta C] (v^{i+1}-v^{i}) + \Delta tKq^i + (\Delta t^2\theta K + \Delta tC) v^i = p^{i+1} \\

And for MMJ:

.. math::

   M(v^{i+1}-v^{i}) + \Delta t\Theta Kq^i + \frac{\Delta t}{2}(\mathcal{I}-\Theta)K(q^{i+1} + q^{i-1}) +\Sigma^*(q^{i+1} - q^{i-1}) = p^{i+1}

With :math:`q^{i+1} = q^{i} + \Delta tv^{i+1}`, we get

.. math::
   
   q^{i+1} - q^{i-1} &= \Delta t(v^{i+1} + v^i) \\
   q^{i+1} + q^{i-1} &= 2q^i + \Delta t(v^{i+1} - v^i) \\

and

.. math::
   
   [M + \frac{\Delta t^2}{2}(\mathcal{I} - \Theta)K + \Delta t\Sigma^*] (v^{i+1}-v^{i}) + \Delta tKq^i + 2\Delta t \Sigma^* v^i = p^{i+1} \\
   

Both discretisations writes
   
.. math::
   
   W(v^{i+1}-v^{i}) = \tilde v_{free}^i + p^{i+1} \\
   or \\
   v^{i+1} = v^i_{free} + W^{-1}p^{i+1} \\

with

.. math::

   \begin{array}{c|c|c}
   & Moreau-Jean                       &         Modal-Moreau-Jean \\
   W & = M + \Delta t^2\theta^2 K + \Delta t\theta C & = M + \frac{\Delta t^2}{2}(\mathcal{I} - \Theta)K + \Delta t\Sigma^*\\
   v_{free}^{i} &= v^i - W^{-1}(\Delta tKq^i + (\Delta t^2\theta K + \Delta tC) v^i) & = v^i - W^{-1}(\Delta tKq^i + 2\Delta t \Sigma^* v^i) \\
   \end{array}

   
Notes, remarks, questions
^^^^^^^^^^^^^^^^^^^^^^^^^

* Quel nom pour "modal" Moreau-Jean? i.e. qui est à la source (ref?) du schéma de Bilbao?
* Vérif comportement de W quand :math:`\Delta t \rightarrow 0`

Non-smooth problem formulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. math::

   \dot y &= S_c \dot q \\
   P &= S_c^T\lambda

.. math::

   \dot y_c^{i+1} &= S_c \dot q^{i+1} \\
   P^{i+1} &= S_c^T\lambda^{i+1}


.. math::
   
   \dot y^{i+1} &= S_cv^{i} - S_cW^{-1}(\Delta tKq^i + 2\Delta t \Sigma^* v^i) + S_cW^{-1}S_c^T\lambda^{i+1} \\
           &= q_{LCP} + M_{LCP}\lambda^{i+1}

with

.. math::

   0 \leq \dot y^{i+1} \perp \lambda^{i+1} \geq 0
   
