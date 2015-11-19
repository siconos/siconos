.. _event_capturing:

Event-Capturing schemes
=======================

General Principle
-----------------

Roughtly speaking, the event-capturing, or time-stepping, method consists in the time-discretisation of the whole system (dynamics + relations + non-smooth laws), 
leading to a so-called one-step non smooth problem (OSNSP) solved at each time step.

Indeed, the main stages of the process are:

* integrate the dynamics without constraints, to get some "free" solutions
* formalize and solve a OSNSP (a LCP for example)
* update the dynamics with the OSNSP solutions to get the full state.

The whole Time-Stepping process is described in the sections below, starting from first order systems and then going on with second order, Lagrangian, systems.

In the following sections, the systems are integrated over a time step :math:`[t_i,t_{i+1}]`  of length :math:`h`.
The approximation of any function :math:`F(t,...)` at the time :math:`t_i` is denoted :math:`F_i`.
Note that in the relations writing, we use upper case letters for all variables related to DynamicalSystem objects.
and thus, :math:`X , Q, \ldots` are concatenation of :math:`x, q,\ldots` of the dynamical systems variables concerned by the relation.

Note also that most of the time, the parameter :math:`z` will be omitted, to lighten notations. 

First order systems
-------------------

Time Discretisation of the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First Order Non Linear Systems
""""""""""""""""""""""""""""""

.. math::
   M\dot x(t) &=& f(x,t,z) + r   \\
   x(t_0) &=& x_0

with :math:`r = r^d = \sum_{\alpha} r^{\alpha}, \alpha \in I_d`, :math:`I_d` being the set of all relations in which the current dynamical system, number :math:`d`, is involved. 
In the following, the index "d" will be omitted to lighten notations.

The integration of the ODE over a time step :math:`[t_i,t_{i+1}]`  of length :math:`h`  is :

.. math::
   
   M\int_{t_i}^{t_{i+1}}\dot x\,dt = \int_{t_i}^{t_{i+1}} f(t,x,z)dt + \int_{t_i}^{t_{i+1}}r\,dt   

The left-hand term is :math:`M(x(t_{i+1})-x(t_i)) \approx M(x_{i+1} - x_i)` .

Right-hand terms are approximated with a :math:`\theta`-method:

.. math::
   \int_{t_i}^{t_{i+1}} f(t,x,z)dt &\approx& h \theta f(t_{i+1},x_{i+1},z) + h (1-\theta) f(t_i,x_i,z) \\
   &\approx& h \theta f_{i+1} + h (1-\theta) f_i

and the third integral is approximated with:

.. math::
   \int_{t_i}^{t_{i+1}}r\,dt \approx h r(t_{i+1}) \approx hr_{i+1}

Then, we get the following "residu"

.. math::
   
   \mathcal R(x_{i+1}) &=& M(x_{i+1}-x_i) - h \theta f_{i+1} - h (1-\theta) f_{i} - hr_{i+1} = 0 \\
	     &=& \mathcal R^{free}(x_{i+1}) - hr_{i+1}

<em> Note: We introduce the "free" notation for terms related to the smooth part of the system. </em>

We apply a Newton method to solve :math:`\mathcal R(x_{k+1}) = 0`. The gradient of the residu according to :math:`x` is:

.. math::
   \nabla_{x}\mathcal R(x) = M - h \theta\cdot\nabla_{x}f(t,x)

And we get (index k corresponds to the Newton iteration number):

.. math::
   W_{i+1}^k\cdot (x_{i+1}^{k+1} - x_{i+1}^k) = - \mathcal R(x_{i+1}^k)

with

.. math::
   W_{i+1}^k = M - h \theta\left[\nabla_{x}f\right](t_{i+1},x_{i+1}^k)

If we assume that :math:`W_{i+1}^k` is invertible, we get the solution at Newton iteration k+1:

.. math::
   x_{i+1}^{k+1} &=& x_{i+1}^k - (W_{i+1}^k)^{-1}\mathcal R^{free}(x_{i+1}^{k}) + h(W_{i+1}^k)^{-1}r_{i+1}^{k+1} \\
	      &=& x^{free,k}_{i+1} + h(W_{i+1}^k)^{-1}r_{i+1}^{k+1}

First Order Linear Systems
""""""""""""""""""""""""""

.. math::

   M\dot x(t) &=& A(t,z)x(t) + b(t) + r   \\
   x(t_0) &=& x_0

For the integration of the ODE over a time-step, we proceed as in the previous section for non-linear systems to get:

.. math::

   \mathcal R(x_{i+1}) &=& M(x_{i+1}-x_i) - h \theta(A_{i+1}x_{i+1} + b_{i+1})- h (1-\theta)(A_{i}x_i + b_i) -  hr_{i+1} = 0 \\
   or \\
   (M - h\theta A_{i+1}) x_{i+1} = (M + h (1-\theta)A_{i})\cdot x_i + h\theta(b_{i+1}-b_i) + hb_i +  hr_{i+1}  

We denote:

.. math::

   W_{i+1} =  (M - h\theta A_{i+1})

and assuming it is invertible, we get:

.. math::

   x_{i+1} &=& W_{i+1}^{-1}\left[(M + h (1-\theta)A_{i})\cdot x_i + h\theta(b_{i+1}-b_i) + hb_i\right] +  hW_{i+1}^{-1}r_{i+1}  \\
   &=& x^{free}_{i+1}  +  hW_{i+1}^{-1}r_{i+1} 

First Order Linear Systems with time invariant coefficients
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. math::

   M\dot x(t) &=& Ax(t) + b + r   \\
   x(t_0) &=& x_0

Using the results of the previous section, the discretisation is straightforward:

.. math::
   x_{i+1} &=& x_i + h W^{-1}(A x_i + b) +  hW^{-1}r_{i+1} \\
   &=& x^{free}_{i}  +  hW^{-1}r_{i+1} 

with a W that do not depend on time:

.. math::
   W =  (M - h\theta A)

Time discretization of the relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the following, :math:`R` represents the concatenation of all :math:`r^{\alpha}` vectors for the DS involved in the present relation.

First Order (non-linear) Relations
""""""""""""""""""""""""""""""""""

.. math::

   y &=& h(X,t,\lambda,Z)\\
   R &=& g(X,t,\lambda,Z)

Then, for the iteration :math:`k+1` of the Newton process, we get:

.. math::
   y_{i+1}^{k+1} &=& h(X_{i+1}^{k+1},t_{i+1},\lambda_{i+1}^{k+1})\\
   R_{i+1}^{k+1} &=& g(X_{i+1}^{k+1},t_{i+1},\lambda_{i+1}^{k+1})

These constraints are linearized around state :math:`(X_{i+1}^{k+1},\lambda_{i+1}^{k+1})`:

.. math::

   y_{i+1}^{k+1} &=& y_{i+1}^k - H_0(S_{i+1}^k)X_{i+1}^{k} - H_1(S_{i+1}^k)\lambda_{i+1}^{k} + H_0(S_{i+1}^k)X_{i+1}^{k+1} + H_1(S_{i+1}^k)\lambda_{i+1}^{k+1}  \\
   \\
   R_{i+1}^{k+1} &=& R_{i+1}^k - G_0(S_{i+1}^k)X_{i+1}^{k} - G_1(S_{i+1}^k)\lambda_{i+1}^{k} + G_0(S_{i+1}^k)X_{i+1}^{k+1} + G_1(S_{i+1}^k)\lambda_{i+1}^{k+1} 

Where :math:`S_{i+1}^k` stands for :math:`(X_{i+1}^{k},t_{i+1},\lambda_{i+1}^{k})` and

.. math::
   H_0(X,t,\lambda)=\nabla_X h(X,t,\lambda)&, &  H_1(X,t,\lambda)=\nabla_{\lambda} h(X,t,\lambda) \\
   \\
   G_0(X,t,\lambda)=\nabla_X g(X,t,\lambda)&, &  G_1(X,t,\lambda)=\nabla_{\lambda} g(X,t,\lambda) 

In the case where :

.. math::
   x_{i+1}^{k+1} = x^{free,k}_{i+1} + (w_{i+1}^k)^{-1}r_{i+1}^{k+1}

We can write

.. math::
   X_{i+1}^{k+1} = X^{free,k}_{i+1} + (W_{i+1}^k)^{-1}R_{i+1}^{k+1}

where :math:`(W_{i+1}^k)^{-1}`, is a diagonal block matrix holding the :math:`(w_{i+1}^k)^{-1}`,
then, if there is one and only one interaction we have:

.. math::

   (1-(W_{i+1}^k)^{-1}G_{0,i+1}^k) X_{i+1}^{k+1} = X_{i+1}^{free,k} + (W_{i+1}^k)^{-1} (R_{i+1}^k - G_{0,i+1}^k X_{i+1}^k - G_{1,i+1}^k \lambda_{i+1}^k + G_{1,i+1}^k \lambda_{i+1}^{k+1})

and finally: 

.. math::
   y_{i+1}^{k+1} &=& M_{lcp}\lambda_{i+1}^{k+1} + q_{lcp} \\
   M_{lcp} &=& H_{1,i+1}^k + H_{0,i+1}^k (1-(W_{i+1}^k)^{-1} G_{0,i+1}^k)^{-1} (W_{i+1}^k)^{-1} G_{1,i+1}^k \\
   q_{lcp} &=& y_{i+1} -H_{0,i+1}^k X_{i+1}^k - H_{1,i+1}^k \lambda_{i+1}^k + H_{0,i+1}^k (1-(W_{i+1}^k)^{-1} G_{0,i+1}^k)^{-1}
   [X_{i+1}^{free,k} + (W_{i+1}^k)^{-1} (R_{i+1}^k - G_{0,i+1}^k X_{i+1}^k - G_{1,i+1}^k \lambda_{i+1}^k)]

First Order Linear Relations
""""""""""""""""""""""""""""

.. math::

 y &=& C(t,Z)X(t) + F(t,Z)Z + D(t,Z)\lambda + e(t,Z) \\
 R &=& B(t,Z) \lambda

<em> Note: for time-invariant relations, B, C, F, D and e are constant vectors and matrices </em>

The Time discretization of the relations is fully implicit and may be written as :

.. math::
   y_{i+1} &=& C(t_{i+1})X_{i+1} + D(t_{i+1})\lambda_{i+1} + e(t_{i+1}) + F(t_{i+1})Z \\	
   \\
   R_{i+1} &=& B(t_{i+1})\lambda_{i+1}

Discretisation of the non-smooth law
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Complementarity Condition
"""""""""""""""""""""""""

The complementarity condition writes:

.. math::
   0 \leq y \, &\perp& \, \lambda \geq 0 

and the discretisation is straightforward:

.. math::
   0 \leq y_{i+1} \, &\perp& \, \lambda_{i+1} \geq 0 

Lagrangian systems
------------------

Time Discretisation of the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lagrangian (second order) Non Linear Systems
""""""""""""""""""""""""""""""""""""""""""""

We provide in the following sections a time discretization method of the Lagrangian dynamical systems, consistent with the non smooth character of the solution.

.. math::
   
   M(q(t),z) dv &=& f_L(t,v^+(t), q(t), z)dt + dr \\
   v^+(t) &=& \dot q^+(t) \\
   q(t_0) &=& q_0 \\
   \dot q(t_0^-) &=& v_0 

with 

.. math::
   q(t) = q_0 + \int_{t_0}^t v^+(t)dt

<em> We recall that :math:`v^+(t)` means :math:`v(t^+)` ie right limit of :math:`v` in t. </em>

Left hand side is discretised by assuming that:

.. math::
   \int_{t_i}^{t_{i+1}} M(q(t),z)dv \approx M(q*,z)(v_{i+1}-v_{i}) 

As for first order non-linear systems, we use a :math:`\theta`-method to integrate the other terms, and obtain:

.. math::

   \int_{t_i}^{t_{i+1}} f_L(t, v^+(t), q(t), z) dt \approx  h\theta f_L(t_{i+1}, v_{i+1}, q_{i+1}, z) + h(1-\theta) f_L(t_{i}, v_{i}, q_{i}, z)

and for the last term, we set a new variable :math:`p_{i+1}` such that:

.. math::
   
   \int_{t_i}^{t_{i+1}} dr \approx p_{i+1}

Finally the full system discretisation results in:

.. math::
   
   \mathcal R(v_{i+1}, q_{i+1}) &=& M(q*,z)(v_{i+1}-v_{i}) - h\theta {f_L}_{i+1} - h(1-\theta) {f_L}_{i} - p_{i+1} = 0 \\	
   &=& \mathcal R^{free}(v_{i+1},q_{i+1}) - p_{i+1} 

The "free" notation still stands for terms related to the smooth part of the system. 
The displacement is integrated through the velocity with :

.. math::
   
   q_{i+1} &\approx& q_i + h\theta v_{i+1} + h(1 - \theta)v_{i}

Substituing this into the residu leads to a function depending only on :math:`v_{i+1}`, since state "i" and "k" are supposed to be known.\n
A Newton method will be applied to solve :math:`\mathcal R(v_{i+1}) = 0`.

That requires to compute the gradient of the residu;
assuming that the mass matrix evolves slowly with the configuration in a single time step, we get:

.. math::
   \nabla_{v_{i+1}}\left[M(q*,z)(v_{i+1}-v_{i})\right] \approx M(q^{*},z)

and denoting:

.. math::

   C_t(t,v,q)=-\left[\frac{\partial{f_L(t,v,q)}}{\partial{v}}\right] \\
   \\
   K_t(t,v,q)=-\left[\frac{\partial{f_L(t,v,q)}}{\partial{q}}\right]

we get (index k corresponds to the Newton iteration number):

.. math::
   W(t_{i+1}^k,v_{i+1}^k,q_{i+1}^k)\cdot (v_{i+1}^{k+1}-v_{i+1}^k) = - \mathcal R(v_{i+1}^k)

with

.. math::
   W(t,v,q) = M(q*,z) + h\theta C_t(t,v,q) + h^2\theta^2 K_t(t,v,q)

As an approximation for :math:`q^*`, we choose:

.. math::
   q^* &\approx& (1-\gamma) q_i  + \gamma q_{i+1}^k \\
   &\approx & q_i + h\gamma\left[ (1-\theta) v_i + \theta v_{i+1}^k\right]

with :math:`\gamma \in \left[0,1\right]`.
Moreover, if :math:`M` is evaluated at the first step of the Newton iteration, with :math:`v_{i+1}^0 = v_i`, we get:

.. math::
   M(q^*) \approx M(q_i + h\gamma v_i)

Finally, if :math:`W` is invertible, the solution at iteration k+1 is given by, 

.. math::
   v_{i+1}^{k+1} &=& v_{i+1}^k - (W_{i+1}^k)^{-1} \mathcal R^{free}(v_{i+1}^k) + (W_{i+1}^k)^{-1} p_{i+1}^{k+1} \\
   &=& v^{free,k}_{i+1} + (W_{i+1}^k)^{-1} p_{i+1}^{k+1}
\f}

Lagrangian (second order) Linear Systems with Time Invariant coefficients
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. math::

   M dv + Cv^+(t) + K q(t) &=& F_{ext}(t,z) + p \\
   q(t_0) &=& q0 \\
   \dot q(t_0^-) &=& v_0 

Proceeding in the same way as in the previous section, with :math:`M` constant and :math:`f_L(t,v^+(t), q(t), z) = F_{ext}(t) - Cv^+(t) - Kq(t)`, integration is straightforward:

.. math::

   \mathcal R(v_{i+1}, q_{i+1}) &=& M(v_{i+1}-v_{i}) - h\theta\left[ F_{ext}(t_{i+1}) - Cv_{i+1} - K q_{i+1}\right] - h(1-\theta)\left[ F_{ext}(t_{i}) - Cv_{i} - K q_{i}\right]  - p_{i+1} = 0  

Using the displacement integration through the velocity,


.. math::
   q_{i+1} = q_{i} +  h\left[\theta v_{i+1}+(1-\theta) v_{i}  \right]\\

we get:

.. math::
   W(v_{i+1}-v_{i}) &=& (- hC - h^2\theta  K )v_{i} - h K q_{i} +  h\left[\theta  F_{ext}(t_{i+1})+(1-\theta)  F_{ext}(t_{i})  \right] + p_{i+1} 

with :math:`W` a constant matrix:

.. math::
   W = \left[M + h\theta C + h^2 \theta^2 K \right]

and if :math:`W` is invertible,

.. math::
   
   v_{i+1} &=& v_{i} + W^{-1}\left[(- hC - h^2\theta  K )v_{i} - h K q_{i}+  h\theta  F_{ext}(t_{i+1})+h(1-\theta)  F_{ext}(t_{i}) \right] + W^{-1} p_{i+1} \\
   &=& v^{free}_i + W^{-1} p_{i+1} 

The free velocity :math:`v^{free} ` correponds to the velocity of the system without any constraints.

Time discretization of the relations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Lagrangian Scleronomous Relations
"""""""""""""""""""""""""""""""""

.. math::
   y &=& h(Q,Z) \\
   \dot y &=& G_0(Q,Z)V \\
   P &=& G_0^t(Q,Z)\lambda  

with

.. math::
   G_0(Q) &=& \nabla_Qh(Q) \\

From now on, to lighten the notations, the parameter :math:`Z` will omitted.

Considering the Newton process introduced above for Lagrangian non linear systems, the constraints write:

.. math::
   \dot y_{i+1}^{k+1} = G_0(Q_{i+1}^{k+1}))V_{i+1}^{k+1} \\
   P_{i+1}^{k+1} = G_0^t(Q_{i+1}^{k+1}))\lambda_{i+1}^{k+1}

To evaluate :math:`G_0` we still use the prediction :math:`Q^*` defined in the previous section:
.. math::

   Q^*( V_{i+1}^{k+1}) = Q_i + h\gamma \left[ (1-\theta) V_i + \theta  V_{i+1}^{k+1} \right]

Then we get:

.. math::

   \dot y_{i+1}^{k+1} = G_0(Q^*(V_{i+1}^{k+1}))V_{i+1}^{k+1} \\
   \\
   P_{i+1}^{k+1} = G_0^t(Q^*(V_{i+1}^{k+1}))\lambda_{i+1}^{k+1} 

These constraints are linearized around the point :math:`V_{i+1}^{k}` and we neglect the second order terms in the computation of the jacobians.
It leads to: 

.. math::

   \dot y_{i+1}^{k+1} = G_0(Q^*(V_{i+1}^k))V_{i+1}^{k+1} \\
   \\
   P_{i+1}^{k+1} = G_0^t(Q^*(V_{i+1}^k))\lambda_{i+1}^{k+1} 

As for the evaluation of the mass, the prediction of the position, :math:`Q^*` can be evaluated at the first iteration of the Newton process,

.. math::
   Q^*(V_{i+1}^0) =  Q_i + h\gamma V_i

Lagrangian Rheonomous Relations
"""""""""""""""""""""""""""""""

.. math::
   
   y &=& h(Q,t) \\
   \dot y &=& G_0(Q,t)V + G_1(Q,t) \\
   P &=& G_0^t(Q,t)\lambda  \\
   with\\
   G_0(Q,t) &=& \nabla_Qh(Q,t) \\
   G_1(Q,t) &=& \frac{\partial{h(Q,t)}}{\partial{t}} \\

As for scleronomous relations, we get:

.. math::
   
   \dot y_{i+1}^{k+1} &=& G_0(Q^*(V_{i+1}^k),t_{i+1})V_{i+1}^{k+1} +  G_1(Q^*(V_{i+1}^k, t_{i+1})) \\
   \\
   P_{i+1}^{k+1} &=& G_0^t(Q^*(V_{i+1}^k),t_{i+1})\lambda_{i+1}^{k+1} 

Lagrangian Compliant Relations
""""""""""""""""""""""""""""""

.. math::

   y &=& h(Q,\lambda(t)) \\
   \dot y &=& G_0(Q,\lambda(t))V + G_1(Q,\lambda(t))\dot\lambda \\
   P &=& G_0^t(Q,\lambda(t))\lambda
   with\\
   G_0(Q,\lambda(t)) &=& \nabla_Qh(Q,\lambda(t)) \\
   G_1(Q,\lambda(t)) &=& \nabla_\lambda h(Q,\lambda(t)) \\

Following the same process as in the paragraph above, it comes: 

.. math::
   \dot y_{i+1}^{k+1} &=& G_0(Q^*(V_{i+1}^k),\lambda_{i+1}^k)V_{i+1}^{k+1} +  G_1(Q^*(V_{i+1}^k, \lambda_{i+1}^k))\lambda_{i+1}^{k+1} \\
   \\
   P_{i+1}^{k+1} &=& G_0^t(Q^*(V_{i+1}^k),\lambda_{i+1}^k)\lambda_{i+1}^{k+1} 

Lagrangian Linear Relations
"""""""""""""""""""""""""""

.. math::

   y &=& HQ + D\lambda + FZ + b \\
   \dot y &=& HV + D\lambda \\
   P &=& H^t\lambda  

The discretisation is straightforward:

.. math::

   \dot y_{i+1} &=& HV_{i+1} + D\lambda_{i+1}
   \\
   P_{i+1} &=& H^t\lambda_{i+1}

Time discretization of the Non Smooth laws
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A natural way of discretizing the unilateral constraint  leads to the following implicit discretization :

.. math::

   0 \leq y_{i+1} \perp  \lambda_{i+1}  \geq 0

In the Moreau's time--stepping, we use a reformulation of the unilateral constraints in terms of velocity:

.. math::

   If y(t) =0, \ then \ 0 \leq \dot y \perp  \lambda  \geq 0

which leads to the following discretisation :

.. math::
   If \ y^{p} \leq 0, \ then \ 0 \leq \dot y_{i+1} \perp  \lambda_{i+1}  \geq 0

 where :math:`y^{p}` is a prediction of the position at time :math:`t_{i+1}`, for instance,    :math:`y^{p} = y_{i} + \frac{h}{2}  \dot y_i`.

If we want to introduce now the Newton impact law, we consider an equivalent velocity defined by

.. math::
   \dot y^{e}_{i+1} = \dot y_{i+1} + e \dot y_{i}

and we apply the constraints directly on this velocity :

.. math::
   
   If \ y^{p} \leq 0, \ then \ 0 \leq \dot y^{e}_{i+1} \perp  \lambda_{i+1}  \geq 0

Summary of the time discretized equations 
"""""""""""""""""""""""""""""""""""""""""

* \b First \b Order \b Non \b Linear \b Systems:

.. math::

   x_{i+1}^{k+1} &=& x^{free,k}_{i+1} + h(W_{i+1}^k)^{-1}r_{i+1}^{k+1} \\
   W_{i+1}^k &=& M - h \theta\cdot\nabla_{x}f(x_{i+1}^k,t_{i+1}) \\
   x^{free,k}_{i+1}  &=& x_{i+1}^k - (W_{i+1}^k)^{-1}\mathcal R^{free}(x_{i+1}^{k}) \\
   \mathcal R^{free}(x_{i+1}^{k}) &=& M(x_{i+1}^k-x_i) - h \theta f(x_{i+1}^k,t_{i+1}) - h (1-\theta) f(x_{i},t_i)

* \b First \b Order \b Linear \b Systems:

.. math::
   
   x_{i+1} &=& x^{free}_{i+1}  +  hW_{i+1}^{-1}r_{i+1} \\
   W_{i+1} &=& (M - h\theta A_{i+1}) \\
   x^{free}_{i+1} &=& W_{i+1}^{-1}\left[(M + h (1-\theta)A_{i})\cdot x_i + h\theta(b_{i+1}-b_i) + hb_i\right]

* <b> First Order Linear Systems with time-invariant coefficients: </b>

.. math::

   x_{i+1} &=& x^{free}_{i}  +  hW^{-1}r_{i+1} \\
   W &=& (M - h\theta A) \\
   x^{free}_i &=&  x_i + h W^{-1}(A x_i + b)

* <b> First Order Non Linear Relations </b>

.. math::

   y_{i+1}^{k+1} &=& y_{i+1}^k - H_0(S_{i+1}^k)X_{i+1}^{k} - H_1(S_{i+1}^k)\lambda_{i+1}^{k} + H_0(S_{i+1}^k)X_{i+1}^{k+1} + H_1(S_{i+1}^k)\lambda_{i+1}^{k+1}  \\
   \\
   R_{i+1}^{k+1} &=& R_{i+1}^k - G_0(S_{i+1}^k)X_{i+1}^{k} - G_1(S_{i+1}^k)\lambda_{i+1}^{k} + G_0(S_{i+1}^k)X_{i+1}^{k+1} + G_1(S_{i+1}^k)\lambda_{i+1}^{k+1} \\
   \\
   S_{i+1}^k \ for \ (X_{i+1}^{k},t_{i+1},\lambda_{i+1}^{k}) \\
   \\
   H_0(X,t,\lambda)=\nabla_X h(X,t,\lambda)&, &  H_1(X,t,\lambda)=\nabla_{\lambda} h(X,t,\lambda) \\
   \\
   G_0(X,t,\lambda)=\nabla_X g(X,t,\lambda)&, &  G_1(X,t,\lambda)=\nabla_{\lambda} g(X,t,\lambda) 

* <b> First Order Linear Relations </b>

.. math::
   
   y_{i+1} &=& C(t_{i+1})X_{i+1} + D(t_{i+1})\lambda_{i+1} + e(t_{i+1}) + F(t_{i+1})Z \\	
   R_{i+1} &=& B(t_{i+1})\lambda_{i+1}

* <b> Lagrangian Non Linear Case </b>

.. math::

   v_{i+1}^{k+1} &=& v^{free,k}_{i+1} + (W_{i+1}^k)^{-1} p_{i+1}^{k+1} \\
   q_{i+1}^{k+1} &=& q_i + h\theta v_{i+1}^{k+1} + h(1 - \theta)v_{i} \\
   v^{free,k}_{i+1} &=& v_{i+1}^k - (W_{i+1}^k)^{-1} \mathcal R^{free}(v_{i+1}^k) \\
   \mathcal R^{free}(v_{i+1}^k) &=& M(q*)(v_{i+1}^k-v_{i}) - h\theta f_L(t_{i+1},v_{i+1}^k,q_{i+1}^k) - h(1-\theta) f_L(t_i,v_i,q_i) \\
   W_{i+1}^k &=& M(q*) + h\theta C_t(t_{i+1},v_{i+1}^k,q_{i+1}^k) + h^2\theta^2 K_t(t_{i+1},v_{i+1}^k,q_{i+1}^k) \\
   q^* &=& q_i + h\gamma v_i \\
   C_t(t,v,q)&=&-\left[\frac{\partial{f_L(t,v,q)}}{\partial{v}}\right] \\
   K_t(t,v,q)&=&-\left[\frac{\partial{f_L(t,v,q)}}{\partial{q}}\right]

* <b> Lagrangian Linear and Time--Invariant Coefficients Case </b>

.. math::
   
   v_{i+1} &=& v^{free}_i + W^{-1} p_{i+1} \\
   q_{i+1} &=& q_{i} +  h\left[\theta v_{i+1}+(1-\theta) v_{i}  \right]\\
   v^{free}_i &=& v_{i} + W^{-1}\left[(- hC - h^2\theta  K )v_{i} - h K q_{i}+  h\theta  F_{ext}(t_{i+1})+h(1-\theta)  F_{ext}(t_{i}) \right] \\
   W &=&   \left[M + h\theta C + h^2 \theta^2 K \right]

* <b> Lagrangian Scleronomous Relations </b>

.. math::
     
   \dot y_{i+1}^{k+1} = G_0(Q^*(V_{i+1}^k))V_{i+1}^{k+1} \\
   P_{i+1}^{k+1} = G_0^t(Q^*(V_{i+1}^k))\lambda_{i+1}^{k+1} 

* <b> Lagrangian Rheonomous Relations </b>

.. math::
   
   \dot y_{i+1}^{k+1} &=& G_0(Q^*(V_{i+1}^k),t_{i+1})V_{i+1}^{k+1} +  G_1(Q^*(V_{i+1}^k, t_{i+1})) \\
   P_{i+1}^{k+1} &=& G_0^t(Q^*(V_{i+1}^k),t_{i+1})\lambda_{i+1}^{k+1} 

* <b> Lagrangian Compliant Relations </b>

.. math::

   \dot y_{i+1}^{k+1} &=& G_0(Q^*(V_{i+1}^k),\lambda_{i+1}^k)V_{i+1}^{k+1} +  G_1(Q^*(V_{i+1}^k, \lambda_{i+1}^k))\lambda_{i+1}^{k+1} \\
   P_{i+1}^{k+1} &=& G_0^t(Q^*(V_{i+1}^k),\lambda_{i+1}^k)\lambda_{i+1}^{k+1} 

* <b> Lagrangian Linear Relations </b>

.. math::

   \dot y_{i+1} &=& HV_{i+1} + D\lambda_{i+1} \\
P_{i+1} &=& H^t\lambda_{i+1}

The LCP formalization (DRAFT)
-----------------------------

Let us now consider the whole NSDS, indeed all the dynamical systems and all the interactions. The objective is to write a one-step non-smooth problem using the discretized equations above. Let us introduce the following variables and definitions:

:math:`\mathcal{I}_a = \{ \alpha / I_{\alpha} \ is \ active \}`, :math:`I_{\alpha}` being the interaction number :math:`\alpha`. The index :math:`a` and the sense of "active" depends on the type problem.
At the time, for first order systems in Siconos, all interactions are active and :math:`\mathcal{I}_a = \mathcal{I}_0` represents all the declared interactions.
For second order system, considering the newton-impact law for example, :math:`\mathcal{I}_a = \mathcal{I}_1 =\{\alpha / y^p_{\alpha} < 0 \}`.

Let :math:`\mathcal{DS}` denotes the set of all dynamical systems declared in the NSDS and :math:`\mathcal{DS}^{\alpha} \subset \mathcal{DS}` the set of systems involved in the interaction :math:`I_{\alpha}`. Finally, :math:`\mathcal{DS}^{\alpha,\beta}` represents the set of systems involved in :math:`I_{\alpha}` and :math:`I_{\beta}`.

Note that lower case letters are used to refer to dynamical systems and Greek letters to interactions.

First Order fully non-linear case
"""""""""""""""""""""""""""""""""

The dynamics write: 

.. math::

   (d)\left\lbrace\begin{array}{ccl}
   x_{i+1}^{k+1} &=& x^{free,k}_{i+1} + h(W_{i+1}^k)^{-1}r_{i+1}^{k+1} \\
   \\
   r_{i+1}^{k+1} &=& \displaystyle\sum_{\beta\in \mathcal{I}_d} {r_{i+1}^{k+1}}^\beta
   \end{array}\right.

The index "(d)" means these relations apply for dynamical system d. 

.. math::
   
   X^{\alpha}=\left[\begin{array}{c}
   x_d \\
   x_e \\
   \ldots
   \end{array}\right]&,&
   R^{\alpha}=
   \left[\begin{array}{c}
   r_d^{\alpha} \\
   r_e^{\alpha} \\
   \ldots
   \end{array}\right],
   with r_d = \sum_{\alpha} r_d^{\alpha}

Let us now consider the first order linear relation case, written for :math:`I_{\alpha}`:

.. math::
   
   y^{\alpha} &=& C^{\alpha}X^{\alpha} + ... \\
   R^{\alpha} &=& D^{\alpha}\lambda^{\alpha}

If :math:`\mathcal{DS}^{\alpha}` holds systems number :math:`d,e,\ldots`, then:

.. math::

   X^{\alpha}=\left[\begin{array}{c}
   x_d \\
   x_e \\
   \ldots
   \end{array}\right]&,&
   R^{\alpha}=
   \left[\begin{array}{c}
   r_d^{\alpha} \\
   r_e^{\alpha} \\
   \ldots
   \end{array}\right],
   with r_d = \sum_{\alpha} r_d^{\alpha}

Considering the discretized constraints and dynamics, one can easily reduce the system to the unknown :math:`(y,\lambda)`. Using the complementarity relation from the non-smooth law, this leads to a LCP formulation.
Before going more into details, let us introduce some new .

For example, considering linear system and linear relations, one gets:

.. math::
   Y_{i+1} = q_{LCP} + M_{LCP}\Lambda_{i+1}

:math:`Y_{i+1} (resp \Lambda_{i+1})` is a collection of all the :math:`y_{i+1}` vectors for all interactions in the non-smooth system. 

and block :math:`\alpha` (ie corresponding to interaction :math:`\alpha`) of :math:`q_{LCP}` is:

.. math::

   q_{LCP}^{\alpha} = \left[CX^{free}_{i+1} + e + FZ \right]_{\alpha}

:math:`M_{LCP}` is a block-matrix, with as diagonal blocks:

.. math::
   M_{LCP}^{\alpha,\alpha} = \left[ hC\hat W B + D \right]_{\alpha}

Lagrangian fully non-linear case
""""""""""""""""""""""""""""""""

The relations for a single interaction :math:`I_{\alpha}`, writes: 

.. math::

   \dot y_{i+1}^{k+1,\alpha} = G_0^{\alpha}(Q^*(V_{i+1}^k))V_{i+1}^{k+1,\alpha} \\
   p_{i+1}^{k+1,\alpha} = {G_0^{\alpha}}^t(Q^*(V_{i+1}^k))\lambda_{i+1}^{k+1,\alpha} 

with, for a single dynamical system, :math:`d`: 

.. math::

   v_{i+1}^{k+1,d} &=& v^{free,k,d}_{i+1} + (W_{i+1}^{k,d})^{-1} p_{i+1}^{k+1,d} \\
   p_{i+1}^{k+1,d} &=& \sum_{\beta \in \mathcal{I}_d} p_{i+1}^{k+1,\beta}

Then if we consider the whole NSDS, we get:

.. math::

   \dot Y_{i+1}^{k+1} = M_{lcp}\Lambda_{i+1}^{k+1} + q_{lcp}

:math:`\dot Y_{i+1}^{k+1}` (resp. :math:`\Lambda_{i+1}^{k+1}`) represents the vector that collects all :math:`\dot y_{i+1}^{k+1,\alpha} :math:`(resp. :math:`\lambda_{i+1}^{k+1,\alpha}`) for all interactions in the NSDS.
:math:`q_{lcp}` is a vector where each block :math:`\alpha` is given by:

.. math::
   q_{lcp}^{\alpha} = G_0^{\alpha}(\ldots)V^{free,k,alpha}_{i+1} + e\cdot \dot Y_i

and :math:`M_{lcp}` is a block matrix with for diagonal blocks:

.. math::
   M_{lcp}^{\alpha,\alpha} = G_0^{\alpha}(\ldots)(W_{i+1}^{k,\alpha})^{-1}{G_0^{\alpha}}^t(\ldots) 

and for extra-diagonal blocks:

.. math::
   M_{lcp}^{\alpha,\beta} = G_0^{\alpha}(\ldots)(W_{i+1}^{k,\alpha})^{-1}{G_0^{\beta}}^t(\ldots) 

TO BE REVIEWED

In Siconos, the class LagrangianDS represents non-linear Lagrangian  ndof-dimensional dynamical systems of the form :

.. math::
   \ddot q(t) &=& f_L(q,\dot q, t) + Tu + p \\
   q(t_0)&=&q_0 \\
   \dot q(t_0)&=&v_0 \\

Time discretisation of the dynamic

Considering a time step :math::`[t_i,t_{i+1}]` of length h, we denote :math:`a(t_i)=a_i`. 

Assuming:

.. math::

   \int_{[t_i,t_{i+1}]} M(q) \ddot q  \, dt \approx M(q*)(\dot q_{i+1}-\dot q{i}) \\
   p_{i+1} = \frac  1 h \int_{[t_i,t_{i+1}]} p \,d\nu

and using a \f$ \theta \f$-method to integrate the other terms, we obtain:

.. math::
   RES_{i+1}=M(q*)(\dot q_{i+1}-\dot q_{i})+h\left[(1-\theta)(Q_i+F_i)+\theta (Q_{i+1}+F_{i+1})-(1-\theta)F^{ext}_i-\theta F^{ext}_{i+1}\right]-hR_{i+1}=0  \\
   q_{i+1} = q_{i} +  h\left[\theta  \dot q_{i+1}+(1-\theta)  \dot q_{i}  \right]

\b Newton \b algorithm

For a known state i, we are looking for \f$ q_{i+1} \f$ solution of \f$ RES_{i+1}=0 \f$. \n
Applying Newton's method, we get: 

.. math::
   \dot q_{i+1}^{k+1}=\dot q_{i+1}^k-W^kRES_{i+1}^k \\
   with \\
   W^{-1}=\left[\frac{\partial{RES}}{\partial{\dot{q}}}\right]_{\dot q_{i+1}^k} 

Assuming: 

.. math::
   \left[\frac{\partial{M(q*)(\dot q_{i+1}^k-\dot q_{i})}}{\partial{\dot q }}\right]_{\dot q_{i+1}^k}&\approx&M(q*)\approx M(q_{i+1}^k)  \\
   \left[\frac{\partial{F^{ext}}}{\partial{\dot{q}}}\right]_{\dot q_{i+1}^k}&=& 0  

and denoting:

.. math::
   C_t^k=\left[\frac{\partial{(F+Q)}}{\partial{\dot{q}}}\right]_{\dot q_{i+1}^k} \\
   K_t^k=\left[\frac{\partial{(F+Q)}}{\partial{q}}\right]_{\dot q_{i+1}^k} 

we get:

.. math::
   W^{-1}=M(q_{i+1}^k)+h\theta C_t^k+h^2\theta^2K_t

and we write: 

.. math::
   \dot q_{i+1}^{k+1}=\dot q_{free}^k+W^khR_{i+1}^k 

with

.. math::
   
   \dot q_{free}^{k}&=&\dot q_{i+1}^k-W^kRES_{free}^k \\
   RES_{free}^k&=&M(q_{i+1}^k)(\dot q_{i+1}^k-\dot q_{i})+h\left[(1-\theta)(Q_i+F_i)+\theta (Q_{i+1}^k+F_{i+1}^k)-(1-\theta)F^{ext}_i-\theta {F^{ext}_{i+1}}^k\right] 

Relations discretization

The Time discretization of the relations is fully implicit and may be written as :

.. math::

   y_{i+1} &=& H^{T} q_{i+1} + b\\
   \dot y_{i+1} &=& H^{T}\dot q_{i+1} \\
   R_{i+1} &=& H \lambda_{i+1}

Time discretization of the Non Smooth laws

A natural way of discretizing the unilateral constraint  leads to the following implicit discretization :

.. math::

   0 \leq y_{i+1} \perp  \lambda_{i+1}  \geq 0

In the Moreau's time--stepping, we use a reformulation of the unilateral constraints in terms of velocity:

.. math::

   If \ y(t) =0, \ then \ 0 \leq \dot y \perp  \lambda  \geq 0

which leads to the following discretisation :

.. math::

   If \ y^{p} \leq 0, \ then \ 0 \leq \dot y_{i+1} \perp  \lambda_{i+1}  \geq 0

where \f$ y^{p} \f$ is a prediction of the position at time \f$ t_{i+1} \f$, for instance, \f$ {y^{p}}^k = y_{i+1}^k + \frac{h}{2}  \dot y_{i+1}^k \f$.

If we want to introduce now the Newton impact law, we consider an equivalent velocity defined by

.. math::

   {\dot y}^{e,k+1}_{i+1} = \dot y_{i+1}^{k+1} + e \dot y_{i+1}^k

finally, we apply the constraints directly on this velocity :

.. math::

   If \ y^{p} \leq 0, \ then \ 0 \leq {\dot y}^{e,k+1}_{i+1} \perp  \lambda_{i+1}^{k+1}  \geq 0

Summary of the time-discretized equations

.. math::
   
   \dot q_{i+1}^{k+1} &=& \dot q_{free}^k  + h W^k R_{i+1}^k \\
   q_{i+1}^{k+1} &=& q_{i} +  h\left[\theta  \dot q_{i+1}^{k+1}+(1-\theta)  \dot q_{i}  \right] \\
   \dot y_{i+1}^{k+1} &=& H^{T}\dot q_{i+1}^{k+1} \\
   R_{i+1}^{k} &=& H \lambda_{i+1}^{k}\\
   {y^{p}}^k &=& y_{i+1}^k + \frac{h}{2}  \dot y_{i+1}^k\\
   If \ &{y^{p}}^k& \leq 0, \ then \ 0 \leq {\dot y}^{e,k+1}_{i+1} \perp  \lambda_{i+1}^{k}  \geq 0

This set of equations can be reduced to a ``condensed'' system in terms of \f$ {\dot y}^{e,k+1}_{i+1} \f$ and \f$ {\lambda_{i+1}^k} \f$ :

.. math::

   {\dot y}^{e,k+1}_{i+1} &=&  H^{T} \dot q_{free}^k + h H^{T} W^k H \lambda_{i+1}^k  + e \dot y_{i+1}^k\\
   {y^{p}}^k &=& y_{i+1}^k + \frac{h}{2}  \dot y_{i+1}^k\\
   If \ &{y^{p}}^k& \leq 0, \ then \ 0 \leq {\dot y}^{e,k+1}_{i+1} \perp  \lambda_{i+1}^{k}  \geq 0

For LCP solver, we set

.. math::

   q^k_{lcp}&=&  H^{T} \dot q_{free}^k + e \dot y_{i+1}^k\\
   M^k_{lcp}&=&  h H^{T} W^k H

and we obtain: 

.. math::

   {\dot y}^{e,k+1}_{i+1} = M^k_{lcp}\lambda_{i+1}^k  + q^k_{lcp}

   
Lagrangian fully linear case
""""""""""""""""""""""""""""

.. math::

   \dot q_{i+1} &=& \dot q^{free}  + h W R_{i+1} \\
   q_{i+1} &=& q_{i} +  h\left[\theta  \dot q_{i+1}+(1-\theta)  \dot q_{i}  \right] \\
   \dot y_{i+1} &=& H^{T}\dot q_{i+1} \\
   R_{i+1} &=& H \lambda_{i+1}\\
   y^{p} &=& y_{i} + \frac{h}{2}  \dot y_i\\
   If &y^{p}& \leq 0, \ then  \ 0 \leq \dot y^{e}_{i+1} \perp  \lambda_{i+1}  \geq 0

This set of equations can be reduced to a "condensed" system in terms of :math:`\dot y^{e}_{i+1}` and :math:`{\lambda_{i+1}}` :

.. math::

   \dot y^{e}_{i+1} &=&  H^{T} \dot q^{free} + h H^{T} W H \lambda_{i+1}  + e \dot y_{i}\\
   y^{p} &=& y_{i} + \frac{h}{2}  \dot y_i\\
   If &y^{p}& \leq 0, \ then \ 0 \leq \dot y^{e}_{i+1} \perp  \lambda_{i+1}  \geq 0

The Simulation process
----------------------

As for Event-Driven, we introduce level index sets, with level = 0 for first order systems and level=1 for second order systems (this is related to the relative degrees but we won't get into details about that here).

:math:`I_0` is the set of all the potential UnitaryRelations (UR).
For second order systems:
:math:`I_1 = \{ ur_\alpha\in I_{0} , y^p_{\alpha} = 0 \}`. 
Thus, the LCP is built only for unitary relations that belongs to :math:`I_level`, level=0 for first order and level=1 for second order systems. 

Then, the steps of a Moreau's Time-Stepping simulation will be:

Knowing all values at the beginning of the time step :math:`[t_i,t_{i+1}]`,

-# compute the free solutions
-# for :math:`ur \in I_level` formalize and solve a LCP 
-# update the state (according to the possibly LCP results)
-# go to next time step

::
   
   SP::TimeStepping s(new TimeStepping(myModel));  
   SP::TimeDiscretisation t(new TimeDiscretisation(timeStep,s));

   s->initialize();

   int N = t->getNSteps(); // Number of time steps
   
   // --- Time loop ---
   while(k < N)// for each time step ...
   {
   // compute xFree, or qFree,vFree
   s->computeFreeStep();
   // Formalize and solve a LCP
   computeOneStepNSProblem("timeStepping");
   // Update state, using last computed values
   s->update(level); // 
   // transfer of state i+1 into state i and time incrementation
   s->nextStep();
   }

Note that all time-independent operators are computed during simulation initialisation.

Misc
----

This section is just a (draft) list of various unsorted topics related to Moreau Time Stepping.

- <b> How to set the behavior of the simulation according to the output value from the non-smooth solver call?</b>\n
Each time ComputeOneStepNS function, ie the Numerics solver, is called, it returns an int,
giving some information about the convergence of the solver:
output = 0 => solver succeeded. Else, the meaning of output depends on the solver called (see Numerics doc).\n
By default, when the convergence is not achieved, an exception is throwed and the process stops.
If you want to change that, it's possible by defining a specific function of the form::
  
  //
  // your inputFile.cpp
  //
  void myF(int info, SP::Simulation s)
  {
  // ...
  }
  
  int main(int argc, char* argv[])
  {
  // 
  // ...
  SP::TimeStepping s = ...
  s->setCheckSolverFunction(&myF);

Then after each call to s->computeOneStepNS(...), your function myF will be called.
That may be usefull to change the solver type, the tolerance or whatever you want. 

