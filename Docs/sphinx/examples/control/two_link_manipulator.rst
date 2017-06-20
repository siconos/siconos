.. _two_link_manipulator:

Two-link manipulator
====================

Authors: Bernard Brogliato, Irinel-Constantin Morarescu.

Tracking control of the trajectory described by the end point of a two-link planar manipulator

*Keywords*: :doxysiconos:`LagrangianLinearTIDS`, :doxysiconos:`LagrangianScleronomousR`, :doxysiconos:`NewtonImpactNSL`, :doxysiconos:`TimeStepping`, :doxysiconos:`EulerMoreauOSI`, :doxysiconos:`LCP`

Files: Control/Two-linkManipulator/Two-linkManipulator.cpp, Control/Two-linkMultiConstraints/Two-linkMulticonstrManip.cpp

Two-link planar manipulator, general framework
----------------------------------------------

The model belongs to the class of complementary Lagrangian systems subject to frictionless unilateral constraints whose dynamics may be expressed as:

.. math::

   \left\{ \begin{array}{ccc}
   M(X)\ddot{X}_{1}+C(X,\dot{X})\dot{X}+G(X) = U +\nabla F(X)\lambda_{X}\\
   0\leq\lambda_X\perp F(X)\geq0,\hspace{2.7cm}\\ 
   \mbox{Collision rule} \hspace{5.5cm} & \end{array} \right. 

where

* :math:`X\in\mathcal{R}^{n}` is the vector of generalized coordinates,
* :math:`M(X)=M^{T}(X)\in\mathcal{R}^{n\times n}` is the positive definite inertia matrix,
* :math:`F(X)\in\mathcal{R}^{m}` represents the distance to the constraints, 
* :math:`C(X,\dot{X})` is the the matrix containing Coriolis and centripetal forces,
* :math:`G(X)` contains conservative forces,
*  :math:`\lambda_{X}\in\mathcal{R}^{m}` are the Lagrangian multipliers associated to the constraints,
* :math:`U\in\mathcal{R}^{n}` is the vector of generalized torque inputs.

For the sake of completeness we precise that :math:`\nabla` denotes the Euclidean gradient :math:`\nabla F(X)=(\nabla F_{1}(X), \ldots,F_{m}(X))\mathcal{R}^{n\times m}` where :math:`\nabla F_{i}(X)\in\mathcal{R}^{n}` represents the vector of partial derivatives of :math:`F_{i}` with respect to the components of :math:`X`. We assume that the functions :math:`F_{i}` are continuously differentiable and that :math:`\nabla F_{i}(X(t_{\ell}))\neq0` for all impact times :math:`t_{\ell}`. It is worth to precise here that for a given function :math:`f` its derivative with respect to the time :math:`t` will be denoted by :math:`\dot{f}`.

Dynamics Equation based on Lagrangian formulation
-------------------------------------------------

Let us introduce the following notations:

* :math:`\theta_{i}` - the joint angle of the joint :math:`i`,
* :math:`m_{i}` - the mass of link :math:`i`,
* :math:`I_{i}` - the moment of inertia of link :math:`i` about the axis that passes through the center of mass and is parallel to the :math:`Z` axis,
* :math:`l_{i}` - the length of link :math:`i`,
* :math:`g` - the gravitational acceleration.

.. image:: /figures/control/Two-linkManipulator/planar_manipulator.*

In order to simplify the presentation we assume that the center of mass of each link is right at the middle of the link. We will also assume that the gravitational force acts in the negative direction of the :math:`Y` axis.

Let us choose :math:`\theta=(\theta_{1},\theta_{2})^{T}` as the generalized coordinates and denote :math:`T` the kinetic energy, :math:`P` the potential energy and :math:`L=T-P` the Lagrangian function. Using the notations above and taking into account the assumption that we have made, the potential and the kinetic energy of the system can be expressed as follows: 

.. math::
   
   P&=\frac{1}{2} m_{1}gl_{1}\sin{\theta_{1}}+m_{2}g\left(l_{1}\sin\theta_{1}+\frac{l_{2}}{2}\sin(\theta_{1}+\theta_{2})\right),\\ 
   T&=\frac{ m_{1}l_{1}^{2}\dot{\theta_{1}}^{2}}{8}+\frac{m_{2}}{2}\left(l_{1}^{2}\dot{\theta}_{1}^{2}+\frac{l_{2}^{2}}{4}(\dot{\theta}_{1}+\dot{\theta}_{2})^{2}+l_{1}l_{2}(\dot{\theta}_{1}^{2}+\dot{\theta}_{1}\dot{\theta}_{2})\cos\theta_{2}\right)+\frac{I_{1}\dot{\theta_{1}}^{2}+I_{2}(\dot{\theta_{1}}+\dot{\theta_{2}})^{2}}{2}. 

Substituting :math:`P` and :math:`T` in the formula of :math:`L` the dynamics equation of the two-links manipulator: 

.. math::
   \frac{\mbox{d}}{\mbox{d}t}\left(\frac{\partial L}{\partial\dot{\theta}}\right)-\frac{\partial L}{\partial\theta}=0 

can be rewritten as 

.. math::

   \label{dynamics of two-links manipulator} 
   M(\theta)\ddot{\theta}+C(\theta,\dot{\theta})\dot{\theta}+G(\theta)=0 

where :math:`M,C\in \cal{R}^{2\times 2},G\in R^{2}` 
have the same meaning as in (general dynamics) and they are explicitly given by: 

.. math::

   \begin{array}{ccc}
   & M=\left[\begin{array}{cc}
   M_{11} & M_{12}\\ 
   M_{21} & M_{22}\\ 
   \end{array}\right],\quad & 
   \left\{ \begin{array}{ccc}
   M_{11} & = & \displaystyle\frac{m_{1}l_{1}^{2}}{4}+m_{2}\left(l_{1}^{2}+\frac{l_{2}^{2}}{4}l_{1}l_{2}\cos\theta_{2}\right)+I_{1}+I_{2}\hspace{0.5cm}\\
   M_{12} & = & M_{21} =\displaystyle\frac{m_{2}l_{2}^{2}}{4}+\frac{m_{2}l_{1}l_{2}}{2}\cos{\theta_{2}}+I_{2}\hspace{2cm}\\
   M_{22} & = & \displaystyle\frac{m_{2}l_{2}^{2}}{4}+I_{2}\hspace{5.8cm}\\\end{array}\right.,\\ &C=\left[\begin{array}{cc} C_{11} & C_{12}\\ C_{21} & C_{22}\\ \end{array}\right],\quad &\left\{\begin{array}{ccc} C_{11} & = & -m_{2}l_{1}l_{2}\dot{\theta}_{2}\sin\theta_{2}\\\
   C_{12} & = & -\displaystyle\frac{m_{2}l_{1}l_{2}}{2}\dot{\theta}_{2}\sin\theta_{2}\\
   C_{21} & = & \displaystyle\frac{m_{2}l_{1}l_{2}}{2}\dot{\theta}_{1}\sin\theta_{2}\\
   C_{22} & = & 0\hspace{2.5cm}\\ \end{array}\right.,\\ &G=\left[\begin{array}{c} G_{1} \\ G_{1} \\ \end{array}\right],\hspace{1cm}\quad &\left\{\begin{array}{ccc} G_{1} & = & \displaystyle\frac{g}{2}[l_{1}(2m_{1}+m_{2})\cos\theta_{1}+m_{2}l_{2}\cos(\theta_{1}+\theta_{2})]\\
   G_{2} & = & \displaystyle\frac{m_{2}gl_{2}}{2}\cos(\theta_{1}+\theta_{2})\hspace{4.5cm}\\ \end{array}\right..
  \end{array}

Constrained problem formulation
-------------------------------

General Consideration
"""""""""""""""""""""

We are interested on the problem of control of the trajectory described by the end point of the manipulator's second link. The constraint surface corresponds to the ground (i.e. :math:`y=0`). Obviously the associated admissible domain is :math:`\Phi=\{(x,y)\mid y\geq0\}`. In order to apply the previous theoretical considerations we must consider a coordinates transformation. Entering into details, if :math:`(x,y)` are the Cartesian coordinates of this point, we will consider the generalized coordinates

.. math::

   q=\left[\begin{array}{c} 
   y \\ 
   x \\ 
   \end{array}\right], y\geq0 

The coordinates transformation is simply given by 

.. math::
   
   \left\{\begin{array}{c} y=l_{1}\sin\theta_{1}+l_{2}\sin(\theta_{1}+\theta_{2})\\
   x=l_{1}\cos\theta_{1}+l_{2}\cos(\theta_{1}+\theta_{2})\\ \end{array}\right. 

and the corresponding Jacobian matrix can be easily derived as 

.. math::

   \label{jacobian matrix} 
   J_{q}(\theta)=\left[\begin{array}{cc} l_{1}\cos\theta_{1}+l_{2}\cos(\theta_{1}+\theta_{2}) & l_{2}\cos(\theta_{1}+\theta_{2})\\ -l_{1}\sin\theta_{1}-l_{2}\cos(\theta_{1}+\theta_{2}) & -l_{2}\sin(\theta_{1}+\theta_{2})\\ \end{array} \right].
   
We call a singular configurations for the system above as those for which the end-effector velocities in a certain direction can not be realized by any finite joint velocity.

Without entering into details, from the mathematical point of view, singular configurations can be characterized using the Jacobian matrix. In the case of the two-link manipulator the singular configurations are given by: :math:`\det J_{q}(\theta)=0\Leftrightarrow l_{1}l_{2}\sin\theta_{2}=0\Leftrightarrow \theta_{2}=0^{\circ}` or :math:`\theta_{2}=180^{\circ}` (see figure below\}, the configurations A and B do not allow a velocity in the direction of the origin (or the opposite) realized by finite joint velocities).

.. image:: /figures/control/Two-linkManipulator/Singular.*
	   
We consider only one unilateral constraint for the system associated to the end point of the manipulator's second link. Therefore, we do not take into account the case when some other parts of the manipulator touch the ground. Let us consider that the system must accomplish a cyclic task consisting of tracking a circle that violates the constraint. In order to track the trajectory the manipulator must follow the ground line from the point where the circle leave the admissible domain to the point where the circle re-enter in it. Thus, there exists a free-motion phase and a constraint motion phase during which a contact force is imposed.

Controller design
"""""""""""""""""

The controller used here consists of different low-level control laws for each phase of the system. More precisely, the controller can be expressed as

.. math::

   T(q)U=\left\{ \begin{array}{cc} U_{nc} & \mbox{for } t\in\Omega_{2k}\\
   U_{t} & \mbox{for } t\in I_{k}\\
   U_{c} & \mbox{for } t\in\Omega_{2k+1} \end{array}\right., 

where :math:`T(q)=\left(\begin{array}{c} T_{1}(q)\\ T_{2}(q) \end{array}\right)\in\mathcal{R}^{n\times n}`. The new coordinates :math:`q` are chosen such that :math:`\Phi=\{q\mid Dq\geq0\}`. The tangent cone :math:`T_{\Phi}(q_{1}=0)=\{v\mid Dv\geq0\}` is the space of admissible velocities on the boundary of :math:`\Phi`.

The controller law used in the following is based on the fixed-parameter scheme presented by J.J. Slotine. First, let us introduce some notations: :math:`\tilde{q}=q-q_{d},\,\bar{q}=q-q_{d}^{*},\,s=\dot{\tilde{q}}+\gamma_{2}\tilde{q}, \,\bar{s}=\dot{\bar{q}}+\gamma_{2}\bar{q},\,q_{r}=\dot{q}_{d}-\gamma_{2}\tilde{q}` where :math:`\gamma_{2}>0` is a scalar gain and :math:`q_{d},\,q_{d}^{*}` will be explicitly defined in the next section. Using the notations above the controller is given by

.. math::

   \label{Slotine scheme}
   T(q)U=\left\{ \begin{array}{ccc} 
   U_{nc} & =  M(q)\ddot{q}_{r}+C(q,\dot{q})\dot{q}_{r}+g(q)-\gamma_{1}s \hspace{3.4cm}\\ 
   U_{t} & =  U_{nc} \mbox{ before the first impact}\hspace{4.2cm}\\ U_{t} & = & M(q)\ddot{q}_{r}+C(q,\dot{q})\dot{q}_{r}+g(q)-\gamma_{1}\bar{s} \mbox{ after the first impact}\\ U_{c} & =  U_{nc}-P_{d}+K_{f}(P_{q}-P_{d})\hspace{4.4cm} \end{array}\right.

where :math:`\gamma_{1}>0` is a scalar gain, :math:`K_{f}>0,\,P_{q}=D^{T}\lambda` and :math:`P_{d}=D^{T}\lambda_{d}` is the desired contact force during constraint motion.

Desired trajectory
""""""""""""""""""

First of all we split the time axis into intervals :math:`\Omega_{k}` and :math:`I_{k}` corresponding to specific phases of the system. Precisely, :math:`\Omega_{2k}` corresponds to free-motion phases (:math:`F(X)>0`) and :math:`\Omega_{2k+1}` corresponds to constrained-motion phases (:math:`F_{i}(X)=0` for some index :math:`i\in\{1,\ldots,m\}`). Therefore, during the :math:`\Omega_{k}` phases no impact can occur. Between a free and a constrained phase the dynamical systems always passes into a transition phase :math:`I_{k}` containing some impacts. Since the dynamics of the system does not change during the transition between a constrained and a free-motion phase, in time domain one gets the following typical task representation): 

.. math::

   \label{task}
   \mathcal{R}^{+}=\Omega_{0}\cup I_{0}\cup\Omega_{1}\cup\Omega_{2}\cup I_{1}\cup\ldots\cup\Omega_{2k}\cup I_{k}\cup\Omega_{2k+1}\cup\ldots 

The sequence :math:`\Omega_{2k}\cup I_{k}\cup\Omega_{2k+1}` will be referred as the cycle :math:`k` of the system evolution.Consider the following notations:

* :math:`t_{0}^{k}` is the first impact during the cycle :math:`k`,
* :math:`t^{*k}` is the time corresponding to :math:`q^{*}_{1d}(t^{*k})=0`,
* :math:`t_{\infty}^{k}` is the accumulation point of the sequence :math:`\{t_{\ell}^{k}\}_{\ell\geq0}` of the impact instants during the cycle :math:`k` (obviously :math:`t_f^{k}\geq t_{\infty}^{k}`), 
* :math:`\tau_{1}^{k}` is such that :math:`q^{*}_{1d}(\tau_{1}^{k})=-\varphi V_{1}(\tau_{0}^{k})` and :math:`\dot{q}^{*}_{1d}(\tau_{1}^{k})=0`, where :math:`\varphi>0` is chosen by the designer in order to impose a closed-loop dynamics with impacts,  
* :math:`t_{d}^{k}` is the detachment instant, therefore :math:`\Omega_{2k+1}=[t_{f}^{k},t_{d}^{k}]`.

It is noteworthy that :math:`t_0^k,\,t_\infty^k,t\,_d^k` are state dependent whereas :math:`t^{*k},\,\tau_1^k` and :math:`\tau_0^k` are exogenous and imposed by the designer. On :math:`[\tau_{0}^{k},t_{0}^{k})` we impose that :math:`q^{*}_{d}(\cdot)` is twice differentiable and :math:`q^{*}_{1d}(t)` decreases towards :math:`-\varphi V_{1}(\tau_{0}^{k})` on :math:`[\tau_{0}^{k},\tau_{1}^{k}]`. For the sake of simplicity, in order to satisfy the previous requirements we define on :math:`[\tau_{0}^{k},\tau_{1}^{k}]` the signal :math:`q^{*}_{1d}(\cdot)` as a degree 3 polynomial function with limit conditions (:math:`t_{ini}=\tau_{0}^{k}` and :math:`t_{end}=\tau_{1}^{k}`). Therefore, 

.. math::

   \label{desired trajectory} 
   \begin{array}{ccc} q^{*}_{1d}&= a_{3}t^{3}+a_{2}t^{2}+a_{1}t+a_{0}\\ \dot{q}^{*}_{1d}&=3a_{3}t^{2}+2a_{2}t+a_{1}\hspace{1cm} \end{array} 

with the coefficients given by: 

.. math::

   \label{desired trajectory coefficients} 
   \begin{array}{ccc} a_{3}&=& 2[q_{1d}(\tau_{0}^{k})+\varphi V_{1}(\tau_{0}^{k})]\hspace{0.3cm}\\ a_{2}&=&-3[q_{1d}(\tau_{0}^{k})+\varphi V_{1}(\tau_{0}^{k})]\\ a_{1}&=&0\hspace{4cm}\\ a_{0}&=&q_{1d}(\tau_{0}^{k})\hspace{3cm} \end{array}.

The signal :math:`q^{*}_{2d}(t)\in C^{2}(\mathcal{R}^{+})` is frozen during the transition phase:

* :math:`q^{*}_{2d}(t)=q^{*}_{2d},\,\dot{q}^{*}_{2d}(t)=0` on :math:`[\tau_{0}^{k},t_{\infty}^{k}]`,
* :math:`q^{*}_{2d}(t)` is defined such that :math:`\dot{q}^{*}_{2d}(t^{*k})=0`.

On :math:`(t_{0}^{k},t_{f}^{k}]` we set :math:`q_{d}` and :math:`q^{*}_{d}` as follows: 

.. math::

   \label{desired trajectory definition on transition phases} 
   q_{d}=\left(\begin{array}{c} 0\\ 
   q^{*}_{2d} \end{array}\right),\quad q^{*}_{d}=\left(\begin{array}{c} -\varphi V_{1}(\tau_{0}^{k})\\ 
   q^{*}_{2d} \end{array}\right), 

and on :math:`[t_{f}^{k},t_{d}^{k}]` we set 

.. math::

   \label{desired trajectory definition on constrained phases} 
   q_{d}=\left(\begin{array}{c} 0\\ 
   q_{2d}(t) \end{array}\right),\quad q^{*}_{1d}=0. 

We note that :math:`q_{d}=q^{*}_{d}` on :math:`(t_{f}^{k},t_{d}^{k})`.

The Formalization of the problem into the class of Lagrangian NSDS
------------------------------------------------------------------

Second order non linear Lagrangian dynamical systems
""""""""""""""""""""""""""""""""""""""""""""""""""""

From the input of the physical data, we construct all of the terms which defined a Lagrangian NSDS. In our special case:

* :math:`M` is given in the previous sections,
* :math:`fGyr=C(q,\dot{q})\dot{q}+G(q)`,
* :math:`F_{int}` is identically zero,
* :math:`F_{ext}` is used to introduce the control law :math:`U`.

Relations
"""""""""

* The unilateral constraint requires that :

.. math::

   y \geq 0 

* Physical consideration impose :math:`0\leq\theta_1\leq\pi`
* In order to avoid the singular cases we impose :math:`-\pi<\theta_2<0`

Non Smooth laws
"""""""""""""""

There exists just one unilateral constraint such that : 

.. math::

   0 \leq y \perp \lambda\geq 0 

The Newton impact law at impact is given by : 

.. math::

   if \ y=0,\quad \dot y(t^+)= -e \dot y(t^-)

Exploitation of the results
---------------------------

We present here just some basic results concerning the tracking of the trajectory. More precisely, we present the variation of :math:`y` in time and the path of the end point of the manipulator's second link in :math:`(x,y)`-plane. The variation of other quantities may be also obtained by the user.

.. image:: /figures/control/Two-linkManipulator/two-linkManipulatorResults1.*

.. image:: /figures/control/Two-linkManipulator/two-linkManipulatorResults2.*
