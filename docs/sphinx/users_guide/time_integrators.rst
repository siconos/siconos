.. _time_integrators:

Time integration of the dynamics
================================

Dynamical systems integration over a time-step or between two events must be defined thanks to
'one-step integrators'.

* Euler-Moreau (:class:`EulerMoreauOSI`)

For first-order dynamical systems, in an 'event-capturing' simulation strategy.
  
.. math::

   M x_{k+1} &= M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1}) + h(1-\gamma)r(t_k) \\
   y_{k+1} &=  h(t_{k+1},x_{k+1},\lambda _{k+1}) \\
   r_{k+1} &= g(x_{k+1},\lambda_{k+1},t_{k+1})\\

with a nonsmooth law linking :math:`y_{k+1}` and :math:`\lambda_{k+1}`,
and :math:`\theta \in [0,1], \gamma \in [0,1]`.

Another variant can also be used (FullThetaGamma scheme)

.. math::

   M x_{k+1} = M x_{k} +h f(x_{k+\theta},t_{k+1}) + h r(t_{k+\gamma}) \\[2mm]
   y_{k+\gamma} =  h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
   r_{k+\gamma} = g(x_{k+\gamma},\lambda_{k+\gamma},t_{k+\gamma})\\[2mm]
   \mbox{nslaw} ( y_{k+\gamma} , \lambda_{k+\gamma})

* Moreau-Jean (:class:`MoreauJeanOSI`)

For mechanical (second-order) systems, in an 'event-capturing' simulation strategy.

.. math::
   
   M (v_{k+1}-v_k) + h K q_{k+\theta} + h C v_{k+\theta} - h F_{k+\theta} = p_{k+1} = G P_{k+1}\\ 
   q_{k+1} = q_{k} + h v_{k+\theta}, \\
   U_{k+1} = G^\top\, v_{k+1}, \\
   \begin{array}{lcl}
   0 \leq U^\alpha_{k+1} + e  U^\alpha_{k} \perp P^\alpha_{k+1}  \geq 0,& \quad&\alpha \in \mathcal I_1, \\
   P^\alpha_{k+1}  =0,&\quad& \alpha \in \mathcal I \setminus \mathcal I_1,\end{array}

with  :math:`\theta \in [0,1]`. The index set :math:`\mathcal I_1` is the discrete equivalent
to the rule that allows us to apply the Signorini  condition at the velocity level.
Numerically, this set is defined as

.. math::

   \mathcal I_1 = \{\alpha \in \mathcal I \mid G^\top (q_{k} + h v_{k}) + w \leq 0\text{ and } U_k \leq 0 \}.

* Schatzman-Paoli (:class:`SchatzmanPaoliOSI`)
* zero-order  (:class:`ZeroOrderHoldOSI`) 
* Lsodar (:class:`LsodarOSI`)

For 'event-driven' simulation strategy. Integrator based on LSODAR (https://computation.llnl.gov/casc/odepack/) rootfinding routine :
"Lsodar solves problems dy/dt = f with full or banded Jacobian and automatic method selection, and at the same time, it finds the roots of any of a set of given functions of the form g(t,y). This is often useful for finding stop conditions or points at which switches are to be made in the function f". 
In Siconos, Lsodar is used for event-driven algorithm, to integrate the dynamics with stops at new non-smooth events (violation of a constraint)

* Hem5 (:class:`Hem5OSI`)

For 'event-driven' simulation strategy. Based on Ernst Hairer HEM5 integrator (http://www.unige.ch/~hairer/software.html)

* Newmark (:class:`NewMarkAlphaOSI`)
