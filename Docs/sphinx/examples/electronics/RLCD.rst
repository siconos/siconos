.. _RLCD_example:

Simulation of an electrical oscillator supplying a resistor through a half-wave rectifier
-----------------------------------------------------------------------------------------

Author: Pascal Denoyelle, September 22, 2005

You may refer to the source code of this example, `found here <https://github.com/siconos/siconos/blob/master/examples/Electronics/CircuitRLCD/CircuitRLCD.cpp>`_.

There is also a PDF version of this document, `viewable online here <https://github.com/siconos/siconos/blob/master/examples/Electronics/CircuitRLCD/Template-CircuitRLCD.pdf>`_.

Description of the physical problem : electrical oscillator with half-wave rectifier
````````````````````````````````````````````````````````````````````````````````````

In this sample, a LC oscillator initialized with a given voltage
across the capacitor and a null current through the inductor provides
the energy to a load resistance through a half-wave rectifier
consisting of an ideal diode (see :ref:`fig-CircuitRLCD`).

.. _fig-CircuitRLCD:

.. figure:: /figures/electronics/CircuitRLCD/SchemaCircuitRLCD.*
   :align: center

   fig 1: Electrical oscillator with half-wave rectifier

Only the positive wave of the oscillating voltage across the LC is
provided to the resistor. The energy is dissipated in the resistor
resulting in a damped oscillation.


Definition of a general abstract class of NSDS : the linear time invariant complementarity system (LCS)
```````````````````````````````````````````````````````````````````````````````````````````````````````

This type of non-smooth dynamical system consists of :

* a time invariant linear dynamical system (the oscillator). The state
  variable of this system is denoted by :math:`x`.

* a non-smooth law describing the behaviour of the diode as a
  complementarity condition between current and reverse voltage
  (variables (:math:`y,\lambda`) )

* a linear time invariant relation between the state variable
  :math:`x` and the non-smooth law variables (:math:`y,\lambda`)

Dynamical system and Boundary conditions
''''''''''''''''''''''''''''''''''''''''

Remark:

.. pull-quote::

   In a more general setting, the system's evolution would be
   described by a DAE :

   .. math::

      G \cdot x' = A \cdot x + E \cdot u + b + r

   with :math:`G , A , E` matrices constant over time (time invariant
   system), :math:`u , b` source terms functions of time and
   :math:`r`, a term coming from the non-smooth law variables :
   :math:`r = B \cdot \lambda + a` with :math:`B , a` constant over
   time.  We will consider here the case of an ordinary differential
   equation :

   .. math::

      x' = A \cdot x + E \cdot u + b + r

   and an initial value problem for which the boundary conditions are
   :math:`t_0 \in \mathbb{R} , x(t_0)= x_0`.

Relation between constrained variables and state variables
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

In the linear time invariant framework, the non-smooth law acts on the
linear dynamical system evolution through the variable :math:`r = B
\cdot \lambda + a`. Reciprocally, the state variable :math:`x` acts on the
non-smooth law through the relation :math:`y = C \cdot x + D \cdot
\lambda + F \cdot u + e` with :math:`C , D , F , e` constant over
time.

Definition of the Non Smooth Law between constrained variables
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

It is a complementarity condition between y and :math:`\lambda` :
:math:`0 \leq y \, \perp \, \lambda \geq 0`. This corresponds to the
behaviour of the rectifying diode, as described in
:ref:`non-smooth-laws`.

The formalization of the electrical oscillator with half-wave rectifier into the LCS
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The equations come from the following physical laws :

* the Kirchhoff current law (KCL) establishes that the sum of the
  currents arriving at a node is zero,
* the Kirchhoff voltage law (KVL) establishes that the sum of the
  voltage drops in a loop is zero,
* the branch constitutive equations define the relation between the
  current through a bipolar device and the voltage across it

Refering to :ref:`fig-CircuitRLCD`, the Kirchhoff laws could be written as:

.. math::

   \begin{array}{l}
   v_L = v_C\\
   v_R + v_D = v_C\\
   i_C + i_L + i_R = 0\\
   i_R = i_D
   \end{array}

while the branch constitutive equations for linear devices are:

.. math::

   \begin{array}{l}
   i_C = C v_C'\\
   v_L = L i_L'\\
   v_R = R i_R
   \end{array}

and last the "branch constitutive equation" of the ideal diode that is
no more an equation but instead a complementarity condition:

.. math::

   0 \leq i_D \, \perp \, -v_D \geq 0


This is illustrated in :ref:`fig-diode-reg`, where the left-hand
sketch displays the ideal diode characteristic and the right-hand
sketch displays the usual exponential characteristic as stated by
Shockley's law.

.. _fig-diode-reg:

.. figure:: /figures/electronics/CircuitRLCD/diode-caract.*
   :align: center
   :width: 13cm

   fig 2: Non-smooth and smooth characteristics of a diode


.. _sec-dyn-eq:

Dynamical equation
''''''''''''''''''

After rearranging the previous equations, we obtain:

.. math::

   \left( \begin{array}{c}
   v_L'\\
   i_L'
   \end{array} \right)
   =
   \left( \begin{array}{cc}
   0 & \frac{-1}{C}\\
   \frac{1}{L} & 0
   \end{array} \right)
   \cdot
   \left( \begin{array}{c}
   v_L\\
   i_L
   \end{array} \right)
   +
   \left( \begin{array}{c}
   \frac{-1}{C}\\
   0
   \end{array} \right)
   \cdot i_D

that fits in the frame of \ref{sec-def-NSDS} with

.. math::

   x =
   \left( \begin{array}{c}
   v_L\\
   i_L
   \end{array} \right)

and

.. math:: \lambda = i_D

Relations
'''''''''

We recall that the :math:`r = B \cdot \lambda + a` equation is expressed with

.. math::

   r =
   \left( \begin{array}{c}
   \frac{-1}{C}\\
   0
   \end{array} \right)
   \cdot i_D

from the dynamical equation (:ref:`sec-dyn-eq`).

Rearranging the initial set of equations yields:

.. math::

   -v_D =
   \left( \begin{array}{cc}
   -1 & 0
   \end{array} \right)
   \cdot
   \left( \begin{array}{c}
   v_L\\
   i_L
   \end{array} \right)
   + R i_D

as the second equation of the linear time invariant relation with

.. math:: y = -v_D



.. _non-smooth-laws:

Non Smooth laws
'''''''''''''''

There is just the complementarity condition resulting from the ideal diode characteristic:

.. math:: 0 \leq i_D \, \perp \, -v_D \geq 0

Description of the numerical simulation: the Moreau's time-stepping scheme
``````````````````````````````````````````````````````````````````````````

Time discretization of the dynamical system
'''''''''''''''''''''''''''''''''''''''''''

The integration of the ODE over a time step :math:`[t_i,t_{i+1}]` of length :math:`h` is:

.. math::

   \int_{t_i}^{t_{i+1}}x'\,dt = \int_{t_i}^{t_{i+1}}
   A \cdot x\,dt + \int_{t_i}^{t_{i+1}}(E \cdot u + b) dt
   + \int_{t_i}^{t_{i+1}}r\,dt

The left-hand term is :math:`x(t_{i+1})-x(t_i)`.

Right-hand terms are approximated this way:

* :math:`\int_{t_i}^{t_{i+1}} A \cdot x\,dt` is approximated using a
  :math:`\theta`-method

  .. math::

     \int_{t_i}^{t_{i+1}} A \cdot x\,dt \approx h \theta
     (A \cdot x(t_{i+1})) + h (1-\theta) (A \cdot x(t_{i}))

* since the second integral comes from independent sources, it can be
  evaluated with whatever quadrature method, for instance a
  :math:`\theta`-method

  .. math::

     \int_{t_i}^{t_{i+1}}(E \cdot u + b) dt \approx h \theta
     (E \cdot u(t_{i+1}) + b(t_{i+1})) +
     h (1-\theta) (E \cdot u(t_{i}) + b(t_{i}))

* the third integral is approximated like in an implicit Euler integration

  .. math:: \int_{t_i}^{t_{i+1}}r\,dt \approx h r(t_{i+1})


By replacing the accurate solution :math:`x(t_i)` by the approximated value :math:`x_i`, we get:

.. math::

   x_{i+1}-x_i = h \theta (A \cdot x_{i+1}) + h (1-\theta) (A \cdot x_{i}) +
                 h \theta (E \cdot u(t_{i+1}) + b(t_{i+1})) + h (1-\theta) (E \cdot u(t_{i}) + b(t_{i})) + h r_{i+1}

Assuming that :math:`I - h \theta A` is invertible, matrix :math:`W`
is defined as :math:`(I - h \theta A)^{-1}`. We get then:

.. math::

   x_{i+1} = W(I + h (1-\theta) A) \cdot x_{i} +
             W (h \theta (E \cdot u(t_{i+1}) + b(t_{i+1})) + h (1-\theta) (E \cdot u(t_{i}) + b(t_{i}))) + h W r_{i+1}

An intermediate variable :math:`x_{free}` related to the smooth part of the system is defined as:

.. math::

   x_{free} = W(I + h (1-\theta) A) \cdot x_{i} +
              W (h \theta (E \cdot u(t_{i+1}) + b(t_{i+1})) + h (1-\theta) (E \cdot u(t_{i}) + b(t_{i})))

Thus the calculus of :math:`x_{i+1}` becomes:

.. math:: x_{i+1} = x_{free} + h W r_{i+1}


Time discretization of the relations
''''''''''''''''''''''''''''''''''''

It comes straightforwardly:

.. math::

   r_{i+1} =& B \cdot \lambda_{i+1} + a

   y_{i+1} =& C \cdot x_{i+1} + D \cdot \lambda_{i+1} + F \cdot u(t_{i+1}) + e


Time discretization of the non-smooth law
'''''''''''''''''''''''''''''''''''''''''

It comes straightforwardly:

.. math:: 0 \leq y_{i+1} \, \perp \, \lambda_{i+1} \geq 0

Summary of the time discretized equations
'''''''''''''''''''''''''''''''''''''''''

These equations are summarized assuming that there is no source term
and simplified relations as for the electrical oscillator with
half-wave rectifier.

.. math::
   
   \begin{array}{ccc}
   W & = & (I - h \theta A)^{-1} \\
   x_{free} & = & W(I + h (1-\theta) A) \cdot x_{i} \\
   x_{i+1} & = & x_{free} + h W r_{i+1} \\
   r_{i+1} & = & B \cdot \lambda_{i+1}  \\
   y_{i+1} & = & C \cdot x_{i+1} + D \cdot \lambda_{i+1}  \\
   & 0 \leq y_{i+1} \, \perp \, \lambda_{i+1} \geq 0 &
   \end{array}

Numerical simulation
''''''''''''''''''''

The integration algorithm with a fixed step is described here :

   **Algorthm 1**: Integration of the electrical oscillator with half-wave
   rectifier through a fixed Moreau time stepping scheme

   **Require:** :math:`R > 0 , L > 0 , C > 0`

   **Require:** Time parameters :math:`h,T,t_0` and :math:`\theta` for the integration

   **Require:** Initial value of inductor voltage :math:`v_L = x_0(0)`

   **Require:** Optional, initial value  of inductor current :math:`i_L = x_0(1)` (default: 0)

   :math:`n_{step} = \frac{T - t_0}{h}`

   // *Dynamical system specification*

   :math:`A = \left( \begin{array}{cc} 0 & \frac{-1}{C}\\\frac{1}{L} & 0 \end{array} \right)`

   // *Relation specification*

   :math:`B = \left( \begin{array}{c}\frac{-1}{C}\\0\end{array} \right)`

   :math:`C = \left( \begin{array}{cc}-1 & 0\end{array} \right)`

   :math:`D = (R)`

   // *Construction of time independent operators*

   **Require:** :math:`I - h \theta A` invertible
   :math:`W = (I - h \theta A)^{-1}`
   :math:`M = D + h C W B`

   // *Non-smooth dynamical system integration*

   **for** :math:`i=0` to :math:`n_{step}-1` **do**

   .. math::

      x_{free} =& W (I + h (1 - \theta) A) x_i && \textrm{// Computation of $x_{free}$}

      q =& C \cdot x_{free} &&  \textrm{// Formalization of the one step LCP}

      (y_{i+1},\lambda_{i+1}) =& \textrm{solveLCP}(M,q) &&  \textrm{// One step LCP solving}

      x_{i+1} =& x_{free} + h W B \lambda_{i+1} &&  \textrm{// Computation of new state}

   **end for**


Comparison with numerical results coming from SPICE models and algorithms
`````````````````````````````````````````````````````````````````````````

We have used the SMASH simulator from Dolphin to perform a simulation
of this circuit with a smooth model of the diode as given by
Shockley's law , with a classical one step solver (Newton-Raphson) and
a choice between backward-Euler and trapezoidal integrators.

Characteristic of the diode in the SPICE model
''''''''''''''''''''''''''''''''''''''''''''''

:ref:`fig-carac-diode` depicts the static :math:`I(V)` characteristic
of two diodes with default SPICE parameters and two values for the
emission coefficient :math:`N`: 1.0 (standard diode) and 0.25 (stiff
diode).

.. _fig-carac-diode:

.. figure:: /figures/electronics/CircuitRLCD/caracdiode.*
   :align: center
   :width: 13cm

   fig 3: Diodes characteristics from SPICE model :math:`(N=0.25, N=1)`

The stiff diode is close to an ideal one with a threshold of 0.2 V.

Simulation results
''''''''''''''''''

:ref:`fig-comp-SMASH-SICONOS-BE10us` displays a comparison of the
SMASH and SICONOS results with a backward Euler integration and a
fixed time step of 10 μs. A stiff diode model was used in SMASH
simulations.  For :ref:`fig-comp-SMASH-SICONOS-TRAP10us`, a
trapezoidal integrator was used, yielding a better accuracy.  One can
notice that the results from both simulators are very close. The
slight differences are due to the smooth model of the diode used by
SMASH, and mainly to the threshold of around 0.2 V. Such a threshold
yields small differences in the conduction state of the diode with
respect to the ideal diode.

.. _fig-comp-SMASH-SICONOS-BE10us:

.. figure:: /figures/electronics/CircuitRLCD/comp_SMASH_SICONOS_BE10us.*
   :align: center
   :width: 13cm

   fig 4: SMASH and SICONOS simulation results with backward Euler
   integration, 10 μs time step

.. _fig-comp-SMASH-SICONOS-TRAP10us:

.. figure:: /figures/electronics/CircuitRLCD/comp_SMASH_SICONOS_TRAP10us.*
   :align: center
   :width: 13cm

   fig 5: SMASH and SICONOS simulation results with trapezoidal
   integration, 10 μs time step
