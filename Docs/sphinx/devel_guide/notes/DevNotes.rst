OneStepNSProblem formalisation for several interactions
=======================================================

+-----------+----------------+
| author    | F. Pérignon    |
+===========+================+
| date      | May 16, 2006   |
+-----------+----------------+
| version   | ?              |
+-----------+----------------+

LinearDS - Linear Time Invariant Relations
------------------------------------------

General notations
~~~~~~~~~~~~~~~~~

We consider :math:`n` dynamical systems of the form:

.. math:: \dot x_i = A_i x_i + R_i

| Each system if of dimension :math:`n_i`, and we denote
:math:`N = \displaystyle{\sum_{i=1}^{n} n_i}`.
| An interaction, :math:`I_{\alpha}` is composed with a non smooth law,
:math:`nslaw_{\alpha}` and a relation:

.. math:: y_{\alpha} = C_{\alpha}X_{\alpha} + D_{\alpha}\lambda_{\alpha}

| The “dimension” of the interaction, ie the size of vector
:math:`y_{\alpha}`, is denoted :math:`m_{\alpha}` and we set:

.. math:: M = \sum_{\alpha=1}^{m} m_{\alpha}

:math:`m` being the number of interactions in the Non Smooth Dynamical
System.
| :math:`X_{\alpha}` is a vector that represents the DS concerned by the
interaction. Its dimension is noted :math:`N_{\alpha}`, this for
:math:`n_{\alpha}` systems in the interaction.
| :math:`C_{\alpha}` is a :math:`m_{\alpha} \times N_{\alpha}`
row-blocks matrix and :math:`D_{\alpha}` a
:math:`m_{\alpha} \times m_{\alpha}` square matrix.

.. math::

   C_{\alpha}=\left[\begin{array}{ccc} 
   C_{\alpha}^i & C_{\alpha}^j & ...\end{array}\right]

| with :math:`i,j,...\in \mathcal{DS}_{\alpha}` which is the set of DS
belonging to interaction :math:`\alpha`.
| We also have the following relation:

.. math::

   \left[\begin{array}{c} 
   R_{\alpha}^i \\
   R_{\alpha}^j \\
   ...  
   \end{array}\right] = B_{\alpha}\lambda_{\alpha}
   =\left[\begin{array}{c} 
   B_{\alpha}^i \\
   B_{\alpha}^j \\
   ...
   \end{array}\right]\lambda_{\alpha}

| :math:`R_{\alpha}^i` represents the contribution of interaction
:math:`\alpha` on the reaction of the dynamical system :math:`i`, and
:math:`B_{\alpha}^i` is a :math:`n_i \times m_{\alpha}` block matrix.
| And so:

.. math:: R_i = \sum_{\beta\in\mathcal{I}_i}R_{\beta}^i=\sum_{\beta\in\mathcal{I}_i}B^i_{\beta} \lambda_{\beta}

| with :math:`\mathcal{I}_i` the set of interactions in which dynamical
system number :math:`i` is involved.
| Introducing the time discretization, we get:

.. math::

   \begin{aligned}
   x_i^{k+1}-x_i^k = h A_i x_i^{k+1} + h R_i^{k+1}  \\
   \nonumber\\
   y_{\alpha}^{k+1} = C_{\alpha}X_{\alpha}^{k+1} + D_{\alpha}\lambda_{\alpha}^{k+1}\\
   \nonumber\\
   R_i^{k+1} = \sum_{\beta\in\mathcal{I}_i}B^i_{\beta} \lambda_{\beta}^{k+1}\end{aligned}

ie, with :math:`W_i = (I-h A_i)^{-1}`:

.. math::

   \begin{aligned}
   x_i^{k+1}&=& W_i x_i^{k} + hW_i R_i^{k+1}  \\
   \nonumber\\
   y_{\alpha}^{k+1} &=& C_{\alpha}W_{\alpha} X_{\alpha}^{k} + C_{\alpha}hW_{\alpha}\sum_{\beta\in\mathcal{I}_i}B^i_{\beta} \lambda_{\beta}^{k+1} + D_{\alpha}\lambda_{\alpha}^{k+1} \\
   &=& C_{\alpha}W_{\alpha} X_{\alpha}^{k} + (C_{\alpha}hW_{\alpha}B_{\alpha} + D_{\alpha}) \lambda_{\alpha}^{k+1} + \sum_{\beta\neq\alpha}(\sum_{i\in\mathcal{DS}_{\alpha}\cap\in\mathcal{DS}_{\beta}} hC_{\alpha}^iW_i B^i_{\beta} \lambda_{\beta}^{k+1})\end{aligned}

with

.. math::

   \label{Walpha}
   W_{\alpha}=\left[\begin{array}{ccc} 
   W_i &  0   & ... \\
   0   &  W_j & ...\\
   0  & ... & ... \\ 
   \end{array}\right]

| the block-diagonal matrix of all the :math:`W` for the dynamical
systems involved in interaction :math:`\alpha`.
| The global-assembled :math:`Y` vector, of dimension M, composed by
:math:`m` :math:`y_{\alpha}` subvectors, is given by:

.. math::

   \begin{aligned}
   Y_{k+1} = q_{OSNSP} + M_{OSNSP}\Lambda_{k+1}\end{aligned}

or,

.. math::

   \begin{aligned}
   Y_{k+1} =\left[\begin{array}{c} 
   y_1 \\
   ...  \\
   y_m
   \end{array}\right]_{k+1}
   &=&\left[\begin{array}{ccc} 
   C_1^1 & \ldots & C_1^n \\
   \vdots & \ldots & \vdots \\
   C_m^1 & \ldots & C_m^n 
   \end{array}\right]\left[\begin{array}{cccc} 
   W_1 & 0 & \ldots &0 \\
   0  & W_2 & \ddots & \vdots \\
   \vdots &\ddots  & \ddots & \vdots \\
   &&0& W_n
   \end{array}\right]
   \left[\begin{array}{c} 
   x_1  \\
   \vdots \\
   \vdots \\
   x_n 
   \end{array}\right]_k \\
   &+&\left[\begin{array}{cccc} 
   D_1+h\sum_{j\in \mathcal{DS}_1}C_1^jW_jB_1^j & h\displaystyle{\sum_{j\in \mathcal{DS}_1\cap\mathcal{DS}_2}C_1^jW_jB_2^j} & \ldots &\\
   \vdots&\ddots& &\\
   & h\displaystyle{\sum_{j\in \mathcal{DS}_m}C_m^jW_jB_{m-1}^j}  & D_m+h\displaystyle{\sum_{j\in \mathcal{DS}_m\cap\mathcal{DS}_{m-1}}C_m^jW_jB_m^j} \\
   \end{array}\right]\left[\begin{array}{c} 
   \lambda_1  \\
   \vdots \\
   \lambda_m 
   \end{array}\right]_{k+1} \nonumber\end{aligned}

To sum it up, the block-diagonal term of matrix :math:`M_{OSNSP}`, for
block-row :math:`\alpha` is:

.. math:: D_{\alpha}+h\sum_{j\in \mathcal{DS}_{\alpha}}C_{\alpha}^jW_jB_{\alpha}^j

This is an :math:`m_{\alpha}\times m_{\alpha}` square matrix. The
extra-diagonal block term, in position (:math:`\alpha,\beta`) is:

.. math:: h\sum_{j\in \mathcal{DS}_{\alpha}\cap\mathcal{DS}_{\beta}}C_{\alpha}^jW_jB_{\beta}^j

| and is a :math:`m_{\alpha}\times m_{\beta}` matrix. This matrix
differs from 0 when interactions :math:`\alpha` and :math:`\beta` are
coupled, ie have common DS.

Or, for the relation l of interaction :math:`\alpha`, we get:

.. math:: D_{\alpha,l}+h\sum_{j\in \mathcal{DS}_{\alpha}}C_{\alpha,l}^jW_jB_{\alpha}^j

for the diagonal, and

.. math:: h\sum_{j\in \mathcal{DS}_{\alpha}\cap\mathcal{DS}_{\beta}}C_{\alpha,l}^jW_jB_{\beta}^j

| for extra-diagonal terms.
| :math:`D_{\alpha,l}`, row number :math:`l` of :math:`D_{\alpha}`, the
same for :math:`C_{\alpha,l}`

Finally, the linked-Interaction map provides, for each interaction
(named “current interaction”), the list of all the interactions (named
“linked interaction”) that have common dynamical system with the
“current interaction”.

A simple example
~~~~~~~~~~~~~~~~

We consider :math:`n=3` dynamical systems and :math:`m=2` interactions:

.. math::

   \begin{aligned}
   I_{\mu}& \rightarrow& \mathcal{DS}_{\mu} = \{DS_1, DS_3\}, m_{\mu} = 3 \\
   I_{\theta}&\rightarrow& \mathcal{DS}_{\theta} = \{DS_2, DS_3\}, m_{\theta} = 1  \\\end{aligned}

The linked-interaction map is :

.. math::

   \begin{aligned}
   I_{\mu} &\rightarrow& I_{\theta}, commonDS = DS_3 \\
   I_{\theta} &\rightarrow&I_{\mu}, commonDS = DS_3 \\\end{aligned}

And:

.. math::

   \begin{aligned}
   M &=& 4, N = \displaystyle{\sum_{i=1}^{3} n_i} \\
   \mathcal{I}_1 &=& \{I_{\mu} \}\\
   \mathcal{I}_2 &=& \{I_{\theta}\} \\
   \mathcal{I}_3 &=& \{I_{\mu}, I_{\theta}\} \\\end{aligned}

.. math::

   \begin{aligned}
   y_1 = \left[\begin{array}{ccc} 
   C_1^1 & C_1^3 \end{array}\right]
   \left[\begin{array}{c}
   x_1 \\
   x_3 
   \end{array}\right]
   + D_1\lambda_1 \\
   y_2 = \left[\begin{array}{ccc} 
   C_2^2 & C_2^3 \end{array}\right]
   \left[\begin{array}{c}
   x_2 \\
   x_3 
   \end{array}\right]
   + D_2\lambda_2 \end{aligned}

.. math::

   \begin{aligned}
   \left[\begin{array}{c}
   R_1 \\
   R_2 \\
   R_3 \end{array}\right]=
   \left[\begin{array}{c}
   B_1^1\lambda_1  \\
   B_2^2\lambda_2  \\
   B_1^3\lambda_1 + B_2^3\lambda_2
   \end{array}\right]\end{aligned}

.. math::

   \begin{aligned}
   M_{OSNSP} &=& \left[\begin{array}{cc} 
   D_1+hC_1^1W_1B_1^1+hC_1^3W_3B_1^3 & hC_1^3W_3B_2^3 \\
   hC_2^3W_3B_1^3 & D_2+hC_2^2W_2B_2^2+hC_2^3W_3B_2^3 
   \end{array}\right]\left[\begin{array}{c} 
   \lambda_1  \\
   \lambda_2
   \end{array}\right]_{k+1} \end{aligned}

relative degree
~~~~~~~~~~~~~~~

Let us consider the global vector

.. math::

   \begin{aligned}
   Y =\left[\begin{array}{c} 
   y_1 \\
   ...  \\
   y_M
   \end{array}\right] = CX + D\Lambda\end{aligned}

We denote by :math:`r_j` the relative degree of equation :math:`j`,
:math:`j\in [1..M]`. We have:

.. math::

   \begin{aligned}
   y_j = \displaystyle{\sum_{i=1}^n C_j^i x_i +D_{j,j}\lambda_j + \sum_{i\neq j, i=1}^m D_{j,i} \lambda_i } \end{aligned}

| :math:`D_{j,i}` a scalar and :math:`C_j^i` a :math:`1 \times n_i`
line-vector.
| If :math:`D_{jj} \neq 0`, then :math:`r_j=0`. Else, we should consider
the first derivative of :math:`y_j`.
| Before that, recall that:

.. math::

   \begin{aligned}
   R_i = \displaystyle{\sum_{k=1}^M B_k^i \lambda_j}\end{aligned}

| Through many of the :math:`B_j^i` are equal to zero, we keep them all
in the following lines.
| Then:

.. math::

   \begin{aligned}
   \dot y_j &=& \displaystyle{\sum_{i=1}^n C_j^i (A_i x_i +  \sum_{k=1}^M B_k^i \lambda_k  ) + f(\lambda_k)_{k\neq j}} \\
   &=& \displaystyle{\sum_{i=1}^n C_j^i (A_i x_i + B_j^i \lambda_j + \sum_{k=1,k\neq j}^M B_k^i \lambda_k  ) + \ldots}\end{aligned}

| So, if :math:`\displaystyle{\sum_{i=1}^n C_j^i B_j^i} \neq 0` (note
that this corresponds to the product between line :math:`j` of :math:`C`
and column :math:`j` of :math:`B`) then :math:`r_j=1` else we consider
the next derivative, and so on.
| In derivative :math:`r`, the coefficient of :math:`\lambda_j` will be:

.. math::

   \begin{aligned}
   coeff_j&=& \displaystyle{\sum_{i=1}^n C_j^i (A_i)^{r-1} B_j^i }\end{aligned}

if :math:`coeff_j\neq 0` then :math:`r_j = r`.

LagrangianDS - Lagrangian Linear Relations
------------------------------------------

General notations
~~~~~~~~~~~~~~~~~

We consider :math:`n` dynamical systems, lagrangian and non linear, of
the form:

.. math:: M_i(q_i) \ddot q_i + N_i(\dot q_i, q_i) = F_{Int,i}(\dot q_i , q_i , t)+F_{Ext,i}(t) + p_i

| Each system if of dimension :math:`n_i`, and we denote
:math:`N = \displaystyle{\sum_{i=1}^{n} n_i}`.
| An interaction, :math:`I_{\alpha}` is composed with a non smooth law,
:math:`nslaw_{\alpha}` and a relation:

.. math:: y_{\alpha} = H_{\alpha}Q_{\alpha} + b_{\alpha}

| The “dimension” of the interaction, ie the size of vector
:math:`y_{\alpha}`, is denoted :math:`m_{\alpha}` and we set:

.. math:: M_y = \sum_{\alpha=1}^{m} m_{\alpha}

:math:`m` being the number of interactions in the Non Smooth Dynamical
System.
| :math:`Q_{\alpha}` is a vector that represents the DS concerned by the
interaction. Its dimension is noted :math:`N_{\alpha}`, this for
:math:`n_{\alpha}` systems in the interaction.
| :math:`H_{\alpha}` is a :math:`m_{\alpha} \times N_{\alpha}`
row-blocks matrix and :math:`b_{\alpha}` a :math:`m_{\alpha}` vector.

.. math::

   H_{\alpha}=\left[\begin{array}{ccc} 
   H_{\alpha}^i & H_{\alpha}^j & ...\end{array}\right]

| with :math:`i,j,...\in \mathcal{DS}_{\alpha}` which is the set of DS
belonging to interaction :math:`\alpha`.
| We also have the following relation:

.. math::

   \left[\begin{array}{c} 
   R_{\alpha}^i \\
   R_{\alpha}^j \\
   ...  
   \end{array}\right] = {}^tH_{\alpha}\lambda_{\alpha}
   =\left[\begin{array}{c} 
   {}^tH_{\alpha}^i \\
   {}^tH_{\alpha}^j \\
   ...
   \end{array}\right]\lambda_{\alpha}

| :math:`R_{\alpha}^i` represents the contribution of interaction
:math:`\alpha` on the reaction of the dynamical system :math:`i`, and
:math:`{}tH_{\alpha}^i` is a :math:`n_i \times m_{\alpha}` block matrix.
| And so:

.. math:: R_i = \sum_{\beta\in\mathcal{I}_i}R_{\beta}^i=\sum_{\beta\in\mathcal{I}_i}{}H^i_{\beta} \lambda_{\beta}

| with :math:`\mathcal{I}_i` the set of interactions in which dynamical
system number :math:`i` is involved.
| Introducing the time dicretisation, we get:

.. math::

   \begin{aligned}
   \dot q_i^{k+1} = \dot q_{free,i} + W_iR_i^{k+1}
   \nonumber\\
   \dot y_{\alpha}^{k+1} = H_{\alpha}\dot Q_{\alpha}^{k+1} \\
   \nonumber\\
   R_i^{k+1} = \sum_{\beta\in\mathcal{I}_i}H^i_{\beta} \lambda_{\beta}^{k+1}\end{aligned}

ie,

.. math::

   \begin{aligned}
     y_{\alpha}^{k+1} &=& H_{\alpha} Q_{\alpha}^{free} + H_{\alpha}W_{\alpha}{}^tH_{\alpha}\lambda_{\alpha}+\sum_{i\in \mathcal{DS}_{\alpha}}\sum_{\beta\in\mathcal{I}_i,\alpha\neq\beta}H_{\alpha}^iW_iH_{\beta}^j\lambda_{\beta}\end{aligned}

with :math:`W_{\alpha}` given by .

The global-assembled :math:`Y` vector, of dimension M, composed by
:math:`m` :math:`y_{\alpha}` subvectors, is given by:

.. math::

   \begin{aligned}
   Y_{k+1} = q_{OSNSP} + M_{OSNSP}\Lambda_{k+1}\end{aligned}

with:

.. math::

   \begin{aligned}
   q_{OSNSP}^{\alpha} = H_{\alpha} Q_{\alpha}^{free}\end{aligned}

and for :math:`M_{OSNSP}`, the block-diagonal term for block-row
:math:`\alpha` is

.. math:: \sum_{j\in \mathcal{DS}_{\alpha}}H_{\alpha}^jW_j{}^tH_{\alpha}^j

an :math:`m_{\alpha}\times m_{\alpha}` square matrix. The extra-diagonal
block term, in position (:math:`\alpha,\beta`) is:

.. math:: \sum_{j\in \mathcal{DS}_{\alpha}\cap\mathcal{DS}_{\beta}}H_{\alpha}^jW_j{}^tH_{\beta}^j

| and is a :math:`m_{\alpha}\times m_{\beta}` matrix. This matrix
differs from 0 when interactions :math:`\alpha` and :math:`\beta` are
coupled, ie have common DS.

Or, for the relation l of interaction :math:`\alpha`, we get:

.. math:: \sum_{j\in \mathcal{DS}_{\alpha}}H_{\alpha,l}^jW_j{}^tH_{\alpha}^j

for the diagonal, and

.. math:: \sum_{j\in \mathcal{DS}_{\alpha}\cap\mathcal{DS}_{\beta}}H_{\alpha,l}^jW_j{}^tH_{\beta}^j

| for extra-diagonal terms.
| :math:`H_{\alpha,l}`, row number :math:`l` of :math:`H_{\alpha}`.

WARNING: depending on linear and non linear case for the DS, there
should be a factor h ahead W. See Bouncing Ball template.

Block matrix approach
---------------------

The built of the OSNSProblem matrix could be computed using block matrix
structure. This section describe these matrices. It is not implemented
in Siconos. Using previous notations, :math:`n` is the number of DS.
:math:`m` the number of interations.

Block matrix of DS
~~~~~~~~~~~~~~~~~~

.. math:: \boldsymbol{M}  \boldsymbol{\dot X}=\boldsymbol{A} \boldsymbol{X} + \boldsymbol{R}

where :math:`\boldsymbol{M}=diag(M_1,...M_n)` and
:math:`\boldsymbol{A}=diag(A_1,..,A_n)`.

.. math:: \boldsymbol{R}=\boldsymbol{B}\boldsymbol{\lambda}

.. math::

   \boldsymbol{B}=\left( \begin{array} {c} B^1_{1}...B^1_m\\.\\.\\
       B^n_1...B^n_m  \end{array}\right)

 Some of :math:`B^i_j` doesn’t exist.

Block matrix of interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \boldsymbol{Y}= \boldsymbol{C}  \boldsymbol{X}+
   \boldsymbol{D} \boldsymbol{\lambda}

with :math:` \boldsymbol{D}=diag(D_1..D_m)` and

.. math::

   \boldsymbol{C}=\left( \begin{array} {c}
       C^1_{1}..C^n_1\\.\\.\\C^1_{m}...C^n_{m} \end{array}\right)

 Some of :math:`C^i_j` doesn’t exist.

OSNSProblem using block matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Matrix of the OSNS Problem could be assembled using the following
block-product-matrices
:math:`\boldsymbol{C}\boldsymbol{W}\boldsymbol{B}`.

Dynamical Systems formulations in Siconos.
==========================================

+-----------+------------------+
| author    | F. Pérignon      |
+===========+==================+
| date      | March 22, 2006   |
+-----------+------------------+
| version   | Kernel 1.1.4     |
+-----------+------------------+

Class Diagram
-------------

There are four possible formulation for dynamical systems in Siconos,
two for first order systems and two for second order Lagrangian systems.
The main class is DynamicalSystem, all other derived from this one, as
shown in the following diagram:

|image| [DSDiagram]

General non linear first order dynamical systems
:math:`\rightarrow` class *DynamicalSystem*
------------------------------------------------

| This is the top class for dynamical systems. All other systems classes
derived from this one.

A general dynamical systems is described by the following set of
:math:`n` equations, completed with initial conditions:

.. math::

   \begin{aligned}
     \dot x &=& f(x,t) + T(x) u(x, \dot x, t) + r \\
     x(t_0)&=&x_0 \end{aligned}

-  :math:`x`: state of the system - Vector of size :math:`n`.

-  :math:`f(x,t)`: vector field - Vector of size :math:`n`.

-  :math:`u(x, \dot x, t)`: control term - Vector of size :math:`uSize`.

-  :math:`T(x)`: :math:`n\times uSize` matrix, related to control term.

-  :math:`r`: input due to non-smooth behavior - Vector of size
   :math:`n`.

| The Jacobian matrix, :math:`\nabla_x f(x,t)`, of :math:`f` according
to :math:`x`, :math:`n\times n` square matrix, is also a member of the
class.

| Initial conditions are given by the member :math:`x_0`, vector of size
:math:`n`. This corresponds to x value when simulation is starting, ie
after a call to strategy->initialize().

There are plug-in functions in this class for :math:`f` (vectorField),
:math:`jacobianX`, :math:`u` and :math:`T`. All of them can handle a
vector of user-defined parameters.

First order linear dynamical systems :math:`\rightarrow` class *LinearDS*
-------------------------------------------------------------------------

Derived from DynamicalSystem, described by the set of :math:`n`
equations and initial conditions:

.. math::

   \begin{aligned}
     \dot x &=& A(t)x(t)+Tu(t)+b(t)+r \\
     x(t_0)&=&x_0 \end{aligned}

With:

-  :math:`A(t)`: :math:`n\times n` matrix, state independent but
   possibly time-dependent.

-  :math:`b(t)`: Vector of size :math:`n`, possibly time-dependent.

| Other variables are those of DynamicalSystem class.
| :math:`A` and :math:`B` have corresponding plug-in functions.

| Warning: time dependence for :math:`A` and :math:`b` is not available
at the time in the simulation part for this kind of dynamical systems.

Links with vectorField and its Jacobian are:

.. math::

   \begin{aligned}
     f(x,t) &=& A(t)x(t)+b(t) \\
     jacobianX&=&\nabla_x f(x,t) = A(t) \end{aligned}

Second order non linear Lagrangian dynamical systems
:math:`\rightarrow` class *LagrangianDS*
----------------------------------------------------

Lagrangian second order non linear systems are described by the
following set of\ :math:`nDof` equations + initial conditions:

.. math::

   \begin{aligned}
    M(q) \ddot q + NNL(\dot q, q) + F_{Int}(\dot q , q , t) &=& F_{Ext}(t) + p \\
    q(t_0) &=& q0 \\
    \dot q(t_0) &=& velocity0 \end{aligned}

With:

-  :math:`M(q)`: :math:`nDof\times nDof` matrix of inertia.

-  :math:`q`: state of the system - Vector of size :math:`nDof`.

-  :math:`\dot q` or :math:`velocity`: derivative of the state according
   to time - Vector of size :math:`nDof`.

-  :math:`NNL(\dot q, q)`: non linear terms, time-independent - Vector
   of size :math:`nDof`.

-  :math:`F_{Int}(\dot q , q , t)`: time-dependent linear terms - Vector
   of size :math:`nDof`.

-  :math:`F_{Ext}(t)`: external forces, time-dependent BUT do not depend
   on state - Vector of size :math:`nDof`.

-  :math:`p`: input due to non-smooth behavior - Vector of size
   :math:`nDof`.

The following Jacobian are also member of this class:

-  jacobianQFInt = :math:`\nabla_q F_{Int}(t,q,\dot q)` -
   :math:`nDof\times nDof` matrix.

-  jacobianVelocityFInt = :math:`\nabla_{\dot q} F_{Int}(t,q,\dot q)` -
   :math:`nDof\times nDof` matrix.

-  jacobianQNNL = :math:`\nabla_q NNL(q,\dot q)` -
   :math:`nDof\times nDof` matrix.

-  jacobianVelocityNNL = :math:`\nabla_{\dot q}NNL(q,\dot q)` -
   :math:`nDof\times nDof` matrix.

| There are plug-in functions in this class for :math:`F_{int}`,
:math:`F_{Ext}`, :math:`M`, :math:`NNL` and the four Jacobian matrices.
All of them can handle a vector of user-defined parameters.

Links with first order dynamical system are:

.. math::

   \begin{aligned}
     n &= &2nDof \\
     x &=&\left[\begin{array}{c}q \\ \dot q \end{array}\right] \\
     f(x,t) &=&  \left[\begin{array}{c} \dot q \\ M^{-1}(F_{Ext}-F_{Int}-NNL) \end{array}\right] \\
     \\
     \nabla_x f(x,t) &=& 
     \left[\begin{array}{cc} 
         0_{nDof\times nDof} & I_{nDof\times nDof} \\
         \nabla_q(M^{-1})(F_{Ext}-F_{Int}-NNL) -M^{-1}\nabla_q(F_{Int}+NNL) &  -M^{-1}\nabla_{\dot q}(F_{Int}+NNL) 
       \end{array}\right] \\
     r &=& \left[\begin{array}{c} 0_{nDof} \\ p \end{array}\right] \\
     u(x,\dot x,t) &=& u_L(\dot q, q, t) \text{  (not yet implemented)} \\
     T(x) &=& \left[\begin{array}{c} 0_{nDof} \\ T_L(q) \end{array}\right] \text{  (not yet implemented)} \\\end{aligned}

| With :math:`0_{n}` a vector of zero of size :math:`n`,
:math:`0_{n\times m}` a :math:`n\times m` zero matrix and
:math:`I_{n\times n}`, identity :math:`n\times n` matrix.

Warning: control terms (:math:`Tu`) are not fully implemented in
Lagrangian systems. This will be part of future version.

Second order linear and time-invariant Lagrangian dynamical systems :math:`\rightarrow` class *LagrangianLinearTIDS*
--------------------------------------------------------------------------------------------------------------------

.. math::

   \begin{aligned}
   M \ddot q + C \dot q + K q =  F_{Ext}(t) + p\end{aligned}

With:

-  :math:`C`: constant viscosity :math:`nDof\times nDof` matrix

-  :math:`K`: constant rigidity :math:`nDof\times nDof` matrix

And:

.. math::

   \begin{aligned}
   F_{Int} &=& C \dot q + K q \\
   NNL &=& 0_{nDof} \end{aligned}

Dynamical Systems implementation in Siconos.
============================================

+-----------+--------------------+
| author    | F. Pérignon        |
+===========+====================+
| date      | November 7, 2006   |
+-----------+--------------------+
| version   | Kernel 1.3.0       |
+-----------+--------------------+

Introduction
------------

| This document is only a sequel of notes and remarks on the way
Dynamical Systems are implemented in Siconos.
| It has to be completed, reviewed, reorganized etc etc for a future
Developpers’Guide.
| See also documentation in Doc/User/DynamicalSystemsInSiconos for a
description of various dynamical systems types.

Class Diagram
-------------

There are four possible formulation for dynamical systems in Siconos,
two for first order systems and two for second order Lagrangian systems.
The main class is DynamicalSystem, all other derived from this one, as
shown in the following diagram:

|image| [DSDiagram]

Construction
------------

Each constructor must:

-  initialize all the members of the class and of the top-class if it
   exists

-  allocate memory and set value for all required inputs

-  allocate memory and set value for optional input if they are given as
   argument (in xml for example)

-  check that given data are coherent and that the system is complete
   (for example, in the LagrangianDS if the internal forces are given as
   a plug-in, their Jacobian are also required. If they are not given,
   this leads to an exception).

| No memory allocation is made for unused members :math:`\Rightarrow`
requires if statements in simulation. (if!=NULL ...).

DynamicalSystem
~~~~~~~~~~~~~~~

| **Required data:**
| n, x0, f, jacobianXF
| **Optional:**
| T,u

| **Always allocated in constructor:**
| x, x0, xFree, r, rhs, jacobianXF

Warning: default constructor is always private or protected and apart
from the others and previous rules or remarks do not always apply to it.
This for DS class and any of the derived ones.

LagrangianDS
~~~~~~~~~~~~

| **Required data:**
| ndof, q0, velocity0, mass
| **Optional:**
| fInt and its Jacobian, fExt, NNL and its Jacobian.

| **Always allocated in constructor:**
| mass, q, q0, qFree, velocity, velocity0, velocityFree, p.
| All other pointers to vectors/matrices are set to NULL by default.
| Memory vectors are required but allocated during call to initMemory
function.

Various rules:

-  fInt (NNL) given as a plug-in :math:`\Rightarrow` check that
   JacobianQ/Velocity are present (matrices or plug-in)

-  any of the four Jacobian present :math:`\Rightarrow` allocate memory
   for block-matrix jacobianX (connectToDS function)

-  

| check: end of constructor or in initialize?
| computeF and JacobianF + corresponding set functions: virtual or not?

Specific flags or members
-------------------------

-  isAllocatedIn: to check inside-class memory allocation

-  isPlugin: to check if operators are computed with plug-in or just
   directly set as a matrix or vector

-  workMatrix: used to save some specific matrices in order to avoid
   recomputation if possible (inverse of mass ...)

plug-in management
------------------

| DynamicalSystem class has a member named parameterList which is a
:math:`map<string, SimpleVector*>`, ie a list of pointers to
SimpleVector\*, with a string as a key to identified them. For example,
:math:`parametersList["mass"]` is a SimpleVector\*, which corresponds to
the last argument given in mass plug-in function.
| By default, each parameters vectors must be initialized with a
SimpleVector of size 1, as soon as the plug-in is declared. Moreover, to
each vector corresponds a flag in isAllocatedIn map, to check if the
corresponding vector has been allocated inside the class or not.
| For example, in DynamicalSystem, if
:math:`isPlugin["vectorField"]==true`, then, during call to constructor
or set function, it is necessary to defined the corresponding parameter:
| :math:`parametersList["vectorField"] = new SimpleVector(1)`
| and to complete the :math:`isAllocatedIn` flag:
| :math:`isAllocatedIn["parameter_for_vectorField"] = true`.

Interactions
============

+-----------+--------------------+
| author    | F. Pérignon        |
+===========+====================+
| date      | November 7, 2006   |
+-----------+--------------------+
| version   | Kernel 1.3.0       |
+-----------+--------------------+

Introduction
------------

| This document is only a sequel of notes and remarks on the way
Interactions are implemented in Siconos.
| It has to be completed, reviewed, reorganized etc etc for a future
Developpers’Guide.
| See also documentation in Doc/User/Interaction.

Class Diagram
-------------

Description
-----------

review of interactions (for EventDriven implementation) 17th May 2006.

| variable *nInter* renamed in *interactionSize*: represents the size of
*y* and *:math:`\lambda`*. NOT the number of relations !!

| add a variable *nsLawSize* that depends on the non-smooth law type.
| Examples:

NewtonImpact -> *nsLawSize* = 1

Friction 2D -> *nsLawSize* = 2

Friction 3D -> *nsLawSize* = 3

...

| *nsLawSize* = n with n dim of matrix D in : :math:`y=Cx+D\lambda`, D
supposed to be a full-ranked matrix.
| Warning: this case is represented by only one relation of size n.

*numberOfRelations*: number of relations in the interaction,
*numberOfRelations* =
:math:`{\displaystyle \frac{{\textit{interactionSize}}}{{\textit{nsLawSize}}}}`.

Notes on the Non Smooth Dynamical System construction
=====================================================

+-----------+--------------------+
| author    | F. Pérignon        |
+===========+====================+
| date      | November 7, 2006   |
+-----------+--------------------+
| version   | Kernel 1.3.0       |
+-----------+--------------------+

Introduction
------------

Class Diagram
-------------

Description
-----------

Objects must be constructed in the following order:

DynamicalSystems

NonSmoothLaw: depends on nothing

Relation: no link with an interaction during construction, this will be
done during initialization.

| Interaction: default constructor is private and copy is forbidden. Two
constructors: xml and from data. Required data are a DSSet, a
NonSmoothLaw and a Relation (+ dim of the Interaction and a number).
| Interaction has an initialize function which allocates memory for y
and lambda, links correctly the relation and initializes it .... This
function is called at the end of the constructor. That may be better to
call it in simulation->initialize? Pb: xml constructor needs memory
allocation for y and lambda if they are provided in the input xml file.

NonSmoothDynamicalSystem: default is private, copy fobidden. Two
constructors: xml and from data. Required data are the DSSet and the
InteractionsSet. The topology is declared and constructed (but empty)
during constructor call of the nsds, but initialize in the Simulation,
this because it can not be initialize until the nsds has been fully
described (ie this to allow user to add DS, Inter ...) at any time in
the model, but before simulation initialization).

misc
----

no need to keep a number for Interactions? Only used in xml for OSI, to
know which Interactions it holds.

pb: the number of saved derivatives for y and lambda in Interactions is
set to 2. This must depends on the relative degree which is computes
during Simulation initialize and thus too late. It is so not available
when memory is allocated (Interaction construction). Problem-> to be
reviewed.

OneStepIntegrator and derived classes.
======================================

+-----------+--------------------+
| author    | F. Pérignon        |
+===========+====================+
| date      | November 7, 2006   |
+-----------+--------------------+
| version   | Kernel 1.3.0       |
+-----------+--------------------+

Introduction
------------

| This document is only a sequel of notes and remarks on the way
OneStepIntegrators are implemented in Siconos.
| It has to be completed, reviewed, reorganized etc etc for a future
Developpers’Guide.
| See also documentation in Doc/User/OneStepIntegrator for a description
of various OSI.

Class Diagram
-------------

Misc
----

OSI review for consistency between Lsodar and Moreau:

-  add set of DynamicalSystem\*

-  add set of Interaction\*

-  add link to strategy that owns the OSI

-  remove td object in OSI -> future: replace it by a set of td (one per
   ds)

-  add strat in constructors arg list

| osi -> strat -> Model -> nsds -> topology
| osi -> strat -> timeDiscretisation

| let a timeDiscretisation object in the OSI? set of td (one per ds)?
| create a class of object that corresponds to DS on the simulation side
?
| will contain the DS, its discretization, theta for Moreau ... ?
| Allow setStrategyPtr operation? Warning: need reinitialisation.

| Required input by user:

-  list of DS or list of Interactions ?

-  pointer to strategy

-  ...

Construction
------------

Each constructor must:

-  

Moreau
~~~~~~

| Two maps: one for W, and one for theta. To each DS corresponds a theta
and a W.
| Strategy arg in each constructor.

| **Required data:**

| **Optional:**

| **Always allocated in constructor:**

Warning: default constructor is always private or protected and apart
from the others and previous rules or remarks do not always apply to it.

Lsodar
~~~~~~

| **Required data:**

| **Optional:**

| **Always allocated in constructor:**

First Order Nonlinear Relation 
===============================

+-----------+----------------+
| author    | 0. Bonnefon    |
+===========+================+
| date      | July, 1 2009   |
+-----------+----------------+
| version   | Kernel 3.0.0   |
+-----------+----------------+

Computation of the number of Index Set and various levels
=========================================================

+-----------+----------------------+
| author    | V. Acary             |
+===========+======================+
| date      | Septembre 16, 2011   |
+-----------+----------------------+
| version   | Kernel 3.3.0         |
+-----------+----------------------+

In this chapter, we give some hints on the automatic computation of the
number of index sets, the number of derivatives in the Interaction and
the levelMin and LevelMax.

Why is the relative degree not relevant ?
-----------------------------------------

In this section, we give a very brief overview of the notion of relative
degree.

First order Linear complementary systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A Linear Complementarity System (LCS) is defined by

.. math::

   \label{eq:LCS-bis}
     \begin{cases}
       \dot x = A x +B \lambda \\
        y = C x + D \lambda\\
       0 \leq  y \perp \lambda \geq 0 \\
     \end{cases}

[Relative degree in the SISO case] Let us consider a linear system in
state representation given by the quadruplet
:math:`(A,B,C,D) \in {\mbox{\rm $I\!\!R$}}^{n\times n}\times{\mbox{\rm $I\!\!R$}}^{n \times m}\times {\mbox{\rm $I\!\!R$}}^{m\times n}\times{\mbox{\rm $I\!\!R$}}^{m\times m} `:

.. math::

   \label{eq:LS}
           \begin{cases}
             \dot x = A x +B \lambda \\
             y = C x + D \lambda
           \end{cases}

-  In the Single Input/ Single Output (SISO) case (:math:`m=1`), the
   relative degree is defined by the first non zero Markov parameters :

   .. math::

      \label{eq:Markov-Parameter}
                D, CB, CAB, CA^2B, \ldots, CA^{r-1}B, \ldots

-  In the multiple input/multiple output (MIMO) case (:math:`m>1`), an
   *uniform* relative degree is defined as follows. If :math:`D` is non
   singular, the relative degree is equal to :math:`0`. Otherwise, it is
   assumed to be the first positive integer :math:`r` such that

   .. math::

      \label{eq:mimo-r}
                CA^{i}B =0, \quad i=0\ldots q-2

   while

   .. math::

      \label{eq:mimo-r2}
                CA^{r-1}B \text{ is non singular}.

The Markov parameters arise naturally when we derive with respect to
time the output :math:`y`,

.. math::

   \begin{aligned}
         \label{eq:y-derive}
         y &=& C x + D \lambda \\
         \dot y &=& CA x + CB \lambda, \text{ if } D= 0  \\
         \ddot y &=& CA^2 x + CAB \lambda, \text{ if }  D=0, CB=0\\
         &\ldots& \\
         y^{(r)} &=& CA^{r} x + CA^{r-1}B \lambda, \text{ if } D=0, CB=0, CA^{r-2}B=0, r=1\ldots r-2 \\
         &\ldots&
       \end{aligned}

and the first non zero Markov parameter allows us to define the output
:math:`y` directly in terms of the input :math:`\lambda`.

In continuous time, the nature of solutions depends strongly on the
relative degree. When we want to perform the time–integration of such
systems, we need also to reduce the relative degree or to known it to
correctly operate.

We can observe that the relative degree :math:`0` is well defined only
by the relation (:math:`D` nonsingular) and by the nonsmooth law.
Indeed, let us imagine that the nonsmooth law is defined by
:math:`0\leq\dot y \perp \lambda \geq 0 `. We can easily see that the
relative degree is reduced.

In the MIMO, the computation of non uniform relative degree is hard
task. This is also the case for nonlinear systems.

Second order Lagrangian systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us consider a second order linear and time-invariant Lagrangian
dynamical system (see § [Sec:LagrangianLineatTIDS])

.. math::

   \label{eq:rd1}
     \begin{cases}
       M \dot v + C v + K q = F_{Ext}(t) + p \\
       \dot q = v
     \end{cases}

together with a Lagrangian linear relation

.. math::

   y= Cq + e + D \lambda + Fz,
     \label{eq:rd2}

.. math::

   p = C^t \lambda
   \label{eq:rd3}

and a simple nonsmooth law,

.. math::

   0\leq y \perp \lambda \geq 0
   \label{eq:rd4}

If :math:`D>0`, the relative degree is uniformly zero and the system can
be solved without deriving the output ([eq:rd2]). Indeed, we known that
the solution of the LCP

.. math::

   0\leq Cq + e + D \lambda + Fz, \perp \lambda \geq 0
   \label{eq:rd5}

is unique and Lipschitz with respect to :math:`q`. It can be denoted as
:math:`\lambda(q) = \mbox{SOL}(D,Cq + e +Fz)`. Therefore, the
differential equation ([eq:rd1]) reduces to a standard ODE with a
Lipschitz RHS

.. math::

   \label{eq:rd6}
     \begin{cases}
       M \dot v + C v + K q = F_{Ext}(t) + C^t \lambda(q)  \\
       \dot q = v
     \end{cases}

In the case that we deal with unilateral contact, we usually have
:math:`D=0` and the relative degree of the system is :math:`2`. In this
case, the output has to be differentiated as

.. math::

   \label{eq:rd7}
      \dot y= C \dot q,

and an impact law has to added, for instance the newton’s impact law

.. math::

   \label{eq:rd8}
     \text{ if } y=0, \text{ when } \dot y^+= -e y^-

In the same vein, the equations of motion ([eq:rd1]) is not sufficient
since the velocity may encounter jumps. The dynamics is usually replaced
by a measure differential equation of the form

.. math::

   \label{eq:rd10}
     \begin{cases}
       M dv + C v^+(t) dt + K q(t) dt = F_{Ext}(t)dt + di \\
       \dot q = v
     \end{cases}

where :math:`di` is the measure that can be related to :math:`p` thanks
to

.. math::

   \label{eq:rd11}
     di = p dt + \sigma \delta _{t^*}

is only one jump is expected at :math:`{t^*}`.

Conclusion for the implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the continuous time mathematical analysis, the relative degree is
very important to know if we have to compute the derivatives of the
output :math:`y^{(n)}` and to consider various levels for the input
:math:`p , \sigma, ....`

However in the numerical practice, the time –discretization makes an
assumption on the relative degree and treats the nonsmooth law at
different levels. The resulting time discretized system posseses more or
less variables.

Consider for instance ([eq:rd1]) in the case of the Moreau scheme

| [eq:MoreauTS] M(v\ :sub:`k+1`-v:sub:`k`) + h (K q\ :sub:`k+`\ + C
v\ :sub:`k+`) = p\ :sub:`k+1` = G(q\ :sub:`k+1`) :sub:`k+1`,
| q\ :sub:`k+1` = q\ :sub:`k` + h v\ :sub:`k+`,
| y\ :sub:`k+1` = G\ :sup:``\ (q:sub:`k+1`) v\ :sub:`k+1`

| l \|y\ :sup:``\ :sub:`k+1` 0 0 y\ :sup:``\ :sub:`k+1` + e
y\ :sup:``\ :sub:`k` :sup:``\ :sub:`k+1` 0,
|  :sup:``\ :sub:`k+1` =0.[eq:MoreauTSd]

, I

and the Schatzman–Paoli scheme

| M(q\ :sub:`k+1`-2q:sub:`k`\ +q\ :sub:`k-1`) + h\ :sup:`2` (K
q\ :sub:`k+`\ + C v\ :sub:`k+`) = p\ :sub:`k+1`,
| v\ :sub:`k+1`\ =,
| y\ :sub:`k+1` = h()
| p\ :sub:`k+1`\ = G() :sub:`k+1`
| 0 y\ :sub:`k+1` :sub:`k+1` 0 .

We can see easily that the number of derivatives (or levels) that we
store for :math:`y` and :math:`\lambda` is independent on the relative
degree but is chosen by the OneStepIntegrator with respect to the type
of systems.

How to define and compute the various levels and the number of indexSets 
-------------------------------------------------------------------------

:math:`y` related variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The size of the vector y in the Interaction depends on

-  the OneStepIntegrator type.

   -  see the difference between the Moreau and Schatzman Paoli scheme,

   -  plan the time–discontinuous Galerkin scheme

   -  plan the Higher Order Moreau sweeping process (HOSP)

-  the Simulation type.

   -  In Timestepping or Event-driven we do not need the same number of
      stored :math:`y`

-  the NonSmoothLaw type.

   -  If we consider some cases with or without friction in Timestepping
      or Event-driven, we need to adapt the number of stored :math:`y`

Since the various levels of y are used to build the index sets we will
need from :math:`0` to a computed size that depends on the previous
criteria. Only a part will be used in the OneStepNSProblem.

:math:`\lambda` related variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The size of the vector lambda in the Interaction depends on the same
criteria than in the previous section. Only, the number of lambda is not
the same as y since a multiplier lambda[i] is not necessarily related to
y[i]

Rules for implementation
------------------------

We can define new members in Interaction:

-  \_lowerlevelForOutput, this value is to :math:`0` by default

-  \_upperlevelForOutput, this value must be computed at initialization
   with respect to the previous criteria

-  \_lowerlevelForInput, this value must be computed at initialization
   with respect to the previous criteria

-  \_upperlevelForInput, this value must be computed at initialization
   with respect to the previous criteria

This level are computed in Simulation::ComputeLevelsForInputAndOutput. A
visitor is used for the OneStepIntegrator. Furthermore, four global
levels are computed

-  \_levelMinForOutput this value is the minimum level for the output
   Interaction::\_lowerlevelForOutput for all the interactions

-  \_levelMaxForOutput this value is the maximum level for the output
   Interaction::\_upperlevelForOutput for all the interactions

-  \_levelMinForInput this value is the minimum level for the output
   Interaction::\_lowerlevelForInput for all the interactions

-  \_levelMaxForInput this value is the maximum level for the output
   Interaction::\_upperlevelForInput for all the interactions

-  the values y[i] must be initialized from \_lowerlevelForOutput to
   \_upperlevelForOutput.

-  the values lamdba[i] must be initialized from \_lowerlevelForInput to
   \_upperlevelForInput.

-  the values y[i] in Interaction must be used in priority to store the
   i-th derivative of :math:`y`. When it is needed, higher index
   :math:`i` can be used for other triggering variables. For instance,
   for an Event–Driven scheme with a Lagrangian systems with friction,
   sliding velocity must be stored.

-  the values of lamdba[i] must stored the various multiplier for the
   nonsmooth law. We affect the same index :math:`i` as for the level of
   y[i] present in the corresponding nonsmooth law.

-  The number of IndexSets should follows \_levelMaxForY.

For the dynamical systems :

-  The number of levels for \_r and \_p in the DS should follow
   \_lowerlevelForInput and \_upperlevelForOutput of the associated
   interactions. This is done in Interaction::initialize.

-  A new variable should be added in the LagrangianDS to store the
   multiplier at the position level (\_tau ?) to avoid the use of
   \_p[0]. Indeed, we will continue to assume that \_p is the input in
   the equation of motion. For lambda we can use lambda[0]

TODO LIST AND QUESTIONS

-  What about the case of multiples interactions on a DS with various
   \_lowerlevelForInput and \_upperlevelForOutput ? Normally, all the
   levels should be correctly initialized thanks to the proposed
   implementation (r2821)

-  DynamicalSystem::\_r should be a VectorOfVectors

-  DynamicalSystem::\_r is split in LagrangianDS. a first part is
   LagrangianDS::\_p. The other is not implemented !!
   LagrangianDS::\_tau ?

Newton’s linearization for First Order Systems
==============================================

+---------------+------------------------+
| author        | O.Bonnefon, V. Acary   |
+===============+========================+
| date          | Sept, 07, 2007         |
+---------------+------------------------+
| last update   | Feb, 2011              |
+---------------+------------------------+
|               | April, 2014            |
+---------------+------------------------+
| version       |                        |
+---------------+------------------------+

This section is devoted to the implementation and the study of the
algorithm. The interval of integration is :math:`[0,T]`, :math:`T>0`,
and a grid :math:`t_{0}=0`, :math:`t_{k+1}=t_{k}+h`, :math:`k \geq 0`,
:math:`t_{N}=T` is constructed. The approximation of a function
:math:`f(\cdot)` on :math:`[0,T]` is denoted as :math:`f^{N}(\cdot)`,
and is a piecewise constant function, constant on the intervals
:math:`[t_{k},t_{k+1})`. We denote :math:`f^{N}(t_{k})` as
:math:`f_{k}`. The time-step is :math:`h>0`.

Various first order dynamical systems with input/output relations
-----------------------------------------------------------------

FirstOrderR. Fully nonlinear case
'''''''''''''''''''''''''''''''''

Let us introduce the following system,

.. math::

   \begin{array}{l}
   M \dot{x}(t) = f(x(t),t) + r(t)  \\[2mm]
   y(t) = h(t,x(t),\lambda (t)) \\[2mm]
   r(t) = g(t,x(t),\lambda (t) ) \\[2mm]
   \end{array}
   \label{first-DS}

where :math:`\lambda(t) \in {\mbox{\rm $I\!\!R$}}^m` and
:math:`y(t) \in {\mbox{\rm $I\!\!R$}}^m` are complementary variables
related through a multi-valued mapping. According to the class of
systems, we are studying, the function :math:`f` and :math:`g` are
defined by a fully nonlinear framework or by affine functions. We have
decided to present the time-discretization in its full generality and
specialize the algorithms for each cases in Section [Sec:Spec]. This
fully nonlinear case is not implemented in Siconos yet. This fully
general case is not yet implemented in Siconos.

This case is implemented in Siconos with the relation FirstOrderR using
the subtype NonLinearR

FirstOrderType1R
''''''''''''''''

Let us introduce a new notation,

.. math::

   \begin{array}{l}
   M \dot{x}(t) = f(x(t),t) + r(t)  \\[2mm]
   y(t) = h(t,x(t)) \\[2mm]
   r(t) = g(t,\lambda (t) ) \\[2mm]
   \end{array}
   \label{first-DS1}

This case is implemented in Siconos with the relation FirstOrderType1R.

FirstOrderType2R
''''''''''''''''

Let us introduce a new notation,

.. math::

   \begin{array}{l}
   M \dot{x}(t) = f(x(t),t) + r(t)  \\[2mm]
   y(t) = h(t,x(t),\lambda (t)) \\[2mm]
   r(t) = g(t,\lambda (t) ) \\[2mm]
   \end{array}
   \label{first-DS2}

This case is implemented in Siconos with the relation FirstOrderType2R.

Linear case 
''''''''''''

Let us introduce a new notation,

.. math::

   \begin{array}{l}
   M \dot{x}(t) = Ax(t) + r(t)  +b(t)\\[2mm]
   y(t) = h(x(t),\lambda (t),z) = Cx + Fz + D \lambda  \\[2mm]
   r(t) = g(t,\lambda (t) ) = B \lambda \\[2mm]
   \end{array}
   \label{first-DS3}

Time–discretizations
--------------------

Standard :math:`\theta-\gamma` scheme.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now proceed with the time discretization of ([first-DS3]) by a
fully implicit scheme :

.. math::

   \begin{array}{l}
       \label{eq:toto1}
        M x_{k+1} = M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1})
        + h(1-\gamma)r(t_k)  \\[2mm]
        y_{k+1} =  h(t_{k+1},x_{k+1},\lambda_{k+1}) \\[2mm]
        r_{k+1} =  g(t_{k+1},x_{k+1},\lambda_{k+1})\\[2mm]
        \mbox{NsLaw} ( y_{k+1} , \lambda_{k+1})
     \end{array}

where :math:`\theta = [0,1]` and :math:`\gamma \in [0,1]`. As in , we
call the problem the “one–step nonsmooth problem”.

In the Siconos/Kernel module, the use of :math:`\gamma` is activated in
the class EulerMoreauOSI by the boolean \_useGamma.

This time-discretization is slightly more general than a standard
implicit Euler scheme. The main discrepancy lies in the choice of a
:math:`\theta`-method to integrate the nonlinear term. For
:math:`\theta=0`, we retrieve the explicit integration of the smooth and
single valued term :math:`f`. Moreover for :math:`\gamma =0`, the term
:math:`g` is explicitly evaluated. The flexibility in the choice of
:math:`\theta` and :math:`\gamma` allows the user to improve and control
the accuracy, the stability and the numerical damping of the proposed
method. For instance, if the smooth dynamics given by :math:`f` is
stiff, or if we have to use big step sizes for practical reasons, the
choice of :math:`\theta > 1/2` offers better stability with the respect
to :math:`h`.

Full :math:`\theta-\gamma` scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another possible time–discretization is as follows.

.. math::

   \begin{array}{l}
       \label{eq:toto1-ter}
       M x_{k+1} = M x_{k} + h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h r(t_{k+\gamma}) \\[2mm]
       y_{k+\gamma} = h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
       r_{k+\gamma} = g(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma})\\[2mm]
       \mbox{NsLaw} ( y_{k+\gamma} , \lambda_{k+\gamma})
     \end{array}

We call the scheme ([eq:toto1-ter]) the full :math:`\theta-\gamma`
scheme since it uses also the evaluation at :math:`t_{k+\gamma}` for the
relation.

In the Siconos/Kernel module, the time–stepping scheme is activated in
the class EulerMoreauOSI by the boolean \_useGammaForRelation.

Another possibility for the time discretization in the nonlinear case
would be

.. math::

   \begin{array}{l}
       \label{eq:toto1-quat}
       M x_{k+1} = M x_{k} +h f(x_{k+\theta},t_{k+\theta}) + h r(t_{k+\gamma}) \\[2mm]
       y_{k+\gamma} =  h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
       r_{k+\gamma} = g(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma})\\[2mm]
       \mbox{NsLaw} ( y_{k+\gamma} , \lambda_{k+\gamma})
     \end{array}

This scheme has not been yet implemented in Siconos/Kernel.

Newton’s linearization of ([eq:toto1])
--------------------------------------

Due to the fact that two of the studied classes of systems that are
studied in this paper are affine functions in terms of :math:`f` and
:math:`g`, we propose to solve the "one–step nonsmooth problem”
([eq:toto1]) by performing an external Newton linearization.

Newton’s linearization of the first line of ([eq:toto1])
''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The first line of the problem ([eq:toto1]) can be written under the form
of a residue :math:`\mathcal R` depending only on :math:`x_{k+1}` and
:math:`r_{k+1}` such that

.. math:: \label{eq:NL3}  \mathcal R (x_{k+1},r _{k+1}) =0

with

.. math:: \mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h\gamma r- h(1-\gamma)r_k.

The solution of this system of nonlinear equations is sought as a limit
of the sequence
:math:`\{ x^{\alpha}_{k+1},r^{\alpha}_{k+1} \}_{\alpha \in {\mbox{\rm $I\!\!N$}}}`
such that

.. math:: \label{eq:NL7}   \begin{cases}     x^{0}_{k+1} = x_k \\ \\     r^{0}_{k+1} = r_k \\ \\     \mathcal R_L( x^{\alpha+1}_{k+1},r^{\alpha+1}_{k+1}) = \mathcal     R(x^{\alpha}_{k+1},r^{\alpha}_{k+1})  + \left[ \nabla_{x} \mathcal     R(x^{\alpha}_{k+1},r^{\alpha}_{k+1})\right] (x^{\alpha+1}_{k+1}-x^{\alpha}_{k+1} ) +     \left[ \nabla_{r} \mathcal R(x^{\alpha}_{k+1},r^{\alpha}_{k+1})\right] (r^{\alpha+1}_{k+1} - r^{\alpha}_{k+1} ) =0 \end{cases}

What about :math:`r^0_{k+1}` ?

The residu :math:`\mathcal R _{\free}` is also defined (useful for
implementation only):\ 

.. math:: \mathcal R _{\free}(x) \stackrel{\Delta}{=}  M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k),

\ which yields\ 

.. math:: \mathcal R (x,r) = \mathcal R _{\free}(x)   - h\gamma r - h(1-\gamma)r_k.

.. math:: \mathcal R (x^{\alpha}_{k+1},r^{\alpha}_{k+1}) = \fbox{$\mathcal R^{\alpha}_{k+1} \stackrel{\Delta}{=}  \mathcal R_{\free}(x^{\alpha}_{k+1})  - h\gamma r^{\alpha}_{k+1} - h(1-\gamma)r_k$}\label{eq:rfree-1}

.. math:: \mathcal R_{\free}(x^{\alpha}_{k+1},r^{\alpha}_{k+1} )=\fbox{$ \mathcal R _{\free, k+1} ^{\alpha} \stackrel{\Delta}{=}  M(x^{\alpha}_{k+1} - x_{k}) -h\theta f( x^{\alpha}_{k+1} , t_{k+1}) - h(1-\theta)f(x_k,t_k)$}

The special case of Newton’s linearization of ([eq:toto1]) with FirstOrderType2R ([first-DS2])
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now proceed with the time discretization of ([eq:toto1]) with
FirstOrderType2R ([first-DS2]) by a fully implicit scheme :

.. math:: \begin{array}{l}    \label{eq:mlcp2-toto1-DS2}     M x_{k+1} = M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1})     + h(1-\gamma)r(t_k)  \\[2mm]     y_{k+1} =  h(t_{k+1},x_{k+1},\lambda _{k+1}) \\[2mm]     r_{k+1} = g(t_{k+1},\lambda_{k+1})\\[2mm]  \end{array}

Newton’s linearization of the first line of ([eq:mlcp2-toto1-DS2])
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The linearization of the first line of the
problem ([eq:mlcp2-toto1-DS2]) is similar to the previous case so that
([eq:rfree-2]) is still valid.

Newton’s linearization of the second line of ([eq:mlcp2-toto1-DS2])
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The linearization of the second line of the
problem ([eq:mlcp2-toto1-DS2]) is similar to the previous case so that
([eq:NL11y]) is still valid.

Newton’s linearization of the third line of ([eq:mlcp2-toto1-DS2])
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Since
:math:` K^{\alpha}_{k+1} = \nabla_xg(t_{k+1},\lambda ^{\alpha}_{k+1}) = 0 `,
the linearization of the third line of ([eq:mlcp2-toto1-DS2]) reads as

.. math:: \label{eq:mlcp2-rrL}  \begin{array}{l}    \boxed{r^{\alpha+1}_{k+1} = g(t_{k+1},\lambda ^{\alpha}_{k+1})     + B^{\alpha}_{k+1} ( \lambda^{\alpha+1}-  \lambda^{\alpha}_{k+1} )}         \end{array}

Reduction to a linear relation between :math:`x^{\alpha+1}_{k+1}` and\ :math:`\lambda^{\alpha+1}_{k+1}`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Inserting ([eq:mlcp2-rrL]) into ([eq:rfree-11]), we get the following
linear relation between :math:`x^{\alpha+1}_{k+1}`
and\ :math:`\lambda^{\alpha+1}_{k+1}`, we get the linear relation

.. math:: \label{eq:mlcp2-rfree-13}  \begin{array}{l} \boxed{   x^{\alpha+1}_{k+1}\stackrel{\Delta}{=} x^\alpha_p + \left[ h \gamma (W^{\alpha}_{k+1})^{-1}    B^{\alpha}_{k+1} \lambda^{\alpha+1}_{k+1}\right]}   \end{array}

with

.. math:: \boxed{x^\alpha_p \stackrel{\Delta}{=}  h\gamma(W^{\alpha}_{k+1} )^{-1}\left[g(t_{k+1},\lambda^{\alpha}_{k+1})     -B^{\alpha}_{k+1} (\lambda^{\alpha}_{k+1}) \right ] +x^\alpha_{\free}}

and

.. math:: \label{eq:mlcp2-NL9}  \begin{array}{l}    W^{\alpha}_{k+1} \stackrel{\Delta}{=} M-h\theta A^{\alpha}_{k+1}\\  \end{array}

Reduction to a linear relation between :math:`y^{\alpha+1}_{k+1}` and\ :math:`\lambda^{\alpha+1}_{k+1}`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Inserting ([eq:mlcp2-rfree-13]) into ([eq:NL11y]), we get the following
linear relation between :math:`y^{\alpha+1}_{k+1}` and
:math:`\lambda^{\alpha+1}_{k+1}`,

.. math:: \begin{array}{l} y^{\alpha+1}_{k+1} = y_p + \left[ h \gamma C^{\alpha}_{k+1} ( W^{\alpha}_{k+1})^{-1}  B^{\alpha}_{k+1} + D^{\alpha}_{k+1} \right]\lambda^{\alpha+1}_{k+1}   \end{array}

with

.. math:: \boxed{y_p = y^{\alpha}_{k+1} -\mathcal R^{\alpha}_{yk+1} + C^{\alpha}_{k+1}(x_q) -D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1} }

The special case of Newton’s linearization of ([eq:toto1]) with FirstOrderType1R ([first-DS1])
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now proceed with the time discretization of ([eq:toto1]) with
FirstOrderType1R ([first-DS1]) by a fully implicit scheme :

.. math:: \begin{array}{l}    \label{eq:mlcp3-toto1-DS1}     M x_{k+1} = M x_{k} +h\theta f(x_{k+1},t_{k+1})+h(1-\theta) f(x_k,t_k) + h \gamma r(t_{k+1})     + h(1-\gamma)r(t_k)  \\[2mm]     y_{k+1} =  h(t_{k+1},x_{k+1}) \\[2mm]     r_{k+1} = g(t_{k+1}\lambda_{k+1})\\[2mm]  \end{array}

The previous derivation is valid with :math:` D^{\alpha}_{k+1} =0`.

Time–discretization of the linear case ([first-DS3]) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now proceed with the time discretization of ([eq:toto1]) with
FirstOrderLinearR ([first-DS3]) by a fully implicit scheme :

.. math:: \begin{array}{l}    \label{eq:toto1-DS3}     M x^{\alpha+1}_{k+1} = M x_{k} +h\theta A x^{\alpha+1}_{k+1}+h(1-\theta) A x_k + h \gamma r^{\alpha+1}_{k+1}+ h(1-\gamma)r(t_k)  +hb\\[2mm]     y^{\alpha+1}_{k+1} =  C x^{\alpha+1}_{k+1} + D \lambda ^{\alpha+1}_{k+1} +Fz +e\\[2mm]     r^{\alpha+1}_{k+1} = B \lambda ^{\alpha+1}_{k+1} \\[2mm]  \end{array}

.. math:: R_{\free} = M(x^{\alpha}_{k+1} - x_{k}) -h\theta A x^{\alpha}_{k+1} - h(1-\theta) A x_k -hb_{k+1}

\ 

.. math:: R_{\free} = W(x^{\alpha}_{k+1} - x_{k}) -h A x_{k} -hb_{k+1}

Resulting Newton step (only one step)
'''''''''''''''''''''''''''''''''''''

For the sake of simplicity, let us assume that :math:`\gamma =1`

.. math:: \begin{array}{l}     (M -h\theta A)x^{\alpha+1}_{k+1} = M x_{k} +h(1-\theta) A x_k + hr^{\alpha+1}_{k+1} + hb\\[2mm]     y^{\alpha+1}_{k+1} =  C x^{\alpha+1}_{k+1} + D \lambda ^{\alpha+1}_{k+1} +Fz + e \\[2mm]     r^{\alpha+1}_{k+1} = B \lambda ^{\alpha+1}_{k+1}\\[2mm]  \end{array}

that lead to with: :math:` (M -h\theta A) = W`

.. math:: \begin{array}{l}     x^{\alpha+1}_{k+1} = W^{-1}(M x_{k} +h(1-\theta) A x_k + r^{\alpha+1}_{k+1} +hb) = x\free + W^{-1}(r^{\alpha+1}_{k+1})\\[2mm]     y^{\alpha+1}_{k+1} =  ( D+hCW^{-1}B) \lambda ^{\alpha+1}_{k+1} +Fz + CW^{-1}(M     x_k+h(1-\theta)Ax_k + hb) +e \\[2mm]  \end{array}

with
:math:`x_{\free} = x^{\alpha}_{k+1} + W^{-1}(-R_{\free})= x^{\alpha}_{k+1} - W^{-1}(W(x^{\alpha}_{k+1}- x_k) -hAx_k-hb_{k+1} )= W^{-1}(Mx_k +h(1-\theta)Ax_k +h b_{k+1})`

.. math:: \begin{array}{l}     y^{\alpha+1}_{k+1} =  ( D+hCW^{-1}B) \lambda ^{\alpha+1}_{k+1} +Fz + Cx_{\free}+e\\[2mm]     r^{\alpha+1}_{k+1} = B \lambda ^{\alpha+1}_{k+1}\\[2mm]  \end{array}

Coherence with previous formulation
'''''''''''''''''''''''''''''''''''

.. math:: y_p = y^{\alpha}_{k+1} -\mathcal R^{\alpha}_{yk+1} + C^{\alpha}_{k+1}(x_p -x^{\alpha}_{k+1}) -D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1}

\ 

.. math:: y_p = Cx_k + D \lambda _k  + C(\tilde x_{\free}) -D \lambda_k +Fz + e

\ 

.. math:: y_p = Cx_k   + C(\tilde x_{\free})  +Fz + e

\ 

.. math:: y_p = Cx_k   + C(\tilde x_{\free})  +Fz + e

\ 

.. math:: y_p = C(x_{\free})  +Fz + e

Newton’s linearization of  ([eq:toto1-ter]) 
--------------------------------------------

In this section, we deal with only with the FirstOrderType2R case.

.. math::

   \begin{array}{l}
         \label{eq:full-toto1-ter}
         M x_{k+1} = M x_{k} +h \theta f(x_{k+1},t_{k+1}) +h(1-\theta)f(x_{k},t_{k}) + h r_{k+\gamma} \\[2mm]
         y_{k+\gamma} =  h(t_{k+\gamma},x_{k+\gamma},\lambda _{k+\gamma}) \\[2mm]
         r_{k+\gamma} = g(t_{k+\gamma},\lambda_{k+\gamma})\\[2mm]
       \end{array}

Newton’s linearization of the first line of ([eq:full-toto1-ter])
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The first line of the problem ([eq:full-toto1-ter]) can be written under
the form of a residue :math:`\mathcal R` depending only on
:math:`x_{k+1}` and :math:`r_{k+\gamma}` such that

.. math::

   \label{eq:full-NL3}
     \mathcal R (x_{k+1},r _{k+\gamma}) =0

with

.. math:: \mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r.

The solution of this system of nonlinear equations is sought as a limit
of the sequence
:math:`\{ x^{\alpha}_{k+1},r^{\alpha}_{k+\gamma} \}_{\alpha \in {\mbox{\rm $I\!\!N$}}}`
such that

.. math::

   \label{eq:full-NL7}
      \begin{cases}
        x^{0}_{k+1} = x_k \\ \\
        r^{0}_{k+\gamma} = (1-\gamma ) r_{k} + \gamma r^0_{k+1}  = r_k \\ \\     
        \mathcal R_L( x^{\alpha+1}_{k+1},r^{\alpha+1}_{k+\gamma}) = \mathcal
        R(x^{\alpha}_{k+1},r^{\alpha}_{k+\gamma})  + \left[ \nabla_{x} \mathcal
        R(x^{\alpha}_{k+1},r^{\alpha}_{k+\gamma})\right] (x^{\alpha+1}_{k+1}-x^{\alpha}_{k+1} ) + \\[2mm]
        \qquad\qquad\qquad\qquad\qquad\qquad\left[ \nabla_{r} \mathcal R(x^{\alpha}_{k+1},r^{\alpha}_{k+\gamma})\right] (r^{\alpha+1}_{k+\gamma} - r^{\alpha}_{k+\gamma} ) =0
    \end{cases}

What about :math:`r^0_{k+\gamma}` ?

The residu free is also defined (useful for implementation only):

.. math:: \mathcal R _{\free}(x) \stackrel{\Delta}{=}  M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k).

We get

.. math:: \mathcal R (x^{\alpha}_{k+1},r^{\alpha}_{k+\gamma}) = \fbox{$\mathcal R^{\alpha}_{k+1} \stackrel{\Delta}{=}  \mathcal R_{\free}(x^{\alpha}_{k+1} )  - h r^{\alpha}_{k+\gamma}$}\label{eq:full-rfree-1}

.. math::

   \mathcal R
   _{\free}(x^{\alpha}_{k+1} )=\fbox{$ \mathcal R _{\free, k+1} ^{\alpha} \stackrel{\Delta}{=}  M(x^{\alpha}_{k+1} - x_{k}) -h\theta f( x^{\alpha}_{k+1} , t_{k+1}) - h(1-\theta)f(x_k,t_k)$}

The computation of the Jacobian of :math:`\mathcal R` with respect to
:math:`x`, denoted by :math:`   W^{\alpha}_{k+1}` leads to

.. math::

   \label{eq:full-NL9}
      \begin{array}{l}
       W^{\alpha}_{k+1} \stackrel{\Delta}{=} \nabla_{x} \mathcal R (x^{\alpha}_{k+1})= M - h  \theta \nabla_{x} f(  x^{\alpha}_{k+1}, t_{k+1} ).\\
    \end{array}

At each time–step, we have to solve the following linearized problem,

.. math::

   \label{eq:full-NL10}
       \mathcal R^{\alpha}_{k+1} + W^{\alpha}_{k+1} (x^{\alpha+1}_{k+1} -
       x^{\alpha}_{k+1}) - h  (r^{\alpha+1}_{k+\gamma} - r^{\alpha}_{k+\gamma} )  =0 ,

By using ([eq:full-rfree-1]), we get

.. math::

   \label{eq:full-rfree-2}
     \mathcal R _{\free}(x^{\alpha}_{k+1})  - h  r^{\alpha+1}_{k+\gamma}   + W^{\alpha}_{k+1} (x^{\alpha+1}_{k+1} -
       x^{\alpha}_{k+1})  =0

.. math:: \boxed{ x^{\alpha+1}_{k+1} = h(W^{\alpha}_{k+1})^{-1}r^{\alpha+1}_{\gamma+1} +x^\alpha_{\free}}

with :

.. math:: \boxed{x^\alpha_{\free}\stackrel{\Delta}{=}x^{\alpha}_{k+1}-(W^{\alpha}_{k+1})^{-1}\mathcal R_{\free,k+1}^{\alpha} \label{eq:full-rfree-12}}

The matrix :math:`W` is clearly non singular for small :math:`h`.

Note that the linearization is equivalent to the case ([eq:rfree-2]) and
([eq:rfree-12]) with :math:`\gamma=1` and replacing :math:`r_{k+1}` by
:math:`r_{k+\gamma}`.

Newton’s linearization of the second line of ([eq:full-toto1-ter])
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The same operation is performed with the second equation of
([eq:full-toto1-ter])

.. math::

   \begin{array}{l}
       \mathcal R_y(x,y,\lambda)=y-h(t_{k+\gamma},\gamma x + (1-\gamma) x_k ,\lambda) =0\\ \\
     \end{array}

which is linearized as

.. math::

   \label{eq:full-NL9}
     \begin{array}{l}
       \mathcal R_{Ly}(x^{\alpha+1}_{k+1},y^{\alpha+1}_{k+\gamma},\lambda^{\alpha+1}_{k+\gamma}) = \mathcal
       R_{y}(x^{\alpha}_{k+1},y^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma}) +
       (y^{\alpha+1}_{k+\gamma}-y^{\alpha}_{k+\gamma})- \\[2mm] \qquad  \qquad \qquad \qquad  \qquad \qquad
       \gamma C^{\alpha}_{k+1}(x^{\alpha+1}_{k+1}-x^{\alpha}_{k+1}) - D^{\alpha}_{k+\gamma}(\lambda^{\alpha+1}_{k+\gamma}-\lambda^{\alpha}_{k+\gamma})=0
     \end{array}

This leads to the following linear equation

.. math::

   \boxed{y^{\alpha+1}_{k+\gamma} =  y^{\alpha}_{k+\gamma}
     -\mathcal R^{\alpha}_{y,k+1}+ \\
     \gamma C^{\alpha}_{k+1}(x^{\alpha+1}_{k+1}-x^{\alpha}_{k+1}) +
     D^{\alpha}_{k+\gamma}(\lambda^{\alpha+1}_{k+\gamma}-\lambda^{\alpha}_{k+\gamma})}. \label{eq:full-NL11y}

with,

.. math::

   \begin{array}{l}
     C^{\alpha}_{k+\gamma} = \nabla_xh(t_{k+1}, x^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma} ) \\ \\
     D^{\alpha}_{k+\gamma} = \nabla_{\lambda}h(t_{k+1}, x^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma})
    \end{array}

and

.. math::

   \fbox{$
   \mathcal R^{\alpha}_{yk+1} \stackrel{\Delta}{=} y^{\alpha}_{k+\gamma} - h(x^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma})$}

Note that the linearization is equivalent to the case ([eq:NL11y]) by
replacing :math:`\lambda_{k+1}` by :math:`\lambda_{k+\gamma}` and
:math:`x_{k+1}` by :math:`x_{k+\gamma}`.

Newton’s linearization of the third line of ([eq:full-toto1-ter])
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The same operation is performed with the third equation of
([eq:full-toto1-ter])

.. math::

   \begin{array}{l}
       \mathcal R_r(r,\lambda)=r-g(\lambda,t_{k+1}) =0\\ \\  \end{array}

which is linearized as

.. math::

   \label{eq:full-NL9}
     \begin{array}{l}
         \mathcal R_{L\lambda}(r^{\alpha+1}_{k+\gamma},\lambda^{\alpha+1}_{k+\gamma}) = \mathcal
         R_{r,k+\gamma}^{\alpha} + (r^{\alpha+1}_{k+\gamma} - r^{\alpha}_{k+\gamma}) - B^{\alpha}_{k+\gamma}(\lambda^{\alpha+1}_{k+\gamma} -
         \lambda^{\alpha}_{k+\gamma})=0
       \end{array}

.. math::

   \label{eq:full-rrL}
     \begin{array}{l}
       \boxed{r^{\alpha+1}_{k+\gamma} = g(\lambda ^{\alpha}_{k+\gamma},t_{k+\gamma}) -B^{\alpha}_{k+\gamma}
         \lambda^{\alpha}_{k+\gamma} + B^{\alpha}_{k+\gamma} \lambda^{\alpha+1}_{k+\gamma}}       
     \end{array}

with,

.. math::

   \begin{array}{l}
     B^{\alpha}_{k+\gamma} = \nabla_{\lambda}g(\lambda ^{\alpha}_{k+\gamma},t_{k+\gamma})
    \end{array}

and the residue for :math:`r`:

.. math::

   \boxed{\mathcal
         R_{rk+\gamma}^{\alpha} = r^{\alpha}_{k+\gamma} - g(\lambda ^{\alpha}_{k+\gamma},t_{k+\gamma})}

Note that the linearization is equivalent to the case ([eq:rrL]) by
replacing :math:`\lambda_{k+1}` by :math:`\lambda_{k+\gamma}` and
:math:`x_{k+1}` by :math:`x_{k+\gamma}`.

Reduction to a linear relation between :math:`x^{\alpha+1}_{k+1}` and :math:`\lambda^{\alpha+1}_{k+\gamma}`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Inserting ([eq:full-rrL]) into ([eq:full-rfree-12]), we get the
following linear relation between :math:`x^{\alpha+1}_{k+1}` and
:math:`\lambda^{\alpha+1}_{k+1}`,

.. math::

   \begin{array}{l}
        x^{\alpha+1}_{k+1} = h(W^{\alpha}_{k+1} )^{-1}\left[g(\lambda^{\alpha}_{k+\gamma},t_{k+\gamma}) +
       B^{\alpha}_{k+\gamma} (\lambda^{\alpha+1}_{k+\gamma} - \lambda^{\alpha}_{k+\gamma}) \right ] +x^\alpha_{free}
   \end{array}

that is

.. math::

   \begin{array}{l}
   \boxed{x^{\alpha+1}_{k+1}=x_p + h (W^{\alpha}_{k+1})^{-1}    B^{\alpha}_{k+\gamma} \lambda^{\alpha+1}_{k+\gamma}}
      \end{array}
     \label{eq:full-rfree-13}

with

.. math:: \boxed{x_p \stackrel{\Delta}{=}  h(W^{\alpha}_{k+1} )^{-1}\left[g(\lambda^{\alpha}_{k+\gamma},t_{k+\gamma}) -B^{\alpha}_{k+\gamma} (\lambda^{\alpha}_{k+\gamma}) \right ] +x^\alpha_{free}}

Reduction to a linear relation between :math:`y^{\alpha+1}_{k+\gamma}` and :math:`\lambda^{\alpha+1}_{k+\gamma}`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Inserting ([eq:full-rfree-13]) into ([eq:full-NL11y]), we get the
following linear relation between :math:`y^{\alpha+1}_{k+1}` and
:math:`\lambda^{\alpha+1}_{k+1}`,

.. math::

   \begin{array}{l}
    y^{\alpha+1}_{k+1} = y_p + \left[ h \gamma C^{\alpha}_{k+\gamma} ( W^{\alpha}_{k+1})^{-1}  B^{\alpha}_{k+1} + D^{\alpha}_{k+1} \right]\lambda^{\alpha+1}_{k+1}
      \end{array}

with

.. math:: y_p = y^{\alpha}_{k+1} -\mathcal R^{\alpha}_{yk+1} + \gamma C^{\alpha}_{k+1}(x_q) - D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1}

that is

.. math::

   \boxed{
   y_p =  h(x^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma}) + \gamma C^{\alpha}_{k+1}(x_q) - D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1} }

The linear case
'''''''''''''''

.. math::

   \begin{array}{lcl}
       y_p &=&  h(x^{\alpha}_{k+\gamma},\lambda^{\alpha}_{k+\gamma}) + \gamma C^{\alpha}_{k+1}(x_q) - D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1}\\
           &=&  C^{\alpha}_{k+1} x^{\alpha}_{k+\gamma} + D^{\alpha}_{k+1}\lambda^{\alpha}_{k+\gamma}  + \gamma C^{\alpha}_{k+1}(x_q) - D^{\alpha}_{k+1} \lambda^{\alpha}_{k+1} \\
           &=& C^{\alpha}_{k+1}  (x^{\alpha}_{k+\gamma} + \gamma x_p - \gamma x^{\alpha}_{k+1} ) \\
           &=& C^{\alpha}_{k+1}  ((1-\gamma) x_{k} + \gamma x_{free} ) \text {since } x_p =x_{free} 
   \end{array}

Implementation details
''''''''''''''''''''''

For the moment (Feb. 2011), we set
:math:`x_q=(1-\gamma) x_{k} + \gamma x_{free} ` in the linear case. The
nonlinear case is not yet implemented since we need to change the
management of `` Halpha`` Relation to be able to compute the mid–point
values.

Newton’s linearization for Lagrangian systems
=============================================

+-----------+------------------+
| author    | V. Acary         |
+===========+==================+
| date      | Sept, 20, 2011   |
+-----------+------------------+
| version   |                  |
+-----------+------------------+

This section is devoted to the implementation and the study of the
algorithm. The interval of integration is :math:`[0,T]`, :math:`T>0`,
and a grid :math:`t_{0}=0`, :math:`t_{k+1}=t_{k}+h`, :math:`k \geq 0`,
:math:`t_{N}=T` is constructed. The approximation of a function
:math:`f(\cdot)` on :math:`[0,T]` is denoted as :math:`f^{N}(\cdot)`,
and is a piecewise constant function, constant on the intervals
:math:`[t_{k},t_{k+1})`. We denote :math:`f^{N}(t_{k})` as
:math:`f_{k}`. The time-step is :math:`h>0`.

Various second order dynamical systems with input/output relations
------------------------------------------------------------------

Lagrangian dynamical systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class LagrangianDS defines and computes a generic ndof-dimensional
Lagrangian Non Linear Dynamical System of the form :

.. math::

   \begin{cases}
       M(q,z) \dot v + N(v, q, z) + F_{Int}(v , q , t, z) = F_{Ext}(t, z) + p \\
       \dot q = v
     \end{cases}

where

-  :math:`q \in R^{ndof} ` is the set of the generalized coordinates,

-  :math:` \dot q =v \in R^{ndof} ` the velocity, i. e. the time
   derivative of the generalized coordinates (Lagrangian systems).

-  :math:` \ddot q =\dot v \in R^{ndof} ` the acceleration, i. e. the
   second time derivative of the generalized coordinates.

-  :math:` p \in R^{ndof} ` the reaction forces due to the Non Smooth
   Interaction.

-  :math:` M(q) \in R^{ndof \times ndof}
      ` is the inertia term saved in the SiconosMatrix mass.

-  :math:`
      N(\dot q, q) \in R^{ndof}` is the non linear inertia term saved in
   the SiconosVector \_NNL.

-  :math:` F_{Int}(\dot q , q , t) \in
      R^{ndof} ` are the internal forces saved in the SiconosVector
   fInt.

-  :math:` F_{Ext}(t) \in R^{ndof} ` are the external forces saved in
   the SiconosVector fExt.

-  :math:` z \in R^{zSize}` is a vector of arbitrary algebraic
   variables, some sort of discrete state.

The equation of motion is also shortly denoted as:

.. math:: M(q,z) \dot v = F(v, q, t, z) + p

where :math:`F(v, q, t, z) \in R^{ndof} ` collects the total forces
acting on the system, that is

.. math:: F(v, q, t, z) =  F_{Ext}(t, z) -  NNL(v, q, z) + F_{Int}(v, q , t, z)

This vector is stored in the SiconosVector \_Forces

Fully nonlinear case
~~~~~~~~~~~~~~~~~~~~

Let us introduce the following system,

.. math::

   \label{eq:FullyNonLinear}
     \begin{cases}
       M(q,z) \dot v = F(v, q, t, z) + p  \\
       \dot q = v \\
       y = h(t,q,\lambda) \\
       p = g(t,q,\lambda)
     \end{cases}

where :math:`\lambda(t) \in {\mbox{\rm $I\!\!R$}}^m` and
:math:`y(t) \in {\mbox{\rm $I\!\!R$}}^m` are complementary variables
related through a multi-valued mapping. According to the class of
systems, we are studying, the function :math:`F` , :math:`h` and
:math:`g` are defined by a fully nonlinear framework or by affine
functions. This fully nonlinear case is not implemented in Siconos yet.
This fully general case is not yet implemented in Siconos.

Lagrangian Rheonomous relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \label{eq:RheonomousNonLinear}
     \begin{cases}
       M(q,z) \dot v = F(v, q, t, z) + p \\
       \dot q = v \\
       y = h(t,q) \\
       p = G(t,q)\lambda)
     \end{cases}

Lagrangian Scleronomous relations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \label{eq:ScleronomousNonLinear}
     \begin{cases}
       M(q,z) \dot v  = F(v, q, t, z) + p  \\
       \dot q = v \\
       y = h(q) \\
       p = G(q)\lambda
     \end{cases}

Fully Linear case
'''''''''''''''''

.. math::

   \label{eq:FullyLinear}
     \begin{cases}
       M \dot v   +C v + Kq = F_{Ext}(t, z) + p  \\
       \dot q = v \\
       y = C q + e + D\lambda  + F z \\
       p = C^T\lambda
     \end{cases}

Moreau–Jean event-capturing scheme
----------------------------------

In this section, a time-discretization method of the Lagrange dynamical
equation ([eq:11]), consistent with the nonsmooth character of the
solution, is presented. It is assumed in this section, as in the other
sections, that :math:`v^+(\cdot)=\dot{q}^{+}(\cdot)` is a locally
bounded variation function. The equation of motion reads as,

.. math::

   \label{eq:11-b}
     \begin{cases}
       M(q(t)) {dv} +N(q(t),v^{+}(t)) dt+  F_{\mathrm{int}}(t, q(t), v^+(t))\,dt = F_{\mathrm{ext}}(t)\,dt + dr \\ \\
      v^+(t)=\dot{q}^+(t) \\ \\
     q(0)=q_{0} \in {\mathcal C}(0),\;\dot{q}(0^{-})=\dot{q}_{0}
     \end{cases}

We also assume that :math:`F_{\mathrm{int}}(\cdot)` and
:math:`F_{\mathrm{ext}}(\cdot)` are continuous with respect to time.
This assumption is made for the sake of simplicity to avoid the notation
:math:`F^+_{\mathrm{int}}(\cdot)` and :math:`F^+_{\mathrm{ext}}(\cdot)`.
Finally, we will condense the nonlinear inertia terms and the internal
forces to lighten the notation. We obtain

.. math::

   \label{eq:11-c}
     \begin{cases}
       M(q(t)) {dv} + F(t, q(t), v^+(t))\,dt = F_{\mathrm{ext}}(t)\,dt + dr \\ \\
      v^+(t)=\dot{q}^+(t) \\ \\
     q(0)=q_{0} \in {\mathcal C}(0),\;\dot{q}(0^{-})=\dot{q}_{0}
     \end{cases}

The NSCD method, also known as the Contact Dynamics (CD) is due to the
seminal works of J.J.  and M.  (See also ). A lot of improvements and
variants have been proposed over the years. In this Section, we take
liberties with these original works, but we choose to present a version
of the NSCD method which preserves the essential of the original work.
Some extra developments and interpretations are added which are only
under our responsibility. To come back to the source of the NSCD method,
we encourage to read the above references.

The Linear Time-invariant NonSmooth Lagrangian Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the sake of simplicity of the presentation, the linear
time-invariant case is considered first. The nonlinear case will be
examined later in this chapter.

.. math::

   \label{eq:11-a}
     \begin{cases}
       M dv + (K q(t) + C v^+(t))\,dt = F_{\mathrm{ext}}(t)\,dt + dr  \\ \\
       v^+(t)=\dot{q}^+(t)
     \end{cases}

Time–discretization of the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Integrating both sides of this equation over a time step
:math:`(t_k,t_{k+1}]` of length :math:`h>0`, one obtains

.. math::

   \begin{aligned}
     \begin{cases}
       \displaystyle \int_{(t_k,t_{k+1}]} M dv + \int_{t_k}^{t_{k+1}} (C v^+(t)
         + K q(t)) \,dt = \displaystyle \int_{t_k}^{t_{k+1}} F_{\mathrm{ext}}\,dt +
           \displaystyle \int_{(t_k,t_{k+1}]} dr \:, \\ \\
        q(t_{k+1}) = q(t_{k}) + \displaystyle \int_{t_k}^{t_{k+1}} v^+(t)\,dt 
      \end{cases}\end{aligned}

By definition of the differential measure :math:`dv`, we obtain

.. math::

   \begin{aligned}
   \label{eq:19}
   &  \displaystyle \int_{(t_k,t_{k+1}]} M \,dv = M \int_{(t_k,t_{k+1}]}\,dv = M\,(v^+(t_{k+1})-v^+(t_{k}))  &\end{aligned}

Note that the right velocities are involved in this formulation. The
impulsion :math:`\displaystyle \int_{(t_k,t_{k+1}]} dr` of the reaction
on the time interval :math:`(t_k,t_{k+1}]` emerges as a natural unknown.
The equation of the nonsmooth motion can be written under an integral
form as:

.. math::

   \begin{aligned}
     \begin{cases}
        M\,(v(t_{k+1})-v(t_{k})) =   \displaystyle   \int_{t_k}^{t_{k+1}} (- C v^+(t)
         - K q(t) +  F_{\mathrm{ext}}(t))\,dt +
           \displaystyle \int_{(t_k,t_{k+1}]} dr \:, \\ \\
        q(t_{k+1}) = q(t_{k}) + \displaystyle \int_{t_k}^{t_{k+1}} v^+(t)\,dt 
      \end{cases}\end{aligned}

Choosing a numerical method boils down to choose a method of
approximation for the remaining integral terms. Since discontinuities of
the derivative :math:`v(\cdot)` are to be expected if some shocks are
occurring, i.e.. :math:`dr` has some atoms within the interval
:math:`(t_k,t_{k+1}]`, it is not relevant to use high order
approximations integration schemes for :math:`dr` (this was pointed out
in remark [remark1023]). It may be shown on some examples that, on the
contrary, such high order schemes may generate artefact numerical
oscillations (see ).

The following notation will be used:

-  :math:`q_{k}` is an approximation of :math:`q(t_{k})` and
   :math:`q_{k+1}` is an approximation of :math:`q(t_{k+1})`,

-  :math:`v_{k}` is an approximation of :math:`v^+(t_{k})` and
   :math:`v_{k+1}` is an approximation of :math:`v^+(t_{k+1})`,

-  :math:`p_{k+1}` is an approximation of
   :math:` \displaystyle \int_{(t_k,t_{k+1}]} \,dr`.

A popular first order numerical scheme, the so called
:math:`\theta`-method, is used for the term supposed to be sufficiently
smooth:

.. math::

   \begin{aligned}
     \displaystyle \int_{t_k}^{t_{k+1}} C v + K q \,dt  &\approx& 
     h \left[ \theta (C v_{k+1}+K q_{k+1}) + (1-\theta) (C v_{k}+K q_{k}) \right]   \nonumber \\
     \displaystyle \int_{t_k}^{t_{k+1}} F_{\mathrm{ext}}(t) \,dt &\approx& 
     h\left[\theta  (F_{\mathrm{ext}})_{k+1}+(1-\theta)  (F_{\mathrm{ext}})_{k}  \right]  \nonumber \end{aligned}

The displacement, assumed to be absolutely continuous, is approximated
by:

.. math::

   \begin{aligned}
   &  q_{k+1} = q_{k} +  h\,\left[\theta v_{k+1}+(1-\theta) v_{k}  \right] & \nonumber\end{aligned}

Taking into account all these discretizations, the following
time-discretized equation of motion is obtained:

.. math::

   \label{eq:NSCD-discret}
   \begin{cases}
       M (v_{k+1}-v_{k}) + h\left[\theta  (C  v_{k+1}+K q_{k+1}) + (1-\theta) (C v_{k}+K q_{k})  \right] = \\ \\
       \quad\quad\quad\quad\quad = h\left[\theta  (F_{\mathrm{ext}})_{k+1}+(1-\theta)  (F_{\mathrm{ext}})_{k}  \right] + p_{k+1} \\  \\
       q_{k+1} = q_{k} +  h\left[\theta v_{k+1}+(1-\theta) v_{k} \right]
   \end{cases}

Finally, introducing the expression of :math:`q_{k+1}` in the first
equation of ([eq:NSCD-discret]), one obtains:

.. math::

   \begin{aligned}
     \label{eq:23}
   &  \left[M+h\theta C + h^2 \theta^2 K\right] (v_{k+1} -v_{k}) = - h  C v_{k} - h K q_{k} - h^2 \theta  K v_{k} & \nonumber \\ \nonumber \\
   &+  h\left[\theta  (F_{\mathrm{ext}})_{k+1})+(1-\theta)  (F_{\mathrm{ext}})_{k}  \right]  + p_{k+1}  \:, &\end{aligned}

which can be written as:

.. math::

   \begin{aligned}
     \label{eq:24}
      v_{k+1} = v_{\mathrm{free}}  + \widehat{M}^{-1} p_{k+1}\end{aligned}

where,

-  the matrix

   .. math:: \widehat{M} = \left[M+h\theta C + h^2 \theta^2 K \right]  \label{eq:2002}

   is usually called the iteration matrix.

-  The vector

   .. math::

      \label{eq:2003}
      \begin{array}{ll}
      v_{\mathrm{free}}  & = v_{k} + \widehat{M}^{-1} \left[   - h  C v_{k} - h K q_{k} - h^2 \theta  K v_{k} \right. \\ \\ 
      & \left. +  h\left[ \theta  (F_{\mathrm{ext}})_{k+1})+(1-\theta)  (F_{\mathrm{ext}})_{k} \right] \right] 
      \end{array}

   is the so-called “free” velocity, i.e., the velocity of the system
   when reaction forces are null.

Comments
^^^^^^^^

Let us make some comments on the above developments:

-  The iteration matrix
   :math:` \widehat{M} = \left[M+h\theta C + h^2 \theta^2 K \right] ` is
   supposed to be invertible, since the mass matrix :math:`M` is usually
   positive definite and :math:`h` is supposed to be small enough. The
   matrices :math:`C` and :math:`K` are usually semi-definite positive
   since rigid motions are allowed to bodies.

-  When :math:`\theta=0`, the :math:`\theta`-scheme is the explicit
   Euler scheme. When :math:`\theta=1`, the :math:`\theta`-scheme is the
   fully implicit Euler scheme. When dealing with a plain ODE

   .. math:: M\ddot{q}(t)  + C \dot{q}(t) + K q(t)  = F(t)

   the :math:`\theta-`\ scheme is unconditionally stable for
   :math:`0.5 < \theta \leq 1`. It is conditionally stable otherwise.

-  The equation ([eq:24]) is a linear form of the dynamical equation. It
   appears as an affine relation between the two unknowns,
   :math:`v_{k+1}` that is an approximation of the right derivative of
   the Lagrange variable at time :math:`t_{k+1}`, and the impulse
   :math:`p_{k+1}`. Notice that this scheme is fully implicit. Nonsmooth
   laws have to be treated by implicit methods.

-  From a numerical point of view, two major features appear. First, the
   different terms in the numerical algorithm will keep finite values.
   When the time step :math:`h` vanishes, the scheme copes with finite
   jumps. Secondly, the use of differential measures of the time
   interval :math:`(t_k,t_{k+1}]`, i.e..,
   :math:`dv((t_{k},t_{k+1}])=v^+(t_{k+1})-v^+(t_{k})` and
   :math:`dr((t_{k},t_{k+1}])`, offers a rigorous treatment of the
   nonsmooth evolutions. It is to be noticed that approximations of the
   acceleration are ignored.

These remarks on the contact dynamics method might be viewed only as
some numerical tricks. In fact, the mathematical study of the second
order MDI by Moreau provides a sound mathematical ground to this
numerical scheme. It is noteworthy that convergence results have been
proved for such time-stepping schemes , see below.

The Nonlinear NonSmooth Lagrangian Dynamics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time–discretization of the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Starting from the nonlinear dynamics ([eq:11-c]), the integration of
both sides of this equation over a time step :math:`(t_k,t_{k+1}]` of
length :math:`h>0` yields

.. math::

   \begin{aligned}
     \begin{cases}
       \displaystyle \int_{(t_k,t_{k+1}]} M(q) dv + \int_{t_k}^{t_{k+1}} F(t, q(t), v^+(t)) \,dt = \displaystyle \int_{t_k}^{t_{k+1}} F_{\mathrm{ext}}(t)\,dt +
           \displaystyle \int_{(t_k,t_{k+1}]} dr \:, \\ \\
        q(t_{k+1}) = q(t_{k}) + \displaystyle \int_{t_k}^{t_{k+1}} v^+(t)\,dt 
      \end{cases}\end{aligned}

The first term is generally approximated by

.. math::

   \label{eq:19-NL}
     \displaystyle \int_{(t_k,t_{k+1}]} M(q) \,dv \approx  M(q_{k+\gamma})\,(v_{k+1}-v_{k})

where :math:`q_{k+\gamma}` generalizes the standard notation for
:math:`\gamma \in [0,1]` such that

.. math::

   \label{eq:NL1}
     q_{k+\gamma} = (1-\gamma) q_{k} + \gamma\,  q_{k+1}

The *a priori* smooth terms are evaluated with a :math:`\theta`-method,
chosen in this context for its energy conservation ability,

.. math::

   \begin{aligned}
     \displaystyle \int_{t_k}^{t_{k+1}} F(t,q,v) \,dt  &\approx& 
     h  \tilde F_{k+\theta} \end{aligned}

where :math:`\tilde F_{k+\theta}` is an approximation with the following
dependencies

.. math:: \tilde F(t_k,q_k,v_k,t_{k+1},q_{k+1},v_{k+1},t_{k+\theta},q_{k+\theta},v_{k+\theta})

The mid-values :math:`t_{k+\theta},q_{k+\theta},v_{k+\theta}` are
defined by

.. math::

   \label{eq:NSCD-discret-b}
     \left\{\begin{array}{l}
     t_{k+\theta} = \theta t_{k+1}+(1-\theta) t_{k}\\
     q_{k+\theta} = \theta q_{k+1}+(1-\theta) q_{k}\\
     v_{k+\theta} = \theta v_{k+1}+(1-\theta) v_{k}
     \end{array}\right.,\quad  \theta \in [0,1]

[eq:Simo] The choice of the approximated function
:math:`\tilde F(\cdot)` strongly depends on the nature of the internal
forces that are modeled. For the linear elastic behavior of homogeneous
continuum media, this approximation can be made by:

.. math:: \tilde F_{k+\theta} = \frac 1 2 K{{\,:\,}}\left[E(q_{k})+E(q_{k+1})\right] {{\,:\,}}F(q_{k+1/2})

where :math:`E(:cdot)` is the Green-Lagrange strain tensor, which leads
to an energy conserving algorithm as in . For nonlinear elastic other
smooth nonlinear behaviors, we refer to the work of and references
therein for the choice of the discretization and the value of
:math:`\theta`.

The displacement, assumed to be absolutely continuous is approximated
by:

.. math::

   \begin{aligned}
   &  q_{k+1} = q_{k} +  h\,v_{k+\theta}  & \nonumber\end{aligned}

The following nonlinear time–discretized equation of motion is obtained:

.. math::

   \label{eq:NSCD-discret-nl}
   \begin{cases}
       M(q_{k+\gamma}) (v_{k+1}-v_{k}) + h \tilde F_{k+\theta} = p_{k+1} \\  \\
       q_{k+1} = q_{k} +  h v_{k+\theta}
   \end{cases}

In its full generality and at least formally, substituting the
expression of :math:`q_{k+\gamma},q_{k+1}` and :math:`q_{k+\theta}`, the
first line of the problem can be written under the form of a residue
:math:`\mathcal R` depending only on :math:`v_{k+1}` such that

.. math::

   \label{eq:NL3}
     \mathcal R (v_{k+1}) = p_{k+1}

In the last expression, we have omitted the dependence to the known
values at the beginning the time–step, i.e. :math:`q_k` and :math:`v_k`.

Linearizing the Dynamics
^^^^^^^^^^^^^^^^^^^^^^^^

The system of equations ([eq:NL3]) for :math:`v_{k+1}` and
:math:`p_{k+1}` can be linearized yielding a Newton’s procedure for
solving it. This linearization needs the knowledge of the Jacobian
matrix :math:`\nabla \mathcal R (\cdot)` with respect to its argument to
construct the tangent linear model.

Let us consider that the we have to solve the following equations,

.. math::

   \label{eq:NL4}
     \mathcal R (u) = 0

by a Newton’s method where

.. math::

   \label{eq:NL6}
       \mathcal R (u) =   M(q_{k+\gamma} ) (v_{k+1}-v_{k}) + h \tilde F_{k+\theta}

The solution of this system of nonlinear equations is sought as a limit
of the sequence :math:`\{ u^{\tau}_{k+1}\}_{\tau \in \nbN}` such that

.. math::

   \label{eq:NL7}
      \begin{cases}
        u^{0}_{k+1} = v_k \\ \\
        \mathcal R_L( u^{\tau+1}_{k+1}) =  \mathcal R (u^{\tau}_{k+1}) + \nabla \mathcal R (u^{\tau}_{k+1} )(u^{\tau+1}_{k+1}-u^{\tau}_{k+1} ) =0
    \end{cases}

In practice, all the nonlinearities are not treated in the same manner
and the Jacobian matrices for the nonlinear terms involved in the
Newton’s algorithm are only computed in their natural variables. In the
following, we consider some of the most widely used approaches.

The Nonlinear Mass Matrix
'''''''''''''''''''''''''

The derivation of the Jacobian of the first term of
:math:`\mathcal R (\cdot)` implies to compute

.. math::

   \label{eq:NL2000}
      \nabla_u  \left(M(q_{k+\gamma}(u) ) (u-v_{k})\right) \text{ with } q_{k+\gamma}(u) = q_k + \gamma h[(1-\theta) v_k+ \theta u].

One gets

.. math::

   \label{eq:NL8}
     \begin{array}{ll}
       \nabla_u  \left(M(q_{k+\gamma}(u) ) (u-v_{k})\right) &=   M(q_{k+\gamma}(u))  + \left[ \nabla_u M(q_{k+\gamma}(u) ) \right] (u-v_{k}) \\ \\
                                            &=    M(q_{k+\gamma}(u)) + \left[h \gamma\theta \nabla_{q} M(q_{k+\gamma}(u))\right]  (u-v_{k}) 
   \end{array}

The notation :math:`\nabla_{u}M(q_{k+\gamma}(u))(u-v_{k})` is to be
understood as follows:

.. math:: \nabla_{u}M(q_{k+\gamma}(u))(u-v_{k})=\frac{\partial}{\partial u}[M(q_{k+\gamma}(u))(u-v_{k})]

which is denoted as
:math:`\frac{\partial M_{ij}}{\partial q^{l}}(q_{k+\gamma}(u))(u^{l}-v_{k}^{l})`
in tensorial notation. [remarkBABAS]

A very common approximation consists in considering that the mass matrix
evolves slowly with the configuration in a single time–step, that is,
the term :math:`\nabla_{q} M(q_{k+\gamma})` is neglected and one gets,

.. math::

   \label{eq:NL9}
       \nabla_u  (M(q_{k+\gamma}(u) ) (u-v_{k})) \approx  M(q_{k+\gamma}(u) )

The Jacobian matrix :math:`\nabla \mathcal R (\cdot)` is evaluated in
:math:`u^{\tau}_{k+1}` which yields for the equation ([eq:NL9])

.. math::

   \label{eq:NL10}
       \nabla_u  (M(q_{k+\gamma} ) (u^{\tau}_{k+1}-v_{k})) \approx  M(q_k + \gamma h [(1-\theta)v_k+\theta u^{\tau}_{k+1}] ) )

The prediction of the position which plays an important role will be
denoted by

.. math::

   \label{eq:NL555}
     \tilde q^{\tau}_{k+1}= q_k + \gamma h [(1-\theta)v_k+\theta u^{\tau}_{k+1}]

Very often, the matrix :math:`M(q_{k+\gamma})` is only evaluated at the
first Newton’s iteration with :math:`u^{0}_{k+1}= v_k` leading the
approximation for the whole step:

.. math::

   M(q_k + \gamma h [(1-\theta)v_k+\theta u^{\tau}_{k+1}] ) )\approx M(q_k + h \gamma v_k)
   \label{eq:NL11}

Another way to interpret the approximation ([eq:NL11]) is to remark that
this evaluation is just an explicit evaluation of the predictive
position ([eq:NL555]) given by :math:`\theta=0`:

.. math::

   \label{eq:NL5}
     \tilde q_{k+1}= q_k + h \gamma v_k

Using this prediction, the problem ([eq:NSCD-discret-nl]) is written as
follows:

.. math::

   \label{eq:NSCD-discret2}
   \begin{cases}
       M(\tilde q_{k+1}) (v_{k+1}-v_{k}) + h \tilde F_{k+\theta} = p_{k+1} \\  \\
       q_{k+1} = q_{k} +  h v_{k+\theta} \\ \\
       \tilde q_{k+1}= q_k + h \gamma v_k
   \end{cases}

The Nonlinear Term :math:`F(t,q,v)`
'''''''''''''''''''''''''''''''''''

The remaining nonlinear term is linearized providing the Jacobian
matrices of :math:`F(t,q,v)` with respect to :math:`q` and :math:`v`.
This expression depends strongly on the choice of the approximation
:math:`\tilde F_{k+\theta}`. Let us consider a pedagogical example,
which is not necessarily the best as the Remark [eq:Simo] suggests but
which is one of the simplest,

.. math::

   \label{eq:NL13}
     \tilde F_{k+\theta} = (1-\theta) F(t_k,q_k,v_k) + \theta F(t_{k+1},q_{k+1},v_{k+1})

The computation of the Jacobian of
:math:`  \tilde F_{k+\theta}(t,q(u),u)` for

.. math:: q(u) = q_k+h[(1-\theta)v_k+\theta u]

is given for this example by

.. math::

   \label{eq:NL12}
     \begin{array}{ll}
       \nabla_u  \tilde F_{k+\theta}(t,q,u) &= \theta \nabla_u  F(t,q(u),u) \\ \\
       &= \theta \nabla_q F(t_{k+1},q(u)   ,u) \nabla_{u} q(u) + \theta \nabla_{u} F(t,q(u),u)    \\ \\
       &= h \theta^2 \nabla_q F(t, q(u)   ,u) + \theta \nabla_{u} F(t,q(u),u) \\   
     \end{array}

The standard tangent stiffness and damping matrices :math:`K_t` and
:math:`C_t` are defined by

.. math::

   \label{eq:NL14}
     \begin{array}{ll}
     K_t(t,q,u) &= \nabla_q F(t, q   ,u) \\ \\
     C_t(t,q,u) &= \nabla_u F(t, q   ,u) \\
   \end{array}

In this case, the Jacobian of :math:`  \tilde F_{k+\theta}(t,q(u),u)`
may be written as

.. math::

   \label{eq:NL15}
     \begin{array}{ll}
       \nabla_u  \tilde F_{k+\theta}(t,q,u) &=  h \theta^2  K_t(t,q,u) + \theta C_t(t, q   ,u)  \\   
     \end{array}

The complete Newton’s iteration can then be written as

.. math::

   \label{eq:NL16}
      \widehat M^{\tau+1}_{k+1} (u^{\tau+1}_{k+1}-u^{\tau}_{k+1})  =  \mathcal R (u^{\tau}_{k+1}) +p^{\tau+1}_{k+1}

where the iteration matrix is evaluated as

.. math:: \widehat M^{\tau+1}_{k+1} = (M(\tilde q^{\tau}_{k+1}) +  h^2 \theta^2  K_t(t_{k+1},q^{\tau}_{k+1},u^{\tau}_{k+1}) + \theta h C_t(t, q^{\tau}_{k+1}   ,u^{\tau}_{k+1}))\label{eq:NL17}

(compare with ([eq:2002])).

The choice of :math:`\theta=0` leads to an explicit evaluation of the
position and the nonlinear forces terms. This choice can be interesting
if the time–step has to be chosen relatively small due to the presence a
very rapid dynamical process. This can be the case in crashes
applications or in fracture dynamics . In this case, the iteration
matrix reduces to
:math:`\widehat M^{\tau+1}_{k+1} = M(\tilde q^{\tau}_{k+1})` avoiding
the expensive evaluation of the tangent operator at each time–step.

This choice must not be misunderstood. The treatment of the nonsmooth
dynamics continues to be implicit.

Schatzman–Paoli ’scheme and its linearizations
----------------------------------------------

The scheme
~~~~~~~~~~

| M(q\ :sub:`k`)(q\ :sub:`k+1`-2q:sub:`k`\ +q\ :sub:`k-1`) - h\ :sup:`2`
F(v\ :sub:`k+`, q\ :sub:`k+`, t\ :sub:`k+theta`) = p\ :sub:`k+1`,
| v\ :sub:`k+1`\ =,
| y\ :sub:`k+1` = h()
| p\ :sub:`k+1`\ = G() :sub:`k+1`
| 0 y\ :sub:`k+1` :sub:`k+1` 0 .

Should we have

.. math:: v_{k+1}={\displaystyle \frac{q_{k+1}-q_{k-1}}{2h}}

or

.. math:: v_{k+1}={\displaystyle \frac{q_{k+1}-q_{k}}{h}}

 ? This question is particularly important for the initialization and
the proposed :math:`\theta`-scheme

The Newton linearization
~~~~~~~~~~~~~~~~~~~~~~~~

Let us define the residu on :math:`q`

.. math::

   \label{eq:residu}
     \mathcal R(q) =   M(q)(q-2q_{k}+q_{k-1})  + h^2 F( (\theta v(q)+ (1-\theta) v_k),\theta q+ (1-\theta) q_k),  t_{k+\theta})  -  p_{k+1}

with

.. math::

   \label{eq:residu-linq1}
     v(q) = {\displaystyle \frac{q-q_{k-1}}{2h}}

that is

.. math::

   \label{eq:residu-linq2}
     \mathcal R(q) =   M(q)(q-2q_{k}+q_{k-1})  + h^2 F( (\theta {\displaystyle \frac{q-q_{k-1}}{2h}} + (1-\theta) v_k),\theta q+ (1-\theta) q_k),  t_{k+\theta})   -  p_{k+1}

Neglecting :math:`\nabla_q  M(q)` we get

.. math::

   \label{eq:iterationmatrix}
    \nabla_q \mathcal R(q^\nu) =   M(q^\nu) + h^2  \theta K(q^\nu,v^\nu) + {\displaystyle \frac{1}{2}} h  \theta C(q^\nu,v^\nu)

and we have to solve

.. math::

   \label{eq:iterationloop}
    \nabla_q \mathcal R(q^\nu)(q^{\nu+1}-q^\nu) = -  \mathcal R(q^\nu) .

Linear version of the scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| M(q\ :sub:`k+1`-2q:sub:`k`\ +q\ :sub:`k-1`) + h\ :sup:`2` (K
q\ :sub:`k+`\ + C v\ :sub:`k+`) = p\ :sub:`k+1`,
| v\ :sub:`k+1`\ =,
| y\ :sub:`k+1` = h()
| p\ :sub:`k+1`\ = G() :sub:`k+1`
| 0 y\ :sub:`k+1` :sub:`k+1` 0 .

Let us define the residu on :math:`q`

.. math::

   \label{eq:residu-linq}
     \mathcal R(q) =   M(q-2q_{k}+q_{k-1})  + h^2 (K(\theta q+ (1-\theta) q_k))+ C (\theta v(q)+ (1-\theta) v_k))  -  p_{k+1}

with

.. math::

   \label{eq:residu-linq1}
     v(q) = {\displaystyle \frac{q-q_{k-1}}{2h}}

that is

.. math::

   \label{eq:residu-linq2}
     \mathcal R(q) =   M(q-2q_{k}+q_{k-1})  + h^2 (K(\theta q+ (1-\theta) q_k)))+  h^2 C (\theta {\displaystyle \frac{q-q_{k-1}}{2h}}+ (1-\theta) v_k))  -  p_{k+1}

In this linear case, assuming that :math:`q^0=q^\nu = q_k`, we get

.. math::

   \label{eq:residu-linq2}
     \mathcal R(q^\nu) =   M(-q_{k}+q_{k-1})  + h^2 (K q_k)+  h^2 C (\theta {\displaystyle \frac{q_k-q_{k-1}}{2h}}+ (1-\theta) v_k))  -  p_{k+1}

What about mixing OnestepIntegrator in Simulation?
--------------------------------------------------

Let us consider that we have two simple linear Lagrangian Dynamical
systems

.. math::

   \label{eq:FullyLinear1}
     \begin{cases}
       M_1 \dot v_1  = F_{1,Ext}(t) + p_1   \\
       \dot q_1 = v_1 
     \end{cases}

and

.. math::

   \label{eq:FullyLinear1}
     \begin{cases}
       M_2 \dot v_2   = F_{2,Ext}(t) + p_2  \\
       \dot q_2 = v_2 \\
     \end{cases}

These Dynamical systems ([eq:FullyLinear1]) and ([eq:FullyLinear1])
might numerically solved by choosing two different time–stepping
schemes. Let us choose for instance Moreau’s scheme
for([eq:FullyLinear1])

.. math::

   \label{eq:FullyLinear1-TS}
     \begin{cases}
       M_1 (v_{1,k+1}-v_{1,k})  = F_{1,Ext}(t_{k+1}) + p_{1,k+1}   \\
       q_{1,k+1} = q_{k}+ h  v_{1,k+\theta} 
     \end{cases}

and Schatzman–Paoli’s sheme for ([eq:FullyLinear1])

.. math::

   \label{eq:FullyLinear1-TS}
     \begin{cases}
       M_2(q_{2,k+1}-2q_{2,k}+q_{2,k-1})  = F_{2,Ext}(t_{k+1}) + p_{2,k+1}  \\
       v_{2,k+1} = {\displaystyle \frac{q_{2,k+1}-q_{2,k-1}}{2h}} \\
     \end{cases}

Let us consider known that we have a LagrangianLinearTIR between this
two DSs such that

.. math::

   \label{eq:LTIR-2DS}
     \begin{array}{l}
     y = q_1-q_2 \geq 0 \\ \\
     p = \left[
     \begin{array}{c}
       1 \\
       -1
     \end{array}\right] \lambda
   \end{array}

and a complementarity condition

.. math::

   \label{eq:CP}
     0\leq y \perp \lambda \geq 0

Many questions are raised when we want to deal with the discrete
systems:

-  Which rules should we use for the discretization of ([eq:CP]) ?

   .. math::

      \label{eq:CP-TS1}
          \text{ if } \bar y_{k+1}\leq 0, \text{ then }  0\leq \dot y _{k+1} + e \dot y_{k} \perp \hat \lambda_{k+1}\geq 0

   or

   .. math::

      \label{eq:CP-TS2}
          0\leq y _{k+1} + e y_{k-1} \perp \tilde \lambda_{k+1}\geq 0

-  Should we assume that :math:`y_{k+1} = q_{1,k+1}-q_{2,k+1}` and
   :math:`\dot y_{k+1} = v_{1,k+1}-v_{2,k+1}`

-  How can we link :math:`\hat \lambda_{k+1}` and
   :math:`\tilde \lambda_{k+1}` with :math:`p_{1,k+1}` and
   :math:`p_{2,k+1}` ?

The third is the more difficult question and is seems that it is not
reasonable to deal with two DS related by one interaction with different
osi.In practice, this should be avoided in Siconos.

NewtonEuler Dynamical Systems
=============================

+------------+-------------------------------------------------------------------+--------------+
| Author     | O. Bonnefon                                                       | 2010         |
+------------+-------------------------------------------------------------------+--------------+
| Revision   | section [Sec:NE:sub:`m`\ otion] to [Sec:NE:sub:`T`\ D] V. Acary   | 05/09/2011   |
+------------+-------------------------------------------------------------------+--------------+
| Revision   | section [Sec:NE:sub:`m`\ otion] V. Acary                          | 01/06/2016   |
+------------+-------------------------------------------------------------------+--------------+
| Revision   | complete edition V. Acary                                         | 06/01/2017   |
+------------+-------------------------------------------------------------------+--------------+

The equations of motion
-----------------------

In the maximal coordinates framework, the most natural choice for the
kinematic variables and for the formulation of the equations of motion
is the Newton/Euler formalism, where the equation of motion describes
the translational and rotational dynamics of each body using a specific
choice of parameters. For the translational motion, the position of the
center of mass :math:`x_{\cg}\in {\mbox{\rm $I\!\!R$}}^3` and its
velocity :math:`v_{\cg} = \dot x_{\cg} \in {\mbox{\rm $I\!\!R$}}^3` is
usually chosen. For the orientation of the body is usually defined by
the rotation matrix :math:`R` of the body-fixed frame with respect to a
given inertial frame.

For the rotational motion, a common choice is to choose the rotational
velocity :math:`\Omega \in {\mbox{\rm $I\!\!R$}}^3` of the body
expressed in the body–fixed frame. This choice comes from the
formulation of a rigid body motion of a point :math:`X` in the inertial
frame as

.. math::

   \label{eq:1}
     x(t) = \Phi(t,X) = x_{\cg}(t) + R(t) X.

The velocity of this point can be written as

.. math::

   \label{eq:2}
     \dot x(t) = v_{\cg}(t) + \dot R(t) X

Since :math:`R^\top R=I`, we get
:math:`R^\top \dot R + \dot R^\top R =0`. We can conclude that it exists
a matrix :math:`\tilde \Omega := R^\top \dot R ` such that
:math:`\tilde \Omega + \tilde \Omega^\top=0`, i.e. a skew symmetric
matrix. The notation :math:`\tilde \Omega` comes from the fact that
there is a bijection between the skew symmetric matrix in
:math:`{\mbox{\rm $I\!\!R$}}^{3\times3}` and
:math:`{\mbox{\rm $I\!\!R$}}^3` such that

.. math::

   \label{eq:3}
     \tilde \Omega x  = \Omega \times x, \quad \forall x\in {\mbox{\rm $I\!\!R$}}^3.

The rotational velocity is then related to the :math:`R` by :

.. math::

   \label{eq:angularvelocity}
     \widetilde \Omega = R^\top \dot R, \text { or equivalently, } \dot R  = R \widetilde \Omega

Using these coordinates, the equations of motion are given by

.. math::

   \label{eq:motion-NewtonEuler}
     \left\{\begin{array}{rcl}
         m \;\dot v_{\cg}  & = &f(t,x_{\cg}, v_{\cg},  R,  \Omega) \\
         I \dot \Omega + \Omega \times I \Omega &= & M(t,x_{\cg}, v_{\cg}, R, \Omega) \\
         \dot x_{\cg}&=& v_{\cg}\\
         \dot R  &=& R \widetilde \Omega
       \end{array}
     \right.

where :math:`m> 0` is the mass,
:math:`I\in {\mbox{\rm $I\!\!R$}}^{3\times 3}` is the matrix of moments
of inertia around the center of mass and the axis of the body–fixed
frame.

The vectors :math:`f(\cdot)\in {\mbox{\rm $I\!\!R$}}^3` and
:math:`M(\cdot)\in {\mbox{\rm $I\!\!R$}}^3` are the total forces and
torques applied to the body. It is important to outline that the total
applied forces :math:`f(\cdot)` has to be expressed in a consistent
frame w.r.t. to :math:`v_{\cg}`. In our case, it hae to be expressed in
the inertial frame. The same applies for the moment :math:`M` that has
to be expressed in the body-fixed frame. If we consider a moment
:math:`m(\cdot)` expressed in the inertial frame, then is has to be
convected to the body–fixed frame thanks to

.. math::

   \label{eq:convected_moment}
     M (\cdot) =R^\top  m (\cdot)

If we perform the time derivation of :math:`RR^\top =I` rather than
:math:`R^\top R=I`, we get :math:`R \dot R^\top + \dot R R^\top =0`. We
can conclude that it exists a matrix
:math:`\tilde \omega := \dot R R^\top ` such that
:math:`\tilde \omega + \tilde \omega^\top=0`, i.e. a skew symmetric
matrix. Clearly, we have

.. math::

   \label{eq:4}
      \tilde \omega = R \tilde \Omega R^\top

and it can be proved that is equivalent to :math:` \omega =R \Omega`.
The vector :math:`\omega` is the rotational velocity expressed in the
inertial frame. The equation of motion can also be expressed in the
inertial frame as follows

.. math::

   \label{eq:motion-NewtonEuler-inertial}
     \left\{\begin{array}{rcl}
         m \;\dot v_{\cg}  & = &f(t,x_{\cg}, v_{\cg},  R,  R^T \omega) \\
         J(R) \dot \omega + \omega \times J(R) \omega &= & m(t,x_{\cg}, v_{\cg}, R, \omega) \\
         \dot x_{\cg}&=& v_{\cg}\\
         \dot R  &=& \widetilde \omega R
       \end{array}
     \right.

where the matrix :math:`J(R) = R I R^T` is the inertia matrix in the
inertial frame. Defining the angular momentum with respect to the
inertial frame as

.. math::

   \label{eq:1}
     \pi(t) = J(R(t)) \omega(t)

the equation of the angular motion is derived from the balance equation
of the angular momentum

.. math::

   \label{eq:5}
     \dot \pi(t) = m(t,x_{\cg}, v_{\cg}, R, \omega)).

For a given constant (time invariant) :math:`\tilde \Omega`, let us
consider the differential equation

.. math::

   \label{eq:5}
     \begin{cases}
       \dot R(t) = R \tilde \Omega\\
       R(0) = I
     \end{cases}

Let us recall the definition of the matrix exponential,

.. math::

   \label{eq:6}
     \exp(A) = \sum_{k=0}^{\infty} \frac {1}{k!} A^k

A trivial solution of is :math:`R(t) = \exp(t\tilde\Omega) ` since

.. math::

   \label{eq:7}
     \frac {d}{dt}(\exp(At)) = \exp(At) A.

More generally, with the initial condition :math:`R(t_0)= R_0`, we get
the solution

.. math:: R(t) = R_0 \exp((t-t_0)\tilde\Omega)\label{eq:8}

Another interpretation is as follows. From a (incremental) rotation
vector, :math:`\Omega` and its associated matrix :math:`\tilde \Omega`,
we obtain a rotation matrix by the exponentation of
:math:`\tilde \Omega`:

.. math::

   \label{eq:9}
     R = \exp(\tilde\Omega).

Since we note that :math:`\tilde \Omega^ 3 = - \theta^2 \tilde \Omega`
with :math:`\theta = \|\Omega\|`, it is possible to get a closed form of
the matrix exponential of :math:`\tilde \Omega`

.. math::

   \label{eq:10}
     \begin{array}[lcl]{lcl}
       \exp(\tilde \Omega) &=& \sum_{k=0}^{\infty} \frac {1}{k!} (\tilde \Omega)^k \\
                           &=&  I_{3\times 3} + \sum_{k=1}^{\infty} \frac {(-1)^{k-1}}{(2k-1)!}  \theta ^{2k-1} \tilde \Omega + (\sum_{k=0}^{\infty} \frac {(-1)^{k-1}}{(2k)!} \theta)^{2k-2} \tilde \Omega^2\\[2mm]
                           &=&  I_{3\times 3} + \frac{\sin{\theta}} {\theta} \tilde \Omega +  \frac{(\cos{\theta}-1)}{\theta^2}\tilde \Omega^2   
     \end{array}

that is

.. math::

   \label{eq:11}
     R =  I_{3\times 3} + \frac{\sin{\theta}} {\theta} \tilde \Omega +  \frac{(\cos{\theta}-1)}{\theta^2}\tilde \Omega^2

The formula is the Euler–Rodrigues formula that allows to compute the
rotation matrix on closed form.

todo :

-  add the formulation in the inertial frame of the Euler equation with
   :math:`\omega =R \Omega`.

-  Check that is the Euler-Rodrigues formula and not the Olinde
   Rodrigues formula. (division by :math:`\theta`)

In the numerical practice, the choice of the rotation matrix is not
convenient since it introduces redundant parameters. Since :math:`R`
must belong to :math:`SO^+(3)`, we have also to satisfy
:math:`\det(R)=1` and :math:`R^{-1}=R^\top`. In general, we use a
reduced vector of parameters :math:`p\in{\mbox{\rm $I\!\!R$}}^{n_p}`
such :math:`R = \Phi(p)` and :math:`\dot p = \psi(p)\Omega `. We denote
by :math:`q` the vector of coordinates of the position and the
orientation of the body, and by :math:`v` the body twist:

.. math::

   q \coloneqq \begin{bmatrix}
       x_{\cg}\\
       p
     \end{bmatrix},\quad 
     v \coloneqq \begin{bmatrix}
        v_{\cg}\\
        \Omega
      \end{bmatrix}.

The relation between :math:`v` and the time derivative of :math:`q` is

.. math::

   \label{eq:TT}
     \dot q = 
     \begin{bmatrix}
        \dot x_{\cg}\\
        \psi(p) \dot p
      \end{bmatrix}
      = 
      \begin{bmatrix}
        I & 0 \\
        0 & \psi(p)
      \end{bmatrix}
      v
      \coloneqq
      T(q) v

with :math:`T(q) \in {\mbox{\rm $I\!\!R$}}^{{3+n_p}\times 6}`. Note that
the twist :math:`v` is not directly the time derivative of the
coordinate vector as a major difference with Lagrangian systems.

The Newton-Euler equation in compact form may be written as:

.. math::

   \label{eq:Newton-Euler-compact}
   \boxed{ \left \{ 
    \begin{aligned}
     &\dot q=T(q)v, \\
     & M \dot v = F(t, q, v)
    \end{aligned}
    \right.}

where :math:`M\in{\mbox{\rm $I\!\!R$}}^{6\times6}` is the total inertia
matrix

.. math::

   M:= \begin{pmatrix}
       m I_{3\times 3} & 0 \\
       0 & I 
     \end{pmatrix},

and :math:`F(t, q, v)\in {\mbox{\rm $I\!\!R$}}^6` collects all the
forces and torques applied to the body

.. math::

   F(t,q,v):= \begin{pmatrix}
       f(t,x_{\cg},  v_{\cg}, R, \Omega ) \\
       I \Omega \times \Omega + M(t,x_{\cg}, v_{\cg}, R, \Omega )
     \end{pmatrix}.

When a collection of bodies is considered, we will use the same notation
as in ([eq:Newton-Euler-compact]) extending the definition of the
variables :math:`q,v` and the operators :math:`M,F` in a straightforward
way.

Basic elements of Lie groups and Lie algebras theory.
-----------------------------------------------------

Let us recall the definitions of the Lie group Theory taken from and .

Differential equation (evolving) on a manifold :math:`\mathcal M`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :math:`d`-dimensional manifold :math:`\mathcal M` is a
:math:`d`-dimensional smooth surface
:math:` M\subset {\mbox{\rm $I\!\!R$}}^n` for some :math:`n\geq d`.

Let :math:`\mathcal M` be a :math:`d`-dimensional manifold and suppose
that :math:`\rho(t) \in\mathcal M` is a smooth curve such that
:math:`\rho(0) = p`. A tangent vector at :math:`p` is defined as

.. math::

   \label{eq:12}
       a = \left. \frac{d}{dt} (\rho(t)) \right|_{t=0}.

The set of all tangents at :math:`p` is called the tangent space at
:math:`p` and denoted by :math:`T\mathcal M|_p`. It has the structure of
a linear space.

A (tangent) vector field on :math:`\mathcal M` is a smooth function
:math:`F : \mathcal M \rightarrow T\mathcal M` such that
:math:`F (p) \in T\mathcal M|_p` for all :math:`p \in \mathcal M`. The
collection of all vector fields on :math:`\mathcal M` is denoted by
:math:`\mathcal X(\mathcal M)`.

[Differential equation (evolving) on :math:`\mathcal M`] Let :math:`F`
be a tangent vector field on :math:`\mathcal M`. By a differential
equation (evolving) on :math:`\mathcal M` we mean a differential
equation of the form

.. math:: \dot y =F(y), t\geq  0, y(0)\in \mathcal M\label{eq:13}

where :math:`F \in \mathcal X(\mathcal M)`. Whenever convenient, we
allow :math:`F` in  to be a function of time, :math:`F = F(t,y)`. The
flow of :math:`F` is the solution operator
:math:`\Psi_{t,F} : \mathcal M \rightarrow  \mathcal M` such that

.. math:: y(t) = \Psi_{t,F} (y0).\label{eq:14}

Lie algebra and Lie group
~~~~~~~~~~~~~~~~~~~~~~~~~

[commutator] Given two vector fields :math:`F, G` on
:math:`{\mbox{\rm $I\!\!R$}}^n` , the commutator :math:`H = [F, G]` can
be computed componentwise at a given point
:math:`y ∈ {\mbox{\rm $I\!\!R$}}^n` as

.. math:: H_i(y)= \sum_{j=1}^n  G_j(y)\frac{\partial F_i(y)}{\partial y_j}   −F_j(y) \frac{\partial G_i(y)}{\partial y_j} .\label{eq:15}

[lemma:LieBracket] The commutator of vector fields satisfies the
identities

.. math::

   \label{eq:16}
     \begin{array}[lclr]{lclr}
       \protect{[}F, G]&=& −\protect{[}G, F ] & (skew symmetry), \\
       \protect{[} \alpha F,G] &=& \alpha \protect{[}F,G], \alpha \in {\mbox{\rm $I\!\!R$}}&  \\
       \protect{[}F + G, H]&=& \protect{[}F, H] + \protect{[}G, H] & (bilinearity),\\
       0 &=&  \protect{[}F,\protect{[}G,H]]+\protect{[}G,\protect{[}H,F]]+\protect{[}H,\protect{[}F,G]] &(Jacobi’s identity).
     \end{array}

A Lie algebra of vector fields is a collection of vector fields which is
closed under linear combination and commutation. In other words, letting
:math:`\mathfrak g` denote the Lie algebra,

.. math::

   \begin{array}[lclr]{l}
       B \in \mathfrak g \implies \alpha B \in \mathfrak  g \text{ for all } \alpha ∈ R .\\
       B_1,B_2 \in\mathfrak g \implies B_1 +B_2, [B_1,B_2]\in\mathfrak g\label{eq:17}
       \end{array}

Given a collection of vector fields :math:`B = {B_1 , B_2 , \ldots}`,
the least Lie algebra of vector fields containing :math:`B` is called
the Lie algebra generated by :math:`B`

A Lie algebra is a linear space :math:`V` equipped with a Lie bracket, a
bilinear, skew-symmetric mapping

.. math::

   \label{eq:18}
       [ \cdot , \cdot ] : V \times V \rightarrow V

that obeys identities from Lemma [lemma:LieBracket]

[(General) Lie algebra] A Lie algebra homomorphism is a linear map
between two Lie algebras,
:math:`\varphi : \mathfrak g \rightarrow \mathfrak h`, satisfying the
identity

.. math:: \varphi ([v, w]_{\mathfrak g}) = [\varphi(v), \varphi(w)]_{\mathfrak h}, v, w in \mathfrak g\label{eq:19}.

An invertible homomorphism is called an isomorphism.

A Lie group is a differential manifold :math:`\mathcal G` equipped with
a product
:math:`\glaw : \mathcal G\times \mathcal G →\rightarrow \mathcal  G`
satisfying

.. math::

   \label{eq:20}
       \begin{array}[lclr]{lr}
         p \glaw(q \glaw r) = (p\glaw q)\glaw r, \forall  p, q, r ∈ \mathcal G &\text{(associativity)}\\
         \exists I \in \mathcal G \text{ such that } I\glaw p = p \glaw I = p,  \forall p \in \mathcal G&\text{(identity element)}\ \\
         \forall p \in \mathcal G, \exists  p^{-1}  \in \mathcal G \text{ such that }  p^{-1}\glaw p = I&\text{(inverse) }\ \\
         \text{ The maps}  (p, r)  \rightarrow p\glaw r \text{ and }  p  \rightarrow p^{-1} \text{are smooth functions }&\text{(smoothness)}\                                                                                                
       \end{array}

[Lie algebra :math:`\mathfrak g ` of a Lie group :math:`\mathcal G`] The
Lie algebra :math:`\mathfrak g` of a Lie group :math:`\mathcal G` is
defined as the linear space of all tangents to :math:`G` at the identity
:math:`I`. The Lie bracket in :math:`\mathfrak g` is defined as

.. math:: [a,b]= \left.\frac{\partial^2 }{\partial s\partial t} \rho(s)\sigma(t)\rho(-s)\right|_{s=t=0}\label{eq:21}

where :math:`\rho(s)` and :math:`\sigma(t)` are two smooth curves on
:math:`\mathcal G` such that :math:`\rho(0) = \sigma(0) = I`, and
:math:`\dot \rho(0) = a` and :math:`\dot \sigma(0) = b`.

Actions of a group :math:`\mathcal G` on manifold :math:`\mathcal M`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A left action of Lie Group :math:`\mathcal G` on a manifold
:math:`\mathcal M` is a smooth map
:math:`\Lambda^l: \mathcal G \times  \mathcal M \rightarrow \mathcal M`
satisfying

.. math::

   \label{eq:22}
     \begin{array}[lcl]{rcl}
       \Lambda^l(I,y) &=& y, \quad \forall y \in \mathcal M \\
       \Lambda^l(p,\Lambda(r,y)) &=& \Lambda^l(p\glaw r, y) , \quad \forall p,r \in \mathcal G,\quad  \forall y \in \mathcal M .
     \end{array}

A right action of Lie Group :math:`\mathcal G` on a manifold
:math:`\mathcal M` is a smooth map
:math:`\Lambda^r: \mathcal M \times \mathcal G   \rightarrow \mathcal M`
satisfying

.. math::

   \label{eq:23}
     \begin{array}[lcl]{rcl}
       \Lambda^r(y,I) &=& y, \quad \forall y \in \mathcal M \\
       \Lambda^r(\Lambda(y,r), p) &=& \Lambda^r(y,  r\glaw p) , \quad \forall p,r \in \mathcal G,\quad  \forall y \in \mathcal M .
     \end{array}

A given smooth curve
:math:`S(\cdot) : t\in {\mbox{\rm $I\!\!R$}}\mapsto S(t)\in \mathcal G`
in :math:`\mathcal G` such that :math:`S(0)= I` produces a flow
:math:`\Lambda^l(S(t),\cdot)` (resp. :math:`\Lambda^r(\cdot, S(t))`) on
:math:`\mathcal M` and by differentiation we find a tangent vector field

.. math::

   \label{eq:24}
     F(y) = \left. \frac{d}{dt} (\Lambda^l(S(t),y) \right|_{t=0}\quad( \text{resp.  }  F(y) = \left. \frac{d}{dt} (\Lambda^r(y,S(t)) \right|_{t=0} )

that defines a ordinary differential equation on a Lie Group

.. math::

   \label{eq:25}
     \dot y(t) = F(y(t)) = \left. \frac{d}{dt} (\Lambda^l(S(t),y) \right|_{t=0}  \quad( \text{resp.  }\dot y(t) = F(y(t)) = \left. \frac{d}{dt} (\Lambda^r(y,S(t)) \right|_{t=0})

Let
:math:`\lambda^l_{*} : \mathfrak g \rightarrow \mathcal X(\mathcal M) `
(resp.
:math:`\lambda^r_{*} : \mathfrak g \rightarrow \mathcal X(\mathcal M) `
be defined as

.. math:: \lambda^l_{*}(a)(y) = \left.\frac{d}{ds}{ \Lambda^l (\rho(s), y)}\right|_{s=0} \quad (\text{ resp. }  \lambda^r_{*}(a)(y) = \left.\frac{d}{ds}{ \Lambda^r (y, \rho(s))}\right|_{s=0})\label{eq:26}

where :math:`\rho(s)` is a curve in :math:`\mathcal G` such that
:math:`\rho(0)=I` and :math:`\dot\rho (0)=a`. Then :math:`\lambda^l_{8}`
is a linear map between Lie algebras such that

.. math:: [a, b]_{\mathfrak g} = [\lambda^l_{*}(a), \lambda^l_{*}(b)]_{\mathcal X(\mathcal M)}.\label{eq:27}

The following product between an element of an algebra
:math:`a \in \mathfrak g` with an element of a group
:math:`\sigma  \in \mathcal G` can be defined. This will served as a
basis for defining the exponential map.

We define the left product
:math:`(\cdot, \cdot)^l : \mathfrak g \times \mathcal G \rightarrow  \mathcal G`
of an element of an algebra :math:`a \in \mathfrak g` with an element of
a group :math:`\sigma  \in \mathcal G` as

.. math:: (a, \sigma)^l = a \cdot \sigma = \left.\frac{d}{ds} \rho(s) \glaw \sigma \right|_{s=0}\label{eq:28}

where :math:`\rho(s)` is a smooth curve such that :math:`\dot\rho(0)=a`
and :math:`\rho(0)=I`. In the same way, we can define the right product
:math:`(\cdot, \cdot)^r : \mathcal G \times \mathfrak g  \rightarrow   \mathcal G`

.. math::

   \label{eq:29}
     (\sigma,a)^r = \sigma \cdot a  = \left.\frac{d}{ds} \sigma \glaw \rho(s)   \right|_{s=0}

Exponential map
~~~~~~~~~~~~~~~

Let :math:`\mathcal G` be a Lie group and :math:`\mathfrak g` its Lie
algebra. The exponential mapping
:math:`exp : \mathfrak g \rightarrow \mathcal G` is defined as
:math:`\exp(a) = \sigma(1)` where :math:`\sigma (t)` satisfies the
differential equation

.. math:: \dot \sigma(t) = a \cdot \sigma(t), \quad \sigma (0) = I.\label{eq:30}

Let us define :math:`a^k` as

.. math::

   \label{eq:31}
     \left\{\begin{array}[l]{l}
       a^k = \underbrace{a\glaw a \glaw \ldots a\glaw a}_{k \text{ times}} \text{ for } k \geq 1 \\
       a^0  = I
     \end{array}\right.

The exponential map can be expressed as

.. math::

   \label{eq:32}
     \exp(at) = \sum_{k=0}^\infty \frac{(ta)^k}{k!}

since it is a solution of . A simple computation allows to check this
claim:

.. math::

   \label{eq:33}
      \frac{d}{dt}\exp(at) = \sum_{k=1}^\infty  k t^{k-1} \frac{a^k}{k!} = a \glaw \sum_{k=0}^\infty  t^{k} \frac{a^k}{k!} = a \glaw \exp(at).

A similar computation gives

.. math::

   \label{eq:34}
     \frac{d}{dt}\exp(at)  = \sum_{k=0}^\infty  t^{k} \frac{a^k}{k!} \glaw a = \exp(at) \glaw a.

The exponential mapping :math:`exp : \mathfrak g \rightarrow \mathcal G`
can also be defined as :math:`\exp(a) = \sigma(1)` where
:math:`\sigma (t)` satisfies the differential equation

.. math::

   \label{eq:35}
     \dot \sigma(t) = \sigma(t) \cdot a, \quad \sigma (0) = I.

[Theorem:solutionofLieODE] Let
:math:`\Lambda^l:\mathcal G\times\mathcal M \rightarrow \mathcal M` be a
left group action and
:math:`\lambda^l_{∗} : \mathfrak g\rightarrow \mathcal X(\mathcal M)`
the corresponding Lie algebra homomorphism. For any
:math:`a \in \mathfrak g` the flow of the vector field
:math:`F = \lambda^l_{a}(a)`, i.e. the solution of the equation

.. math:: \dot y(t) = F(y(t)) = \lambda^l_{*}(a)(y(t)),\quad  t \geq 0, y(0) = y_0 \in \mathcal M,\label{eq:36}

is given as

.. math:: y(t) = \Lambda^l(\exp(ta), y_0).\label{eq:37}

Let :math:`\Lambda^r:\mathcal M\times\mathcal G \rightarrow \mathcal M`
be a right group action and
:math:`\lambda^r_{∗} : \mathfrak g\rightarrow \mathcal X(\mathcal M)`
the corresponding Lie algebra homomorphism. For any
:math:`a \in \mathfrak g` the flow of the vector field
:math:`F = \lambda^r_{*}(a)`, i.e. the solution of the equation

.. math:: \dot y(t) = F(y(t)) = \lambda^r_{*}(a)(y(t)),\quad  t \geq 0, y(0) = y_0 \in \mathcal M,\label{eq:38}

is given as

.. math:: y(t) = \Lambda^r(y_0,\exp(ta)).\label{eq:39}

Translation (Trivialization) maps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The left and right translation maps defined by

.. math::

   \label{eq:148}
     \begin{array}{rcl}
       L_z  : \mathcal G \times \mathcal G &\rightarrow& \mathcal G \quad \text{ (left translation map )} \\
       y &\mapsto&  z \glaw y
     \end{array}

and

.. math::

   \label{eq:149}
     \begin{array}{rcl}
       R_z(y)  :  \mathcal G \times  \mathcal G  & \rightarrow& \mathcal G \quad \text{ (right translation map )} \\
       y  &\mapsto&  y \glaw z 
     \end{array}

If we identify the manifold :math:`\mathcal M` with the group
:math:`\mathcal G`, The left and right translations can be interpreted
as the simplest example of group action on the manifold. Note that the
left translation map can be viewed as a left or right action on the
group.

If we consider :math:`L_z(y)` as a right group action
:math:` L_z(y) = \Lambda^r( z, y) =z \glaw y `, by differentiation we
get a
:math:`L'_z : T \mathfrak g \cong  \mathfrak g \rightarrow T_z\mathcal G`
with :math:`\dot\rho (0)=a` such that

.. math::

   \label{eq:150}
     \lambda^r_{*}(a)(z) = L'_z(a) = \left.\frac{d}{ds}{ \Lambda^r (z, \rho(s))}\right|_{s=0} = z \glaw a

The map

.. math::

   \label{eq:152}
     \begin{array}{rcl}
     L'_z  : \mathfrak g &\rightarrow& T_z\mathcal G  \\
            a &\mapsto&  z \glaw a
     \end{array}

determines an isomorphism of :math:`\mathfrak g` with the tangent space
:math:`T_z\mathcal G`. In other words, the tangent space can be
identified to :math:`\mathfrak g` as

.. math::

   \label{eq:153}
     T_z\mathcal G =\{L'_z(a) = z \glaw a \mid a \in \mathfrak g  \}

Respectively, if we consider :math:`R_z(y)` as a left group action
:math:` R_z(y) = \Lambda^l( y, z) =y \glaw z `, by differentiation we
get a
:math:`R'_z : T \mathfrak g \cong  \mathfrak g \rightarrow T_z\mathcal G`
with :math:`\dot\rho (0)=a` such that

.. math::

   \label{eq:150}
     \lambda^l_{*}(a)(z) = R'_z(a) = \left.\frac{d}{ds}{ \Lambda^l (\rho(s),z)}\right|_{s=0} = a \glaw z

The map

.. math::

   \label{eq:152}
     \begin{array}{rcl}
     R'_z  : \mathfrak g &\rightarrow& T_z\mathcal G  \\
            a &\mapsto&  a \glaw z
     \end{array}

determines an isomorphism of :math:`\mathfrak g` with the tangent space
:math:`T_z\mathcal G`. In other words, the tangent space can be
identified to :math:`\mathfrak g` as

.. math::

   \label{eq:153}
     T_z\mathcal G =\{R'_z(a) = a \glaw z \mid a \in \mathfrak g  \}

Any tangent vector :math:`F : \mathcal G \rightarrow T_z\mathcal G` can
be written in either of the forms

.. math::

   \label{eq:155}
     F(z) = L'_z(f(a)) = R'_z(g(z))

where :math:`f,g \mathcal G \rightarrow \mathfrak g`.

Adjoint representation
~~~~~~~~~~~~~~~~~~~~~~

Let :math:`p \in \mathcal G` and let :math:`\sigma (t)` be a smooth
curve on :math:`\mathcal G` such that :math:`\sigma (0)` = I and
:math:`\dot \sigma(0) = b \in \mathfrak g`. The adjoint representation
is defined as

.. math:: \operatorname{Ad}_p(b) =\left. \frac{d}{dt} p\sigma(t)p^{-1}\right|_{t=0}\label{eq:40}

The derivative of :math:`\operatorname{Ad}` with respect to the first
argument is denoted :math:`\operatorname{ad}`. Let :math:`\rho(s)` be a
smooth curve on :math:`\mathcal G` such that :math:`\rho(0) = I` and
:math:`\dot \rho(0) = a`, it yields:

.. math::

   \label{eq:41}
       \operatorname{ad}_a(b) = \left.\frac{d}{ds} \operatorname{Ad}_{\rho(s)}(b)\right|_{s=0}  = [a, b]

The adjoint representation can also be expressed with the map

.. math::

   \label{eq:154}
     \operatorname{Ad}_p(b)  = (L_p \glaw R_{p^{-1}})' (b) = (L'_p \glaw R'_{p^{-1}}) (b) = p \glaw b \glaw p^{-1}

For a tangent vector given in , we have

.. math::

   \label{eq:151}
     g(z) = Ad_z(f(z))

Another important relation relating :math:`\operatorname{Ad}`,
:math:`\operatorname{ad}` and :math:`\exp` is

.. math::

   \label{eq:164}
     \operatorname{Ad}_{\exp(a)} =\exp{\operatorname{ad}_a}

Differential of the exponential map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are multiple ways to represent the differential of
:math:`\exp(\cdot)` at a point :math:`a\in \mathfrak g`. Let us start by
the following definition of the differential map at
:math:`a\in\mathfrak g`

.. math::

   \label{eq:147}
     \begin{array}{lcl}
       \exp_a' & : & \mathfrak g \rightarrow  T_{exp(a)}\mathcal G\\
               & &  v \mapsto \exp'_a(v)  = \left.\frac{d}{dt} \exp(a+tv)\right|_{t=0}
     \end{array}

The definition is very similar to the definition of the directional
derivative of :math:`\exp` in the direction :math:`v \in \mathfrak g` at
a point :math:`a\in\mathfrak g`. Using the expression of the tangent
space at :math:`\exp(a)`, we can defined another expression of the
differential map denoted as
:math:`\operatorname{d^l exp}_a : \mathfrak g  \rightarrow \mathfrak g`
such that

.. math::

   \label{eq:156}
     \operatorname{d^l exp}_a = L'_{\exp^{-1}(a)} \glaw \exp_a' = L'_{\exp(-a)} \glaw \exp_a'

This expression appears as a trivialization of the differential map
:math:`\exp'_a`. Using the expression of :math:`L'_z` in . In , an
explicit formula relates :math:`\operatorname{d^l exp}_{a}` to the
iteration of the adjoint operator:

.. math::

   \label{eq:43}
     \operatorname{d^l exp}_a(b) = \sum_{k=0}^\infty \frac{(-1)^k}{(k+1)!} (\operatorname{ad}_a(b))^k \coloneqq \frac{e - \exp\glaw\operatorname{ad}_a}{\operatorname{ad}_a}(b)

where :math:`(\operatorname{ad}_a)^k` is the kth iteration of the
adjoint operator:

.. math::

   \label{eq:44}
     \left\{\begin{array}[l]{l}
       (\operatorname{ad}_a)^k(b) = \underbrace{[a, [ a, [ \ldots, a, [ a, b]]]}_{k \text{ times}} \text{ for } k \geq 1 \\
       (\operatorname{ad}_a)^0(b)  = b
     \end{array}\right.

It is also possible to define the right trivialized differential of the
exponential map

.. math::

   \label{eq:162}
     \operatorname{d^r exp}_a = R'_{\exp^{-1}(a)} \glaw \exp_a' = R'_{\exp(-a)} \glaw \exp_a'

that is

.. math::

   \label{eq:163}
     \operatorname{d^r exp}_a(b) = \exp'_a(b) \glaw \exp(-a)

With these expression, we have equivalently for

.. math::

   \label{eq:157}
      \exp_a'(b)  = \exp_a \glaw \operatorname{d^l exp}_a(b)\quad \text{ and } \exp_a'(b)  = \operatorname{d^r exp}_a(b) \glaw   \exp(a)

To avoid to burden to much the notation, we introduced the unified
definition of the differential map that corresponds to
:math:`\operatorname{dexp}=\operatorname{d^r exp}`

The differential of the exponential mapping, denoted by
:math:`\operatorname{dexp}_a : \mathfrak g \times \mathfrak g \rightarrow \mathfrak g`
is defined as the “right trivialized” tangent of the exponential map

.. math::

   \label{eq:42}
     \frac{d}{dt} (\exp(a(t))) = \operatorname{dexp}_{a(t)}(a'(t)) \exp(a(t))

An explicit formula relates :math:`\operatorname{dexp}_{a}` to the
iteration of the adjoint operator:

.. math::

   \label{eq:43}
     \operatorname{dexp}_a(b) = \sum_{k=0}^\infty \frac{1}{(k+1)!} (\operatorname{ad}_a(b))^k \coloneqq \frac{\exp\glaw\operatorname{ad}_a-e}{\operatorname{ad}_a}(b)

Say what is not the Jacobian in :math:`{\mbox{\rm $I\!\!R$}}^4`

As for :math:`\operatorname{Ad}_a` and :math:`\operatorname{ad}_a`, the
mapping :math:`\operatorname{dexp}_{a}(b)` is a linear mapping in its
second argument for a fixed :math:`a`. Using the relation , we can also
relate the right and the lest trivialization tangent

.. math::

   \label{eq:165}
   \operatorname{d^l exp}_a (b) =   (\operatorname{Ad}_{\exp(a)} \glaw \operatorname{dexp}(a))(b) = (\exp(\operatorname{ad}_{-a}) \glaw \frac{e - \exp\glaw\operatorname{ad}_a}{\operatorname{ad}_a})(b) = \frac{e - \exp\glaw\operatorname{ad}_{-a}}{\operatorname{ad}_a}(b) = \operatorname{dexp}_{-a}(b)

It is also possible to define the the “left trivialized” tangent of the
exponential map

.. math::

   \label{eq:46}
      \frac{d}{dt} (\exp(a(t))) =  \exp(a(t)) \operatorname{d^l exp}_{a(t)}(a'(t)) = \exp(a(t)) \operatorname{dexp}_{-a(t)}(a'(t))

other notation and Lie derivative

.. math::

   \label{eq:178}
         Df \cdot \widehat \Omega (p) = (\widehat \Omega^r f )(p)

Inverse of the exponential map
''''''''''''''''''''''''''''''

The function :math:`\operatorname{dexp}_{a}` is an analytical function
so it possible to invert it to get

.. math::

   \label{eq:45}
     \operatorname{dexp}^{-1}_{a} = \sum_{k=0}^\infty \frac{B_k}{(k)!} (\operatorname{ad}_a)^k(b)

where :math:`B_k` are the Bernouilli number.

Differential of a map :math:`f : \mathcal G \rightarrow \mathfrak g`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We follow the notation developed in . Let us first define the
differential of the map :math:`f : \mathcal G \rightarrow \mathfrak g`
as

.. math::

   \label{eq:166}
     \begin{array}[rcl]{rcl}
       f'_z : T_z\mathcal G &\rightarrow&T_{f(z)}\mathfrak g \cong  \mathfrak g\\
       b &\mapsto& \left.\frac{d}{dt} f(z\glaw \exp(t L'_{z^{-1}}(b))) \right|_{t=0}
     \end{array}

The image of :math:`b` by :math:`f'_z` is obtained by first identifying
:math:`b` with an element of :math:`v \in \mathfrak g` thanks to the
left representation of :math:`T_{f(z)}\mathfrak g` view the left
translation map :math:`v= t L'_z(b)`. The exponential mapping transforms
:math:`v` an element :math:`y` of the Lie Group :math:`\mathcal G`. Then
:math:`f'_z` is obtained by

.. math::

   \label{eq:167}
     f'_z(b) = \lim_{t\rightarrow 0} \frac{f(z\glaw y) - f(z)}{t}

As we have done for the exponential mapping, it is possible to get a
left trivialization of

.. math::

   \label{eq:169}
     \operatorname{d}f_z = (f\glaw L_z)' = f'_z \glaw L'_z

thus

.. math::

   \label{eq:170}
     \operatorname{d}f_z (a) =  f'_z \glaw L'_z(a) = f'_z(L'_z(a)) =  \left.\frac{d}{dt} f(z\glaw \exp(t a )) \right|_{t=0}

Newton Method
'''''''''''''

Let us imagine that we want to solve :math:`f(y) = 0 ` for
:math:`y \in \mathcal G`. A newton method can be written as

 Lie group :math:`SO(3)` of finite rotations and Lie algebra :math:`\mathfrak{so}(3)` of infinitesimal rotations
----------------------------------------------------------------------------------------------------------------

The presentation is this section follows the notation and the
developments taken from . For more details on Lie groups and Lie
algebra, we refer to and .

The Lie group :math:`SO(3)` is the group of linear proper orthogonal
transformations in :math:`{\mbox{\rm $I\!\!R$}}^3` that may be
represented by a set of matrices in
:math:`{\mbox{\rm $I\!\!R$}}^{3\times 3}` as

.. math::

   \label{eq:47}
     SO(3) = \{R \in {\mbox{\rm $I\!\!R$}}^{3\times3}\mid R^TR=I , det(R) = +1  \}

with the group law given by :math:`R_1\glaw R_2 = R_1R_2` for
:math:`R_1,R_2\in SO(3)`. The identity element is
:math:`e = I_{3\times 3}`. At any point of :math:`R\in SO(3)`, the
tangent space :math:`T_RSO(3)` is the set of tangent vectors at a point
:math:`R`.

Left representation of the tangent space at :math:`R`, :math:`T_RSO(3)` 
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Let :math:`S(t)` be a smooth curve
:math:`S(\cdot) : {\mbox{\rm $I\!\!R$}}\rightarrow SO(3)` in
:math:`SO(3)`. An element :math:`a` of the tangent space at :math:`R` is
given by

.. math::

   \label{eq:174}
     a  = \left.\frac{d}{dt} S(t)\right|_{t=0}

such that :math:`S(0)= R`. Since :math:`S(t)\in SO(3)`, we have
:math:`\frac{d}{dt} (S(t)) = \dot S(t)S^T(t) +  S(t) \dot S^T(t) =0`. At
:math:`t=0`, we get :math:`a R^T +  R a^T =0`. We conclude that it
exists a skew–symmetric matrix :math:`\tilde \Omega = R^T a` such that
:math:`\tilde \Omega^T + \tilde \Omega =0`. Hence, a possible
representation of :math:`T_RSO(3)` is

.. math::

   \label{eq:49}
     T_RSO(3) = \{ a = R \tilde \Omega \in {\mbox{\rm $I\!\!R$}}^{3\times 3} \mid \tilde \Omega^T + \tilde \Omega =0 \}.

For :math:`R=I`, the tangent space is directly given by the set of
skew–symmetric matrices:

.. math::

   \label{eq:50}
     T_ISO(3) = \{ \tilde \Omega\in {\mbox{\rm $I\!\!R$}}^{3\times 3} \mid \tilde \Omega^T + \tilde \Omega =0 \}.

The tangent space :math:`T_ISO(3)` with the Lie Bracket
:math:`[\cdot,\cdot]` defined by the matrix commutator

.. math::

   \label{eq:51}
     [A,B] = AB-BA

is a Lie algebra that is denoted by

.. math::

   \label{eq:53}
     \mathfrak{so}(3) =\{\Omega\in {\mbox{\rm $I\!\!R$}}^{3\times 3} \mid \Omega + \tilde \Omega^T =0\}.

For skew symmetric matrices, the commutator can be expressed with the
cross product in :math:`{\mbox{\rm $I\!\!R$}}^3`

.. math::

   \label{eq:52}
     [\tilde \Omega, \tilde \Gamma] = \tilde \Omega \tilde \Gamma - \tilde \Gamma \tilde \Omega= \widetilde{\Omega \times \Gamma }

We use :math:`T_ISO(3) \cong  \mathfrak{so}(3)` whenever there is no
ambiguity.

The notation :math:`\tilde \Omega` is implied by the fact that the Lie
algebra is isomorphic to :math:`{\mbox{\rm $I\!\!R$}}^3` thanks to the
operator
:math:`\widetilde{(\cdot)} :{\mbox{\rm $I\!\!R$}}^3 \rightarrow \mathfrak{so}(3)`
and defined by

.. math::

   \label{eq:54}
    \widetilde{(\cdot)}: \Omega \mapsto \tilde \Omega =
     \begin{bmatrix}
       0 & -\Omega_3 & \Omega_2 \\
       \Omega_3 & 0 & -\Omega_1 \\
       -\Omega_2  & \Omega_1 & 0
     \end{bmatrix}

Note that :math:`\tilde \Omega x = \Omega \times x`.

 A special (right) action of Lie Group :math:`\mathcal G` on a manifold :math:`\mathcal M`. 
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Let us come back to the representation of :math:`T_RSO(3)` given in . It
is clear it can expressed with a representation that relies on
:math:`\mathfrak{so}(3)`

.. math::

   \label{eq:58}
      T_RSO(3) = \{ a = R \tilde \Omega \in {\mbox{\rm $I\!\!R$}}^{3\times 3} \mid \tilde \Omega \in \mathfrak{so}(3) \}.

With , we see that there is a linear map that relates :math:`T_RSO(3)`
to :math:`\mathfrak{so}(3)`. This can be formalize by noting that the
left translation map for a point :math:`R \in SO(3)`

.. math::

   \label{eq:59}
     \begin{array}[lcl]{rcl}
       L_R& :&   SO(3)  \rightarrow  SO(3)\\
          & &  S  \mapsto L_R(S) = R \glaw S = RS\\
     \end{array}

which is diffeomorphism on :math:`SO(3)` is a group action. In our case,
we identify the manifold and the group. Hence, the mapping :math:`L_R`
can be viewed as a left or a right group action. We choose a right
action such that :math:`\Lambda^r(R,S) = L_{R}(S) =  R \glaw S `. By
differentiation, we get a mapping
:math:`L'_R: T_I\mathfrak{so(3)} \cong \mathfrak{so(3)} \rightarrow T_R SO(3)`.
For a given :math:`\tilde\Omega \in \mathfrak{so(3)}` and a point
:math:`R`, the differential :math:`L'_R` by computing the tangent vector
field :math:`\lambda^r_{*}(a)(R)` of the group action
:math:`\Lambda^r(R,S)` for a smooth curve
:math:`S(t) : {\mbox{\rm $I\!\!R$}}\rightarrow S0(3)` such that
:math:`\dot S(0) = \tilde\Omega`:

.. math::

   \label{eq:60}
      \lambda^r_{*}(a)(R) \coloneqq  \left. \frac{d}{dt} \Lambda^r(R,S(t)) \right|_{t=0} = \left. \frac{d}{dt} L_{R}(S(t)) \right|_{t=0} =  \left. \frac{d}{dt} R \glaw S(t) \right|_{t=0} =  R \glaw \dot S(0) = R \tilde\Omega \in X(\mathcal M)

Therefore, the vector field in is a tangent vector field that defines a
Lie-Type ordinary differential equation

.. math::

   \label{eq:61}
     \dot R(t) = \lambda^r_{*}(a)(R(t)) = R(t)  \tilde \Omega

In , the linear operator :math:`\lambda^r_{*}(a)` is defined as the
directional derivative with respect to :math:`S` an denoted
:math:`DL_R(S)`. It defines a diffeomorphism between :math:`T_SSO(3)`
and :math:`T_{RS}SO(3)`. In particular, for :math:`S=I_{3\times3}`, we
get

.. math::

   \label{eq:62}
     \begin{array}{rcl}
       DL_R(I_{3\times3}) : \mathfrak{so}(3) & \rightarrow & T_R SO(3) \\
       \tilde \Omega &\mapsto &DL_R(I_{3\times3}). \tilde \Omega = R \tilde \Omega
     \end{array}

We end up with a possible representation of :math:`T_{R} SO(3)` as

.. math::

   \label{eq:63}
     T_{R} SO(3) =\{\tilde \Omega_R \mid \tilde \Omega_R = DL_R(I_{3\times3}). \tilde \Omega = R \tilde \Omega, \tilde \Omega \in\mathfrak{so}(3)  \}.

In other words, a tangent vector
:math:`\tilde \Omega \in \mathfrak{so}(3)` defines a left invariant
vector field on :math:`SO(3)` at the point :math:`R` given by
:math:`R \tilde \Omega`.

what happens at :math:`S(0)=R`, with
:math:` a =R \tilde \Omega =\dot S(0)` and then
:math:`\dot y(t) = F(y(t)) = R \tilde \Omega y(t) =  R\Omega \times y(t)= \dot S(0) y(t) `.
What else ?

Exponential map :math:`\operatorname{expm}\mathfrak{so(3)} \rightarrow SO(3)`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The relations and shows that is possible to define tangent vector field
from a group action. We can directly apply
Theorem [Theorem:solutionofLieODE] and we get that the solution of

.. math::

   \label{eq:130}
     \begin{cases}
     \dot R(t) = \lambda^r_{*}(a)(R(t)) = R(t)  \tilde \Omega \\
     R(0) = R_0
   \end{cases}

is

.. math::

   \label{eq:138}
     R(t) = R_0 \operatorname{expm}(t \tilde \Omega)

Let us do the computation in this case. Let us assume that the solution
can be sought as :math:`R(t) = \Lambda^r(y_o,S(t))`. The initial
condition imposes that
:math:`R(0) = R_0 = \Lambda(R_0,I) = \Lambda(R_0,S(0))` that implies
:math:`S(0)=I`. Since :math:`\Lambda(R_0,S(t))` is the flow that is
produces by :math:`S(t)` and let us try to find the relation satisfied
by :math:`S(\cdot)`. For a smooth curve :math:`T(s) \in SO(3)` such that
:math:`\dot T(0)= \tilde \Omega`, we have

.. math::

   \label{eq:64}
     \begin{array}[lcl]{lcl}
       \dot R(t) = \lambda^r_*(\tilde \Omega)(R(t)) &=& \left. \frac{d}{ds}\Lambda^r(R(t),T(s)) \right|_{s=0} \\
                                   &=& \left. \frac{d}{ds} \Lambda^r(\Lambda(R_0, S(t)),T(s)) \right|_{s=0} \\
                                   &=& \left. \frac{d}{ds} (\Lambda^r(R_0, S(t)\glaw T(s)) \right|_{s=0} \\
                                   &=& D_2 \Lambda^r(R_0, \glaw S(t) \glaw \dot T(0) ) \\
                                   &=& D_2 \Lambda^r(R_0,  S(t)\glaw \tilde \Omega )
     \end{array}

On the other side, the relation :math:`y(t) = \Lambda^r(y_0,S(t))` gives
:math:`\dot y(t) = D_2 \Lambda^r(y_0,S'(t))` and we conclude that

.. math::

   \label{eq:65}
     \begin{cases}
       \dot S(t) =  S(t)\glaw\tilde\Omega    = S(t) \tilde \Omega\\
       S(0) = I.
     \end{cases}

The ordinary differential equation  is a matrix ODE that admits the
following solution

.. math::

   \label{eq:66}
     S(t) = \operatorname{expm}(t\tilde\Omega)

where
:math:`\exp : {\mbox{\rm $I\!\!R$}}^{3\times 3} \rightarrow {\mbox{\rm $I\!\!R$}}^{3\times 3}`
is the matrix exponential defined by

.. math::

   \label{eq:67}
     \begin{array}[lcl]{lcl}
       \operatorname{expm}(A) &=& \sum_{k=0}^{\infty} \frac {1}{k!} (A)^k.
     \end{array}

We conclude that
:math:`R(t) =\Lambda(R_0,S(t)) = R_0\operatorname{expm}(t\tilde\Omega)`
is the solution of .

We can use the closed form solution for the matrix exponential of
:math:`t \tilde\Omega  \in \mathfrak{so}(3)` as

.. math::

   \label{eq:68}
     \operatorname{expm}(t\tilde\Omega) = I_{3\times 3} + \frac{\sin{t\theta}} {\theta}  \tilde\Omega  +  \frac{(\cos{t \theta}-1)}{\theta^2} \tilde\Omega^2

with :math:`\theta = \|\Omega\|`. For given
:math:`\tilde \Omega \in\mathfrak{so}(3)`, we have

.. math::

   \label{eq:69}
     \det(\tilde \Omega) = \det(\tilde \Omega^T) = \det (-\tilde \Omega^T) = (-1)^3 \det(\tilde \Omega ) = - \det (\tilde \Omega )

that implies that :math:`\det(\tilde \Omega ) =0 `. From , we conclude
that

.. math::

   \label{eq:70}
     \det( \operatorname{expm}(t\tilde \Omega)) = 1.

Furthermore, we have
:math:`\operatorname{expm}(t\tilde \Omega)\operatorname{expm}( -t\tilde \Omega) = \operatorname{expm}(t(\tilde \Omega-\tilde \Omega)) = I`.
We can verify that
:math:`\operatorname{expm}(t\tilde \Omega) \in SO(3)`.

Adjoint representation
''''''''''''''''''''''

In the case of :math:`SO(3)`, the definition of the operator
:math:`\operatorname{Ad}` gives

.. math::

   \label{eq:121}
     \operatorname{Ad}_R(\tilde\Omega)  = R \tilde\Omega R^T

and then mapping :math:`\operatorname{ad}_{\tilde\Omega}(\tilde \Gamma)`
is defined by

.. math::

   \label{eq:56}
     \operatorname{ad}_{\tilde\Omega}(\tilde\Gamma) = \tilde \Omega \tilde\Gamma - \tilde \Gamma \tilde\Omega  =  [\tilde \Omega,\tilde \Gamma] = \widetilde{\Omega \times \Gamma}.

Using the isomorphism between :math:`\mathfrak so(3)` and
:math:`{\mbox{\rm $I\!\!R$}}^3`, it possible the define the mapping
:math:`\operatorname{ad}_{\Omega}(\Gamma) : {\mbox{\rm $I\!\!R$}}^3\times{\mbox{\rm $I\!\!R$}}^3 \rightarrow {\mbox{\rm $I\!\!R$}}^3`
with the realization of the Lie algebra in
:math:`{\mbox{\rm $I\!\!R$}}^3` as

.. math::

   \label{eq:55}
     \operatorname{ad}_\Omega(\Gamma) = \tilde\Omega \Gamma = \Omega\times\Gamma

Differential of the exponential map :math:`\operatorname{dexpm}`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The differential of the exponential mapping, denoted by
:math:`\operatorname{dexpm}` is defined as the ’right trivialized’
tangent of the exponential map

.. math::

   \label{eq:71}
     \frac{d}{dt} (\exp(\tilde \Omega(t))) = \operatorname{dexp}_{\tilde\Omega(t)}(\frac{d \tilde{\Omega}(t)}{dt}) \exp(\tilde\Omega(t))

The differential of the exponential mapping, denoted by
:math:`\operatorname{dexpm}` is defined as the ’left trivialized’
tangent of the exponential map

.. math::

   \label{eq:72}
      \frac{d}{dt} (\exp(\tilde \Omega(t))) = \operatorname{dexp}_{\tilde\Omega(t)}(\frac{d \tilde{\Omega}(t)}{dt}) \exp(\tilde\Omega(t))

Using the formula  and the fact that
:math:`\operatorname{ad}_\Omega(\Gamma) = \Tilde\Omega \Gamma`, we can
write the differential as

.. math::

   \label{eq:122}
     \begin{array}{lcl}
       \operatorname{dexp}_{\tilde\Omega}(\tilde\Gamma) &=& \sum_{k=0}^\infty \frac{1}{(k+1)!} \operatorname{ad}_{\tilde \Omega}^k (\tilde\Gamma) \\
                                          &=& \sum_{k=0}^\infty \frac{1}{(k+1)!} \tilde\Omega^k \tilde \Gamma \\
     \end{array}

Using again the fact that
:math:`\tilde\Omega^3 = -\theta^2 \tilde\Omega`, we get

.. math::

   \label{eq:123}
      \begin{array}{lcl}
        \operatorname{dexp}_{\tilde\Omega} &=& \sum_{k=0}^\infty  \frac{1}{(k+1)!} \tilde\Omega^k \\
                             &=& I  + \sum_{k=0}^\infty  \frac{(-1)^k}{((2(k+1))!} \theta^{2k} \tilde\Omega + \sum_{k=0}^\infty  \frac{(-1)^k}{((2(k+1)+1)!} \theta^{2k} \tilde\Omega^2\\
     \end{array}

Hence, we get

.. math::

   \label{eq:124}
      \begin{array}{lcl}
        \operatorname{dexp}_{\tilde\Omega}  &=& I  + \frac{(1-\cos(\theta))}{\theta^2}\tilde\Omega + \frac{(\theta-\sin(\theta))}{\theta^3}\tilde\Omega^2 
     \end{array}

Since :math:`\operatorname{dexp}_{\tilde\Omega}` is a linear mapping
from :math:`\mathfrak{so(3)}` to :math:`\mathfrak{so(3)}`, we will use
the following notation

.. math::

   \label{eq:172}
     \operatorname{dexp}_{\tilde\Omega}\tilde\Gamma  \coloneqq T(\Omega)\tilde\Gamma

with

.. math::

   \label{eq:173}
      T(\Omega) \coloneqq I  + \frac{(1-\cos(\theta))}{\theta^2}\tilde\Omega + \frac{(\theta-\sin(\theta))}{\theta^3}\tilde\Omega^2  \in {\mbox{\rm $I\!\!R$}}^{3\times 3}

Newton method and differential of a map :math:`f : \mathcal G \rightarrow \mathfrak g`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, let us define the differential of the map
:math:`f : SO(3) \rightarrow \mathfrak {so}(3)` as

.. math::

   \label{eq:183}
     \begin{array}[rcl]{rcl}
       f'_R : T_RSO(3) &\rightarrow&T_{f(R)}\mathfrak {so}(3) \cong  \mathfrak {so}(3)\\
              a &\mapsto& \left.\frac{d}{dt} f(R\glaw \operatorname{expm}(t L'_{R^{-1}}(a))) \right|_{t=0}
     \end{array}

The image of :math:`b` by :math:`f'_z` is obtained by first identifying
:math:`a` with an element of :math:`\tilde\Omega \in \mathfrak {so}(3)`
thanks to the left representation of :math:`T_{f(R)}\mathfrak {so}(3)`
view the left translation map :math:`\tilde\Omega= t L'_R(b)`. The
exponential mapping transforms :math:`\tilde\Omega` an element :math:`S`
of the Lie Group :math:`SO(3)`. Then :math:`f'_z` is obtained by

.. math::

   \label{eq:184}
     f'_R(b) = \lim_{t\rightarrow 0} \frac{f(R\glaw S) - f(R)}{t}

As we have done for the exponential mapping, it is possible to get a
left trivialization of

.. math::

   \label{eq:185}
     \operatorname{d}f_R = (f\circ L_R)' = f'_R \circ L'_R

thus

.. math::

   \label{eq:186}
     \operatorname{d}f_R (\tilde\Omega) =  f'_R \circ L'_R(\tilde\Omega) = f'_R(L'_R(\tilde\Omega)) =  \left.\frac{d}{dt} f(R\glaw \operatorname{expm}(t \tilde\Omega )) \right|_{t=0}

The computation of this differential is non linear with respect to
:math:`\tilde\Omega`.

not clear if we write :math:`\operatorname{d}f_R (\tilde\Omega)`. Better
understand the link with
:math:`\operatorname{dexp}_{\tilde \Omega}{\tilde\Gamma}`

Sometimes, it can be formally written as

.. math::

   \label{eq:180}
     \operatorname{d}f_R (\tilde\Omega) = C(\tilde\Omega)\tilde\Omega

Nevertheless, an explicit expression of :math:`C(\cdot)` is not
necessarily trivial.

Let us consider a first simple example of a mapping
:math:`f(R) = \widetilde{R  x}` for a given
:math:`x\in{\mbox{\rm $I\!\!R$}}^3`. The computation yields

.. math::

   \label{eq:181}
     \begin{array}{rcl}
       \operatorname{d}f_R (\tilde\Omega) &=& \widetilde{ \left.\frac{d}{dt} R \exp(t \tilde\Omega) x  \right|_{t=0}} \\
                              &=& \widetilde{R \left.\frac{d}{dt}\exp(t \tilde\Omega)\right|_{t=0}  x} \\
                              &=& \widetilde{R \left. \operatorname{dexp}_{\tilde\Omega}(\tilde\Omega)\exp(t \tilde\Omega) \right|_{t=0}  x} \\
                              &=& \widetilde{R \operatorname{dexp}_{\tilde\Omega}(\tilde\Omega) x} \\
                              &=& \widetilde{R T(\Omega) \tilde\Omega  x} \\
                              &=& \widetilde{-R T(\Omega) \tilde x \Omega } 
     \end{array}

In that case, it is difficult to find a expression as in , but
considering the function :math:`g(R)` such that
:math:`f(R) = \widetilde g(x)` we get

.. math::

   \label{eq:181}
     \begin{array}{rcl}
       \operatorname{d}g_R (\tilde\Omega)  =- R T(\Omega) \tilde x \Omega  = C(\Omega) \Omega
     \end{array}

with

.. math::

   \label{eq:182}
      C(\Omega) = -R T(\Omega) \tilde x

Lie group of unit quaternions :math:`{\mbox{\rm $I\!\!H$}}_1` and pure imaginary quaternions :math:`{\mbox{\rm $I\!\!H$}}_p`.
-----------------------------------------------------------------------------------------------------------------------------

In Siconos we choose to parametrize the rotation with a unit quaternion
:math:`p \in {\mbox{\rm $I\!\!H$}}` such that :math:`R = \Phi(p)`. This
parameterization has no singularity and has only one redundant variable
that is determined by imposing :math:`\|p\|=1`.

Quaternion definition.
''''''''''''''''''''''

There is many ways to define quaternions. The most convenient one is to
define a quaternion as a :math:`2\times 2` complex matrix, that is an
element of
:math:`{\mbox{\rm $~\vrule height6.6pt width0.5pt depth0.25pt\!\!$C}}^{2\times 2}`.
For this end, we write for
:math:`z \in {\mbox{\rm $~\vrule height6.6pt width0.5pt depth0.25pt\!\!$C}}`,
:math:`z=a+ib` with :math:`a,b \in {\mbox{\rm $I\!\!R$}}^2` and
:math:`i^2=-1` and its conjugate :math:`\bar z= a-ib`. Let
:math:`{e, \bf, i, j, k}` the following matrices in
:math:`{\mbox{\rm $~\vrule height6.6pt width0.5pt depth0.25pt\!\!$C}}^{2\times 2}`

.. math::

   \label{eq:127}
     e =
     \begin{bmatrix}
       1 & 0 \\
       0 & 1  \\
     \end{bmatrix},
     \quad   \bf{i} =
     \begin{bmatrix}
       i & 0 \\
       0 & -i  \\
     \end{bmatrix}
     \quad   \bf{j} =
     \begin{bmatrix}
       0 & 1 \\
       -1 & 0  \\
     \end{bmatrix}
      \quad   \bf{k} =
     \begin{bmatrix}
       0 & i \\
       i & 0  \\
     \end{bmatrix}

Let :math:`{\mbox{\rm $I\!\!H$}}` be the set of all matrices of the form

.. math::

   \label{eq:128}
       p_0 e + p_1 {\bf i} + p_2 {\bf j} + p_3 {\bf k}

where :math:`(p_0,p_1,p_2,p_3) \in {\mbox{\rm $I\!\!R$}}^4`. Every
Matrix in :math:`{\mbox{\rm $I\!\!H$}}` is of the form

.. math::

   \label{eq:129}
       \begin{bmatrix}
         x &y  \\
         - \bar y  & \bar x
       \end{bmatrix}

where :math:`x = p_0 + i p_1` and :math:`y = p_2 + i p_3`. The matrices
in :math:`{\mbox{\rm $I\!\!H$}}` are called quaternions.

The null quaternion generated by
:math:`[0,0,0,0] \in {\mbox{\rm $I\!\!R$}}^4` is denoted by :math:`0` .
Quaternions of the form :math:`p_1 \bf {i} + p_2 \bf{j} + p_3 \bf{k}`
are called pure quaternions. The set of pure quaternions is denoted by
:math:`{\mbox{\rm $I\!\!H$}}_p`.

With the definition of :math:`{\mbox{\rm $I\!\!H$}}` as a set of complex
matrices, It can be show that :math:`{\mbox{\rm $I\!\!H$}}` is a real
vector space of dimension :math:`4` with basis
:math:`{e, \bf, i, j, k}`. Furthermore, with the matrix product,
:math:`{\mbox{\rm $I\!\!H$}}` is a real algebra.

Representation of quaternions
'''''''''''''''''''''''''''''

Thanks to the equations ,   and  , we see that there are several manner
to represent a quaternion :math:`p\in {\mbox{\rm $I\!\!H$}}`. It can be
represented as a complex matrix as in . It can also be represented as a
vector in :math:`{\mbox{\rm $I\!\!R$}}^4` , :math:`p= [p_0,p_1,p_2,p_3]`
with the isomorphism . In other words, :math:`{\mbox{\rm $I\!\!H$}}` is
isomorphic to :math:`{\mbox{\rm $I\!\!R$}}^4`. The first element
:math:`p_0` can also be viewed as a scalar and three last ones as a
vector in :math:`{\mbox{\rm $I\!\!R$}}^3` denoted by
:math:`\vv{p} = [p_1,p_2,p_3]`, and in that case,
:math:`{\mbox{\rm $I\!\!H$}}` is viewed as
:math:`{\mbox{\rm $I\!\!R$}}\times {\mbox{\rm $I\!\!R$}}^3`. The
quaternion can be written as :math:`p=(p_0,\vv{p})`.

Quaternion product
''''''''''''''''''

The quaternion product denoted by :math:`p \glaw q ` for
:math:`p,q\in {\mbox{\rm $I\!\!H$}}_1` is naturally defined as the
product of complex matrices. With its representation in
:math:`{\mbox{\rm $I\!\!R$}}\times {\mbox{\rm $I\!\!R$}}^3`, the
quaternion product is defined by

.. math::

   \label{eq:73}
     p \glaw q =
     \begin{bmatrix}
       p_oq_o - \vv{p}\vv{q} \\
       p_0\vv{q}+q_o\vv{p} + \vv{p}\times\vv{q}
     \end{bmatrix}.

Since the product is a matrix product, it is not communicative, but it
is associative. The identity element for the quaternion product is

.. math::

   e=  \begin{bmatrix}
       1 & 0 \\
       0 & 1  \\
     \end{bmatrix} =(1,\vv{0})\label{eq:57}.

Let us note that

.. math::

   \label{eq:74}
     (0,\vv{p})\glaw (0,\vv(q)) = - (0,\vv{q})\glaw (0,\vv{p}).

The quaternion multiplication can also be represented as a matrix
operation in :math:`{\mbox{\rm $I\!\!R$}}^{4\times4}`. Indeed, we have

.. math::

   \label{eq:75}
     p \glaw q  =
     \begin{bmatrix}
       q_0 p_0 -q_1p_1-q_2p_2-q_3p_3\\
       q_0 p_1 +q_1p_0-q_2p_3+q_3p_2\\
       q_0 p_2 +q_1p_3+q_2p_0-q_3p_1\\
       q_0 p_3 -q_1p_2+q_2p_1+q_3p_0\\
     \end{bmatrix}

that can be represented as

.. math::

   \label{eq:76}
     p \glaw q  =
     \begin{bmatrix}
       p_0 & -p_1 & -p_2 & -p_3 \\
       p_1 & p_0 & -p_3 & p_2 \\
       p_2 & p_3 & p_0 & -p_1 \\
       p_3 & -p_2 & p_1 & p_0 \\
     \end{bmatrix}
     \begin{bmatrix}
       q_0\\
       q_1\\
       q_2\\
       q_3
     \end{bmatrix} := [p_\glaw]q

or

.. math::

   \label{eq:77}
     p \glaw q  = 
     \begin{bmatrix}
       q_0 & -q_1 & -q_2 & -q_3 \\
       q_1 & q_0 & q_3 & -q_2 \\
       q_2 & -q_3 & q_0 & q_1 \\
       q_3 & q_2 & -q_1 & q_0 \\
     \end{bmatrix}
     \begin{bmatrix}
       p_0\\
       p_1\\
       p_2\\
       p_3
     \end{bmatrix} := [{}_\glaw q] p

Adjoint quaternion, inverse and norm
''''''''''''''''''''''''''''''''''''

The adjoint quaternion of :math:`p` is denoted by

.. math::

   p^\star = \overline{  \begin{bmatrix}
         x &y  \\
         - \bar y  & \bar x
       \end{bmatrix}}^T
     =\begin{bmatrix}
       \bar x & - y  \\
       \bar y  &  x
     \end{bmatrix} =
     \begin{bmatrix}
       p_0, -p_1, -p_2, -p_3
     \end{bmatrix} = (p_0, - \vv{p})

We note that

.. math::

   \label{eq:131}
     p \glaw p^\star = 
     \begin{bmatrix}
         x &y  \\
         - \bar y  & \bar x
       \end{bmatrix}
       \begin{bmatrix}
       \bar x & - y  \\
       \bar y  &  x
     \end{bmatrix} = \det(\begin{bmatrix}
         x &y  \\
         - \bar y  & \bar x
       \end{bmatrix}) e = (x\bar x + y \bar y) e  = (p^2_0 + p^2_1+ p_2^2 + p_3^2)e

The norm of a quaternion is given by
:math:`|p|^2=p^\top p = p_o^2+p_1^2+p_2^2+p_3^2`. In particular, we have
:math:`p \glaw p^\star = p^\star \glaw p = |p|^2 e`. This allows to
define the reciprocal of a non zero quaternion by

.. math::

   \label{eq:78}
     p ^{-1} = \frac 1 {|p|^2} p^\star

A quaternion :math:`p` is said to be unit if :math:`|p| =1`.

Unit quaternion and rotation
''''''''''''''''''''''''''''

For two vectors :math:`x\in {\mbox{\rm $I\!\!R$}}^3` and
:math:`x'\in {\mbox{\rm $I\!\!R$}}^3`, we define the quaternion
:math:`p_x = (0,x)\in {\mbox{\rm $I\!\!H$}}_p` and
:math:`p_{x'} = (0,x')\in {\mbox{\rm $I\!\!H$}}_p`. For a given unit
quaternion :math:`p`, the transformation

.. math::

   \label{eq:79}
     p_{x'} = p \glaw p_x \glaw  p^\star

defines a rotation :math:`R` such that :math:`x'  = R x` given by

.. math::

   \label{eq:80}
     x' = (p_0^2- p^\top \vv{p}) x +2 p_0(\vv{p}\times x) +  2 (\vv{p}^\top x) p = R x

The rotation matrix may be computed as

.. math::

   \label{eq:81}
     R = \Phi(p) =
     \begin{bmatrix}
       1-2 p_2^2- 2 p_3^2 & 2(p_1p_2-p_3p_0) & 2(p_1p_3+p_2p_0)\\
       2(p_1p_2+p_3p_0) & 1-2 p_1^2- 2 p_3^2 & 2(p_2p_3-p_1p_0)\\
       2(p_1p_3-p_2p_0) & 2(p_2p_3+p_1p_0)  & 1-2 p_1^2- 2 p_2^2\\
     \end{bmatrix}

Computation of the time derivative of a unit quaternion associated with a rotation.
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The derivation with respect to time can obtained as follows. The
rotation transformation for a unit quaternion is given by

.. math::

   \label{eq:82}
     p_{x'}(t) = p(t) \glaw p_x \glaw p^\star(t) =  p(t) \glaw p_x \glaw p^{-1}(t)

and can be derived as

.. math::

   \label{eq:83}
     \begin{array}{lcl}
       \dot p_{x'}(t) &=& \dot p(t) \glaw p_x \glaw p^{-1}(t) + p(t) \glaw p_x \glaw \dot p^{-1}(t) \\
                     &=& \dot p(t) \glaw p^{-1}(t)  \glaw   p_{x'}(t)  +      p_{x'}(t) \glaw p(t)  \glaw \dot p^{-1}(t)    
     \end{array}

From :math:`p(t) \glaw p^{-1}(t) =e`, we get

.. math::

   \label{eq:84}
     \dot p(t) \glaw p^{-1}(t) + p \glaw \dot p^{-1}(t) = 0

so ([eq:82]) can be rewritten

.. math::

   \label{eq:85}
     \begin{array}{lcl}
       \dot p_{x'}(t) = \dot p(t) \glaw p^{-1}(t)   \glaw   p_{x'}(t)  -    p_{x'}(t) \glaw  \dot p(t) \glaw p^{-1}(t)
     \end{array}

The scalar part of :math:`\dot p(t) \glaw p^{-1}(t)` is
:math:`(\dot p(t) \glaw p^{-1}(t))_0 = p_o \dot p_0 + \vv{p}^T\vv{\dot p}`.
Since :math:`p` is a unit quaternion, we have

.. math::

   \label{eq:86}
     |p|=1 \implies \frac{d}{dt} (p^\top p) = 0 =  \dot p^\top p + p^\top \dot p =   2( p_o \dot p_0 + \vv{p}^T\vv{\dot p}).

Therefore, the scalar part :math:`(\dot p(t) \glaw p^{-1}(t))_0 =0`. The
quaternion product :math:`\dot p(t) \glaw p^{-1}(t)` and
:math:`p_{x'}(t)` is a product of quaternions with zero scalar part
(see ), so we have

.. math::

   \label{eq:87}
     \begin{array}{lcl}
       \dot p_{x'}(t) = 2 \dot p(t) \glaw p^{-1}(t)   \glaw p_{x'}(t).
     \end{array}

In terms of vector of :math:`{\mbox{\rm $I\!\!R$}}^3`, this corresponds
to

.. math::

   \label{eq:88}
     \dot x'(t) = 2 \vv{ \dot p(t) \glaw p^{-1}(t) } \times x'(t).

Since :math:`x'(t) = R(t) x`, we have
:math:`\dot x' = \dot R(t) x = \tilde \omega(t) R(t) x  = \tilde \omega(t) x'(t) `.
Comparing and , we get

.. math::

   \label{eq:89}
     \tilde \omega(t)  = 2 \vv{ \dot p(t) \glaw p^{-1}(t) }

or equivalently

.. math::

   \dot p(t) \glaw p^{-1}(t) = (0, \frac{\omega(t)}{2} )
     \label{eq:90}

Finally, we can conclude that

.. math::

   \label{eq:91}
     \dot p(t) = (0, \frac{\omega(t)}2 ) \glaw p(t).

Since :math:`\omega(t)=R(t)\Omega(t)`, we have

.. math::

   \label{eq:92}
     (0, \omega(t) ) = (0, R(t) \Omega(t) ) = p(t) \glaw (0, \Omega(t) ) \glaw \bar p(t) = p(t) \glaw (0, \Omega(t) ) \glaw  p^{-1}(t)

and then

.. math::

   \label{eq:93}
     \dot p(t) =\frac 1 2 p(t) \glaw(0, \Omega(t) ) .

The time derivation is compactly written

.. math::

   \label{eq:94}
     \dot p = \frac  1 2 p  \glaw(0, \frac\Omega 2 ) =  [p_\glaw] p_{\frac \Omega 2} = \Psi(p)\frac \Omega 2,

and using the matrix representation of product of quaternion we get

.. math::

   \label{eq:95}
     \Psi(p) =  \begin{bmatrix}
       -p_1 & -p_2 & -p_3 \\
       p_0 & -p_3 & p_2 \\
       p_3 & p_0 & -p_1 \\
       -p_2 & p_1 & p_0 \\
     \end{bmatrix}

The relation can be also inverted by writing

.. math::

   \label{eq:96}
      (0, \Omega(t) ) = 2 p^{-1}(t) \glaw \dot p(t)

Using again matrix representation of product of quaternion, we get

.. math::

   \label{eq:97}
     \Omega(t)  = 2 \vv{p^{-1}(t) \glaw \dot p(t)}  = 2  \begin{bmatrix}
       -p_1 & p_0 & p_3 & -p_2 \\
       -p_2 & -p_3 & p_0 & p_1 \\
       -p_3 & p_2 & -p_1  & p_0\\
     \end{bmatrix}\dot p(t) = 2 \Psi(p)^\top \dot p(t)

Note that we have :math:`\Psi^\top(p)\Psi(p)= I_{4\times 4 }` and
:math:`\Psi(p)\Psi^\top(p)= I_{3\times 3 }`

Lie group structure of unit quaternions.
''''''''''''''''''''''''''''''''''''''''

In terms of complex matrices, an unit quaternion :math:`p` satisfies

.. math::

   \label{eq:125}
     \det\left(    \begin{bmatrix}
         x &y  \\
         - \bar y  & \bar x
       \end{bmatrix} \right) =1

The set of all unit quaternions that we denote
:math:`{\mbox{\rm $I\!\!H$}}_1` is the set of unitary matrices of
determinant equal to :math:`1`. From , we get that

.. math::

   \label{eq:126}
     p \glaw p^\star = e

It implies that the set :math:`{\mbox{\rm $I\!\!H$}}_1` is the set of
special unitary complex matrices. The set is a Lie group usually denoted
as :math:`SU(2)`. Since we used multiple representation of a quaternion,
we continue to use :math:`{\mbox{\rm $I\!\!H$}}_1 \cong SU(2)` as a
notation but with the Lie group structure implied by :math:`SU(2)`.

Let us compute the tangent vector at a point
:math:`p \in {\mbox{\rm $I\!\!H$}}_1`. Let :math:`q(t)` be a smooth
curve :math:`q(\cdot) : t\in {\mbox{\rm $I\!\!R$}}\mapsto q(t)\in H_1`
in :math:`H_1` such that :math:`q(O)= p`. Since :math:`q(t)\in H_1`, we
have :math:`|q(t)|=1` and then
:math:`\frac{d}{dt} |q(t)| = 2(q_0(0) \dot q_0(0) + \vec{q}^T(0) \vec{\dot q}(0) ) =0`.
At :math:`t=0`, we get

.. math::

   \label{eq:48}
     2(p_0 a_0 + \vec{p}^T \vec{a})= 0.

This relation imposes that the quaternions
:math:`2 p^\star \glaw a \in H_1` and :math:`2 a \glaw p^\star \in H_1`,
that is, have to be pure quaternions. Therefore, it exists
:math:`\omega \in {\mbox{\rm $I\!\!R$}}^3` and
:math:`\Omega \in {\mbox{\rm $I\!\!R$}}^3` such that

.. math::

   \label{eq:132}
     (0, \Omega) = 2 p ^\star \glaw a

and

.. math::

   \label{eq:132}
     (0, \omega) = 2 a\glaw p^\star

In other terms, the tangent vector spaces at
:math:`p \in {\mbox{\rm $I\!\!H$}}_1` can be represented as a left
representation

.. math::

   \label{eq:133}
     T_p{\mbox{\rm $I\!\!H$}}_1 = \{ a \mid a = p \glaw (0, \frac \Omega 2 ), \Omega \in {\mbox{\rm $I\!\!R$}}^3\}

or a right representation

.. math::

   \label{eq:1330}
     T_p{{\mbox{\rm $I\!\!H$}}_1} = \{ a \mid a =  (0, \frac \omega 2 ) \glaw p, \omega \in {\mbox{\rm $I\!\!R$}}^3\}

At :math:`p=e`, we get the Lie algebra defined by

.. math::

   \label{eq:134}
     \mathfrak h_1 =  T_e{{\mbox{\rm $I\!\!H$}}_1} =  \{ a = (0, \frac \Omega 2 ), \Omega \in {\mbox{\rm $I\!\!R$}}^3 \}

equipped with the Lie bracket given by the commutator

.. math::

   \label{eq:135}
     [p,q] = p \glaw q- q\glaw p.

We can easily verify that for
:math:`a = (0, \frac \Omega 2 ), \, b = (0, \frac \Gamma 2 ) \in \mathfrak h_1 `,
we have

.. math::

   \label{eq:136}
     [a,b] = (0, \frac \Omega 2 ) \glaw (0, \frac \Gamma 2 ) - (0, \frac \Gamma 2 ) \glaw  (0, \frac \Omega 2 )  = (0, \frac{\Omega\times \Gamma}{2}) \in \mathfrak h_1

As for :math:`\mathfrak so(3)`, the Lie algebra :math:`\mathfrak h_1` is
isomorphic to :math:`{\mbox{\rm $I\!\!R$}}^3` thanks to the operator
:math:`\widehat{(\cdot)} :{\mbox{\rm $I\!\!R$}}^3 \rightarrow \mathfrak h_1`
and defined by

.. math::

   \label{eq:54}
    \widehat{(\cdot)}: \Omega \mapsto \widehat \Omega = (0, \frac \Omega 2 )

With this operator, the Lie Bracket can be written

.. math::

   \label{eq:137}
     [\widehat{\Omega},\widehat{\Gamma}] = \widehat{\Omega \times \Gamma}

 A special (right) action of Lie Group :math:`\mathcal G` on a manifold :math:`\mathcal M`. 
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Let us come back to the representation of
:math:`T_p{\mbox{\rm $I\!\!H$}}_1` given in . It is clear it can
expressed with a representation that relies on :math:`\mathfrak h_1`

.. math::

   \label{eq:158}
      T_RSO(3) = \{ a = p \glaw  \widehat \Omega \mid \widehat \Omega \in \mathfrak h_1 \}.

With , we see that there is a linear map that relates
:math:`T_p{\mbox{\rm $I\!\!H$}}_1` to :math:`\mathfrak h_1`. This linear
map defines a vector field. A special group action is defined by the
left translation map for a point :math:`p \in {\mbox{\rm $I\!\!H$}}_1`

.. math::

   \label{eq:159}
     \begin{array}[lcl]{rcl}
       L_p& :&   {\mbox{\rm $I\!\!H$}}_1 \rightarrow  {\mbox{\rm $I\!\!H$}}_1\\
          & &  q  \mapsto L_p(q) = p \glaw q\\
     \end{array}

which is diffeomorphism on :math:`{\mbox{\rm $I\!\!H$}}_1`. In that
case, we identify the manifold and the group. So, :math:`L_p` can be
viewed as a left or a right group action. We choose a right action. For
our application where
:math:`\mathcal G = \mathcal M = {\mbox{\rm $I\!\!H$}}_1` and
:math:`\Lambda^r(p,q) = L_{p}(q) =  p \glaw q `, we get

.. math::

   \label{eq:160}
      \lambda^r_{*}(a)(p) = \left. \frac{d}{dt} L_{p}(q(t)) \right|_{t=0}  = \left. \frac{d}{dt} p \glaw q(t) \right|_{t=0} =  p \glaw \dot q(0) = p  \glaw \dot q(0)  \in X(\mathcal M)

for a smooth curve :math:`q(t)` in :math:`{\mbox{\rm $I\!\!H$}}_1`.
Since :math:`q(\cdot)` is a smooth curve in
:math:`{\mbox{\rm $I\!\!H$}}_1`, :math:`\dot q(0)` is a tangent vector
at the point :math:`q(0)=I`, that is an element
:math:`a = \widehat  \Omega  \in \mathfrak h_1 ` defined by the
relation . Therefore, the vector field in is a tangent vector field and
we get

.. math::

   \label{eq:161}
     \dot p(t) = \lambda^r_{*}(a)(p(t)) = p(t)  \glaw \widehat \Omega

Exponential map :math:`\operatorname{expq}: \mathfrak h_1 \rightarrow {\mbox{\rm $I\!\!H$}}_1`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

We can directly apply Theorem [Theorem:solutionofLieODE] and we get that
the solution of

.. math::

   \label{eq:130}
     \begin{cases}
     \dot p(t) = \lambda^r_{*}(a)(p(t)) = p(t) \glaw \widehat \Omega \\
     p(0) = Rp_0
   \end{cases}

is

.. math::

   \label{eq:138}
     p(t) = p_0 \operatorname{expq}(t \widehat \Omega)

The exponential mapping
:math:`\operatorname{expq}: \mathfrak h_1 \rightarrow {\mbox{\rm $I\!\!H$}}_1`
can also be defined as
:math:`\operatorname{expq}(\widehat \Omega) = q(1)` where :math:`q (t)`
satisfies the differential equation

.. math::

   \label{eq:235}
     \dot q(t) = q(t) \cdot \widehat \Omega , \quad q (0) = e.

Using the quaternion product, the exponential map can be expressed as

.. math::

   \label{eq:232}
     \operatorname{expq}(t \widehat \Omega ) = \sum_{k=0}^\infty \frac{(t\widehat \Omega)^k}{k!}

since it is a solution of . A simple computation allows to check this
claim:

.. math::

   \label{eq:233}
      \frac{d}{dt}\operatorname{expq}(t \widehat \Omega ) = \sum_{k=1}^\infty  k t^{k-1} \frac{ \widehat \Omega ^k}{k!} =  \sum_{k=0}^\infty  t^{k} \frac{t \widehat \Omega ^k}{k!}\glaw  \widehat \Omega  =   \operatorname{expq}(t \widehat \Omega ) \glaw \widehat \Omega.

A closed form relation for the form the quaternion exponential can also
be found by noting that

.. math::

   \label{eq:140}
     \widehat \Omega ^2  = - \left(\frac \theta 2 \right)^2 e, \text{ and } \widehat \Omega ^3  = - \left(\frac \theta 2 \right)^2 \widehat \Omega.

A simple expansion of at :math:`t=1` equals

.. math::

   \label{eq:141}
     \begin{array}{lcl}
       \operatorname{expq}(\widehat \Omega ) &=& \sum_{k=0}^\infty \frac{(\widehat \Omega)^k}{k!}\\
                               &=& \sum_{k=0}^\infty \frac{(-1)^k}{(2k)!}\left(\frac \theta 2 \right)^{2k} e + \sum_{k=0}^\infty \frac{(-1)^k}{(2k+1)!} \left(\frac \theta 2 \right)^{2k+1} \widehat \Omega \\
                               &=& \cos(\frac \theta 2) e + \frac{\sin(\frac \theta 2)}{\frac \theta 2} \widehat \Omega \\
     \end{array}

that is

.. math::

   \label{eq:144}
     \operatorname{expq}(\widehat \Omega )  = (\cos(\frac \theta 2), \sin(\frac \theta 2) \frac{\Omega}{\theta}   ).

Adjoint representation
''''''''''''''''''''''

In the case of :math:`{\mbox{\rm $I\!\!H$}}_1`, the definition of the
operator :math:`\operatorname{Ad}` gives

.. math::

   \label{eq:121}
     \operatorname{Ad}_p(\widehat\Omega)  = p\glaw \widehat\Omega p^\star

and then mapping
:math:`\operatorname{ad}_{\widehat\Omega}(\widehat \Gamma)` is defined
by

.. math::

   \label{eq:56}
     \operatorname{ad}_{\widehat\Omega}(\widehat\Gamma) = \widehat \Omega \widehat\Gamma - \widehat \Gamma \widehat\Omega  =  [\widehat \Omega,\widehat \Gamma] = \widehat{\Omega \times \Gamma}.

Using the isomorphism between :math:`\mathfrak h_1` and
:math:`{\mbox{\rm $I\!\!R$}}^3`, we can use the the mapping
:math:`\operatorname{ad}_{\Omega}(\Gamma) : {\mbox{\rm $I\!\!R$}}^3\times{\mbox{\rm $I\!\!R$}}^3 \rightarrow {\mbox{\rm $I\!\!R$}}^3`
given by to get

.. math::

   \label{eq:145}
      \operatorname{ad}_{\widehat\Omega}(\widehat\Gamma) = \widehat{\Omega \times \Gamma} = \widehat{\operatorname{ad}_{\Omega}(\Gamma)} =  \widehat{\tilde \Omega \Gamma}

Differential of the exponential map :math:`\operatorname{dexpq}`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The differential of the exponential mapping, denoted by
:math:`\operatorname{dexpq}` is defined as the ’right trivialized’
tangent of the exponential map

.. math::

   \label{eq:71}
     \frac{d}{dt} (\operatorname{expq}(\widehat \Omega(t))) = \operatorname{dexpq}_{\widehat\Omega(t)}(\frac{d \widehat{\Omega}(t)}{dt}) \operatorname{expq}(\widehat\Omega(t))

An explicit expression of
:math:`\operatorname{dexp}_{\widehat\Omega}(\widehat\Gamma)` can also be
developed either by developing the expansion and .

.. math::

   \label{eq:168}
      \operatorname{dexpq}_{\widehat\Omega}(\Gamma) = \sum_{k=0}^\infty \frac{1}{(k+1)!} \operatorname{ad}_{\widehat\Omega}^k (\widehat\Gamma) = \widehat{T(\Omega)\Gamma}

Note that the time derivative in :math:`{\mbox{\rm $I\!\!R$}}^4` is not
differential mapping. The standard time derivative of
:math:`\operatorname{expq}` in the expression gives

.. math::

   \label{eq:171}
       \frac{d}{dt}\operatorname{expq}(\widehat \Gamma(t)) = (- \frac{\sin(\theta)}{\theta} \Omega^T\Gamma, \frac{\sin(\theta)}{\theta}\Gamma  +\frac{\theta \cos(\theta)-\sin(\theta)}{\theta^3}\Omega^T\Omega \Gamma  )

that can be expressed in :math:`{\mbox{\rm $I\!\!R$}}^4` by

.. math::

   \label{eq:175}
     \frac{d}{dt}\operatorname{expq}(\widehat \Gamma(t))  = \nabla \operatorname{expq}(\widehat\Omega) \widehat{\dot\Omega}

with

.. math::

   \label{eq:176}
     \nabla \operatorname{expq}(\widehat\Omega) =
     \begin{bmatrix}
       - \frac{\sin(\theta)}{\theta} \Omega^T \\
       \frac{\sin(\theta)}{\theta}I  +\frac{\theta \cos(\theta)-\sin(\theta)}{\theta^3}\Omega^T\Omega
     \end{bmatrix}

Clearly, we have

.. math::

   \label{eq:177}
     \nabla \operatorname{expq}(\widehat\Omega) \neq  \operatorname{dexpq}_{\widehat\Omega}

Directional derivative and Jacobians of functions of a quaternion
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

experimental

Let
:math:`f : {\mbox{\rm $I\!\!H$}}_1 \rightarrow {\mbox{\rm $I\!\!R$}}` be
a mapping from the group to :math:`{\mbox{\rm $I\!\!R$}}^3`. The
directional derivative of :math:`f` in the direction
:math:`\widehat \Omega \in \mathfrak h_1` at
:math:`p\in {\mbox{\rm $I\!\!H$}}_1` is defined by

.. math::

   \label{eq:139}
    df_p(\widehat \Omega) =\left. \frac{d}{dt} f(p\glaw \operatorname{expq}(t\widehat \Omega)) \right|_{t=0}

As a first simple example let us choose
:math:`f(p) = \vv{p \glaw p_x \glaw p^\star}` for a given
:math:`x \in {\mbox{\rm $I\!\!R$}}^3 `, we get

.. math::

   \label{eq:142}
     \begin{array}{lcl}
       D Id \cdot \widehat \Omega (p) = (\widehat \Omega^r f )(p) &=& \left. \frac{d}{dt}\vv{p\glaw \operatorname{expq}(t\widehat \Omega) \glaw p_x \glaw (p \glaw \operatorname{expq}(t\widehat \Omega))^\star}  \right|_{t=0}\\
                                                                  & = & \vv{p\glaw \frac{d}{dt}\left. \operatorname{expq}(t\widehat \Omega) \right|_{t=0} \glaw p_x \glaw p^\star +  p \glaw p_x \glaw (p \glaw\frac{d}{dt}\left. \operatorname{expq}(t\widehat \Omega) \right|_{t=0})^\star}\\                                              
     \end{array}

We have form the definition of the time derivative of the exponential

.. math::

   \label{eq:143}
     \begin{array}{lcl}
       \frac{d}{dt}\left. \operatorname{expq}(t\widehat \Omega) \right|_{t=0} &=&  \left. \operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega)\operatorname{expq}(t\widehat \Omega) \right|_{t=0} \\
                                                               &=&  \operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega)
     \end{array}

Then, the directional derivative can be written

.. math::

   \label{eq:146}
     \begin{array}{lcl}
       D Id \cdot \widehat \Omega (p) &=& \vv{p\glaw \operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega)\glaw p_x \glaw p^\star  + p \glaw p_x \glaw (\operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega))^* \glaw  p^\star } \\
     &=& \vv{p\glaw ( \operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega)\glaw p_x +   p_x \glaw (\operatorname{dexpq}_{\widehat\Omega}(\widehat \Omega))^*) \glaw  p^\star } 
     \end{array}

Newton-Euler equation in quaternion form
----------------------------------------

Computation of :math:`T` for unit quaternion
''''''''''''''''''''''''''''''''''''''''''''

The operator :math:`T(q)` is directly obtained as

.. math::

   T(q)=\frac 1 2 \label{eq:98}
     \begin{bmatrix}
       2 I_{3\times 3} & & 0_{3\times 3} & \\
       &   -p_1 & -p_2 & -p_3 \\
       0_{4\times 3}  &  p_0 & -p_3 & p_2 \\
       & p_3 & p_0 & -p_1 \\
       & -p_2 & p_1 & p_0 
     \end{bmatrix}

todo :

-  computation of the directional derivative of
   :math:`R(\Omega)= exp(\tilde \Omega)` in the direction
   :math:`\tilde\Omega`, to get :math:`T(\Omega)`

Quaternion representation
'''''''''''''''''''''''''

If the Lie group is described by unit quaternion, we get

.. math::

   \label{eq:99}
     SO(3) = \{p = (p_0,\vv{p}) \in {\mbox{\rm $I\!\!R$}}^{4}\mid |p|=1  \}

with the composition law :math:`p_1\glaw p_2` given by the quaternion
product.

Note that the concept of exponential map for Lie group that are not
parameterized by matrices is also possible.

Mechanical systems with bilateral and unilateral constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us consider that the system ([eq:Newton-Euler-compact]) is subjected
to :math:`m` constraints, with :math:`m_{e}` holonomic bilateral
constraints

.. math::

   \label{eq:bilateral-constraints}
     h^\alpha(q)=0, \alpha \in \mathcal{E}\subset{\mbox{\rm $I\!\!N$}},  |\mathcal E| = m_e,

and :math:`m_{i}` unilateral constraints

.. math::

   \label{eq:unilateral-constraints}
     g_{\n}^\alpha(q)\geq 0, \alpha \in \mathcal{I}\subset{\mbox{\rm $I\!\!N$}},  |\mathcal I| = m_i.

Let us denote as :math:`J^\alpha_h(q) = \nabla^\top_q h^\alpha(q)  ` the
Jacobian matrix of the bilateral constraint :math:`h^\alpha(q)` with
respect to :math:`q` and as :math:`J^\alpha_{g_\n}(q)` respectively for
:math:`g_{\n}^\alpha(q)` . The bilateral constraints at the velocity
level can be obtained as:

.. math::

   \label{eq:bilateral-constraints-velocity}
    0 = \dot h^\alpha(q)= J^\alpha_h(q)\dot q = J^\alpha_h(q) T(q) v \coloneqq H^\alpha(q)  v,\quad  \alpha \in \mathcal{E}.

By duality and introducing a Lagrange multiplier
:math:`\lambda^\alpha, \alpha \in \mathcal E`, the constraint generates
a force applied to the body equal to
:math:`H^{\alpha,\top}(q)\lambda^\alpha`. For the unilateral
constraints, a Lagrange multiplier
:math:`\lambda_{\n}^\alpha, \alpha \in \mathcal I` is also associated
and the constraints at the velocity level can also be derived as

.. math::

   \label{eq:unilateral-constraints-velocity}
    0 \leq  \dot g_\n^\alpha(q)= J^\alpha_{g_\n}(q) \dot q = J^\alpha_{g_\n}(q)  T(q) v , \text{ if } g_{\n}^\alpha(q) = 0,\quad  \alpha \in \mathcal{I}.

Again, the force applied to the body is given by
:math:`(J^\alpha_{g_\n}(q) T(q))^\top\lambda^\alpha_\n`. Nevertheless,
there is no reason that :math:`\lambda^\alpha_\n =r^\alpha_\n` and
:math:`u_\n = J^\alpha_{g_\n}(q) T(q) v` if the :math:`g_n` is not
chosen as the signed distance (the gap function). This is the reason why
we prefer directly define the normal and the tangential local relative
velocity with respect to the twist vector as

.. math::

   \label{eq:unilateral-constraints-velocity-kinematic1}
      u^\alpha_\n  \coloneqq G_\n^\alpha(q) v, \quad u^\alpha_\t  \coloneqq G_\t^\alpha(q) v, \quad \alpha \in \mathcal{I},

and the associated force as :math:`G_\n^{\alpha,\top}(q) r^{\alpha}_\n `
and :math:`G_\t^{\alpha,\top}(q) r^{\alpha}_\t`. For the sake of
simplicity, we use the notation
:math:`u^\alpha  \coloneqq G^\alpha(q) v` and its associated total force
generated by the contact :math:`\alpha` as
:math:`G^{\alpha,\top}(q) r^{\alpha} \coloneqq G_\n^{\alpha,\top}(q) r^{\alpha}_\n + G_\t^{\alpha,\top}(q) r^{\alpha}_\t `.

The complete system of equation of motion can finally be written as

|    q = T(q)v ,
|    M v = F(t,q,v) + H\ :sup:``\ (q) + G\ :sup:``\ (q) r,
| [0.5ex]   

| ll H\ :sup:``\ (q) v = 0 ,& E
| .

| ll r\ :sup:``\ = 0 , & g\ :sub:``\ :sup:``\ (q) > 0,
| K\ :sup:`,\*` u\ :sup:``  r\ :sup:``, & g\ :sub:``\ :sup:``\ (q) = 0,
| u\ :sub:``\ :sup:`,+` = -e:sub:`r`\ :sup:``\ u\ :sub:``\ :sup:`,-`, &
g\ :sub:``\ :sup:``\ (q) = 0 u\ :sub:``\ :sup:`,-` 0,

} & I [eq:NewtonEuler-uni]

where the definition of the variables
:math:`\lambda\in {\mbox{\rm $I\!\!R$}}^{m_e}, r\in {\mbox{\rm $I\!\!R$}}^{3m_i}`
and the operators :math:`H,G` are extended to collect all the variables
for each constraints.

Note that all the constraints are written at the velocity integrators.
Another strong advantage is the straightforward introduction of the
contact dissipation processes that are naturally written at the velocity
level such as the Newton impact law and the Coulomb friction. Indeed, in
Mechanics, dissipation processes are always given in terms of rates of
changes, or if we prefer, in terms of velocities.

Siconos Notation
''''''''''''''''

In the siconos notation, we have for the applied torques on the system
the following decomposition

.. math::

   F(t,q,v):= \begin{pmatrix}
       f(t,x_{\cg},  v_{\cg}, R, \Omega ) \\
       I \Omega \times \Omega + M(t,x_{\cg}, v_{\cg}, R, \Omega )
     \end{pmatrix}
     := \begin{pmatrix}
       f_{ext}(t)  - f_{int}(x_{\cg},  v_{\cg}, R, \Omega ) \\
       - M_{gyr}(\Omega) + M_{ext}(t) -  M_{int}(x_{\cg}, v_{\cg}, R, \Omega )
     \end{pmatrix}.

with

.. math::

   M_{gyr} := \begin{pmatrix}
        \Omega \times I\Omega
     \end{pmatrix}

In the siconos notation, we have for the relation

.. math::

   \label{eq:100}
      C =   J^\alpha(q) \quad CT = J^\alpha(q)T(q)

Time integration scheme in scheme
---------------------------------

Moreau–Jean scheme based on a :math:`\theta`-method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The complete Moreau–Jean scheme based on a :math:`\theta`-method is
written as follows

.. math::

   \label{eq:Moreau--Jean-theta}
       \begin{cases}
         ~~\begin{array}{l}
           q_{k+1} = q_{k} + h T(q_{k+\theta}) v_{k+\theta} \quad \\[1ex]
           M(v_{k+1}-v_k) - h  F(t_{k+\theta}, q_{k+\theta},v_{k+\theta}) =  H^\top(q_{k+1}) Q_{k+1} + G^\top(q_{k+1}) P_{k+1},\quad\,\\[1ex]
         \end{array}\\
         ~~\begin{array}{lcl}
           \begin{array}{l}
             H^\alpha(q_{k+1}) v_{k+1}  =  0\\
           \end{array} & \left. \begin{array}{l}
             \vphantom{H^\alpha(q_{k+1}) v_{k+1}  =  0}\\[1ex]
           \end{array}\right\}    &\alpha \in \mathcal E  \\[1ex]
         ~~~P_{k+1}^\alpha= 0, &
         \left. \begin{array}{l}
             \vphantom{P_{k+1}^\alpha= 0,  \delta^\alpha_{k+1}=0}\\[1ex]
           \end{array}\right\}   & \alpha \not\in \mathcal I^\nu \\[1ex]
                     \begin{array}{l}
             {K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha \\
         \end{array} &
         \left.\begin{array}{l}
             \vphantom{{K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha} \\
           \end{array}\right\}
         &\alpha \in \mathcal I^\nu\\
     \end{array}
   \end{cases}

where :math:`\mathcal I^\nu` is the set of forecast constraints, that
may be evaluated as

.. math::

   \label{eq:101}
     \mathcal I^\nu = \{\alpha \mid \bar g_\n^\alpha \coloneqq g_\n + \frac h 2 u^\alpha_\n \leq 0\}.

Semi-explicit version Moreau–Jean scheme based on a :math:`\theta`-method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   \label{eq:Moreau--Jean-explicit}
       \begin{cases}
         ~~\begin{array}{l}
           q_{k+1} = q_{k} + h T(q_{k}) v_{k+\theta} \quad \\[1ex]
           M(v_{k+1}-v_k) - h  F(t_{k}, q_{k},v_{k}) =  H^\top(q_{k}) Q_{k+1}+  G^\top(q_{k}) P_{k+1},\quad\,\\[1ex]
         \end{array}\\
         ~~\begin{array}{lcl}
           \begin{array}{l}
             H^\alpha(q_{k+1}) v_{k+1}  =  0\\
           \end{array} & \left. \begin{array}{l}
             \vphantom{H^\alpha(q_{k+1}) v_{k+1}  =  0}\\[1ex]
           \end{array}\right\}    &\alpha \in \mathcal E  \\[1ex]
         ~~P_{k+1}^\alpha= 0, &
         \left. \begin{array}{l}
             \vphantom{P_{k+1}^\alpha= 0,  \delta^\alpha_{k+1}=0}\\[1ex]
           \end{array}\right\}   & \alpha \not\in \mathcal I^\nu \\[1ex]
                     \begin{array}{l}
             {K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha \\
         \end{array} &
         \left.\begin{array}{l}
             \vphantom{{K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha} \\
           \end{array}\right\}
         &\alpha \in \mathcal I^\nu\\
     \end{array}
   \end{cases}

In this version, the new velocity :math:`v_{k+1}` can be computed
explicitly, assuming that the inverse of :math:`M` is easily written, as

.. math::

   \label{eq:Moreau--Jean-theta--explicit-v}
     v_{k+1}   =  v_k + M^{-1} h  F(t_{k}, q_{k},v_{k}) +  M^{-1} (H^\top(q_{k}) Q_{k+1}+  G^\top(q_{k}) P_{k+1})

Nearly implicit version Moreau–Jean scheme based on a :math:`\theta`-method implemented in siconos
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A first simplification is made considering a given value of
:math:`q_{k+1}` in :math:`T()`, :math:`H()` and :math:`G()` denoted by
:math:`\bar q_k`. This limits the computation of the Jacobians of this
operators with respect to :math:`q`.

.. math::

   \label{eq:Moreau--Jean-theta-nearly}
       \begin{cases}
         ~~\begin{array}{l}
           q_{k+1} = q_{k} + h T(\bar q_k) v_{k+\theta} \quad \\[1ex]
           M(v_{k+1}-v_k) - h  \theta F(t_{k+1}, q_{k+1},v_{k+1}) - h (1- \theta) F(t_{k}, q_{k},v_{k})  =  H^\top(\bar q_k) Q_{k+1} + G^\top(\bar q_k) P_{k+1},\quad\,\\[1ex]
         \end{array}\\
         ~~\begin{array}{lcl}
           \begin{array}{l}
             H^\alpha(\bar q_k) v_{k+1}  =  0\\
           \end{array} & \left. \begin{array}{l}
             \vphantom{H^\alpha(q_{k+1}) v_{k+1}  =  0}\\[1ex]
           \end{array}\right\}    &\alpha \in \mathcal E  \\[1ex]
         ~~P_{k+1}^\alpha= 0, &
         \left. \begin{array}{l}
             \vphantom{P_{k+1}^\alpha= 0,  \delta^\alpha_{k+1}=0}\\[1ex]
           \end{array}\right\}   & \alpha \not\in \mathcal I^\nu \\[1ex]
                     \begin{array}{l}
             {K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha \\
         \end{array} &
         \left.\begin{array}{l}
             \vphantom{{K}^{\alpha,*} \ni \widehat u_{k+1}^\alpha~ \bot~ P_{k+1}^\alpha \in {K}^\alpha} \\
           \end{array}\right\}
         &\alpha \in \mathcal I^\nu\\
     \end{array}
   \end{cases}

The nonlinear residu is defined as

.. math::

   \label{eq:Moreau--Jean-theta--nearly-residu}
     \mathcal R(v) =  M(v-v_k) - h  \theta F(t_{k+1}, q(v),v) - h (1- \theta) F(t_{k}, q_{k},v_{k}) - H^\top(\bar q_k) Q_{k+1} - G^\top(\bar q_k) P_{k+1}

with

.. math::

   \label{eq:Moreau--Jean-theta--nearly-residu1}
     q(v) = q_{k} + h T(\bar q_k)) ((1-\theta) v_k + \theta v).

At each time step, we have to solve

.. math::

   \label{eq:Moreau--Jean-theta--nearly-residu2}
     \mathcal R(v_{k+1}) =  0

together with the constraints.

Let us write a linearization of the problem to design a Newton
procedure:

.. math::

   \label{eq:Moreau--Jean-theta--nearly-residu3}
     \nabla^\top_v \mathcal R(v^{\tau}_{k+1})(v^{\tau+1}_{k+1}-v^{\tau}_{k+1}) = -  \mathcal R(v^{\tau}_{k+1}).

The computation of :math:` \nabla^\top_v \mathcal R(v^{\tau}_{k+1})` is
as follows

.. math::

   \label{eq:102}
     \nabla^\top_v \mathcal R(v) = M - h \theta \nabla_v F(t_{k+1}, q(v),v)

with

.. math::

   \label{eq:103}
     \begin{array}{lcl}
       \nabla_v F(t_{k+1}, q(v),v) &=& D_2 F(t_{k+1}, q(v),v) \nabla_v q(v) + D_3 F(t_{k+1}, q(v),v) \\
                                   &=& h \theta D_2 F(t_{k+1}, q(v),v) T(\bar q_k) + D_3 F(t_{k+1}, q(v),v) \\
     \end{array}

where :math:`D_i` denotes the derivation with respect the :math:`i^{th}`
variable. The complete Jacobian is then given by

.. math::

   \label{eq:104}
     \nabla^\top_v \mathcal R(v) = M - h \theta D_3 F(t_{k+1}, q(v),v) - h^2 \theta^2 D_2 F(t_{k+1}, q(v),v) T(\bar q_k)

In siconos, we ask the user to provide the functions
:math:`D_3 F(t_{k+1}, q ,v )` and :math:`D_2 F(t_{k+1}, q,v)`.

Let us denote by :math:`W^{\tau}` the inverse of Jacobian of the residu,

.. math::

   \label{eq:105}
     W^{\tau} = (M - h \theta D_3 F(t_{k+1}, q(v),v) - h^2 \theta^2 D_2 F(t_{k+1}, q(v),v) T(\bar q_k))^{-1}.

and by :math:`\mathcal R_{free}(v)` the free residu,

.. math::

   \label{eq:106}
     \mathcal R_{free}(v) =  M(v-v_k) - h  \theta F(t_{k+1}, q(v),v) - h (1- \theta) F(t_{k}, q_{k},v_{k}).

The linear equation [eq:Moreau–Jean-theta–nearly-residu3] that we have
to solve is equivalent to

.. math::

   \label{eq:107}
     \boxed{v^{\tau+1}_{k+1} = v^{\tau}_{k+1} - W  \mathcal R_{free}(v^\tau_{k+1}) + W   H^\top(\bar q_k) Q^{\tau+1}_{k+1} + W G^\top(\bar q_k) P^{\tau+1}_{k+1}}

We define :math:`v_{free}` as

.. math::

   \label{eq:108}
     v_{free}  = v^{\tau}_{k+1} - W  \mathcal R_{free}(v^\tau_{k+1})

The local velocity at contact can be written

.. math::

   \label{eq:109}
     u^{\tau+1}_{\n,k+1} = G(\bar q_k) [  v_{free}^{\tau} + W   H^\top(\bar q_k) Q^{\tau+1}_{k+1} + W G^\top(\bar q_k) P^{\tau+1}_{k+1}]

and for the equality constraints

.. math::

   \label{eq:110}
     u^{\tau+1}_{k+1} = H(\bar q_k) [  v_{free}^{\tau} + W   H^\top(\bar q_k) Q^{\tau+1}_{k+1} + W G^\top(\bar q_k) P^{\tau+1}_{k+1}]

Finally, we get a linear relation between :math:`u^{\tau+1}_{\n,k+1}`
and the multiplier

.. math::

   \label{eq:111}
    \boxed{ u^{\tau+1}_{k+1} =
     \begin{bmatrix}
       H(\bar q_k) \\
       G(\bar q_k)
     \end{bmatrix} v_{free}^{\tau}
     +
     \begin{bmatrix}
       H(\bar q_k)W   H^\top(\bar q_k) & H(\bar q_k)W   G^\top(\bar q_k) \\
       G(\bar q_k)W   H^\top(\bar q_k) & G(\bar q_k)W   G^\top(\bar q_k) \\
     \end{bmatrix}
     \begin{bmatrix}
       Q^{\tau+1}_{k+1} \\
       P^{\tau+1}_{k+1}
     \end{bmatrix}}

choices for :math:`\bar q_k`
''''''''''''''''''''''''''''

Two choices are possible for :math:`\bar q_k`

#. :math:`\bar q_k = q_k`

#. :math:`\bar q_k = q^{\tau}_{k+1}`

todo list:

-  add the projection step for the unit quaternion

-  describe the computation of H and G that can be hybrid

Computation of the Jacobian in special case
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Moment of gyroscopic forces
'''''''''''''''''''''''''''

Let us denote by the basis vector :math:`e_i` given the :math:`i^{th}`
column of the identity matrix :math:`I_{3\times3}`. The Jacobian of
:math:`M_{gyr}` is given by

.. math::

   \label{eq:112}
     \nabla^\top_\Omega M_{gyr}(\Omega) = \nabla^\top_\Omega (\Omega \times I \Omega) =
     \begin{bmatrix}
       e_i \times I \Omega + \Omega \times I e_i, i =1,2,3
     \end{bmatrix}

Linear internal wrench
''''''''''''''''''''''

If the internal wrench is given by

.. math::

   \label{eq:113}
     F_{int}(t,q,v) =
     \begin{bmatrix}
       f_{int}(t,q,v)\\
       M_{int}(t,q,v)
     \end{bmatrix}
     = C v + K q, \quad C \in {\mbox{\rm $I\!\!R$}}^{6\times 6}, \quad K \in {\mbox{\rm $I\!\!R$}}^{6\times 7 }

we get

.. math::

   \label{eq:114}
     \begin{array}{lcl}
       \nabla_v F(t_{k+1}, q(v),v)  &=& h \theta K T(\bar q_k) + C \\
       \nabla^\top_v \mathcal R(v) &=& M - h \theta C - h^2 \theta^2 K T(\bar q_k)
     \end{array}

External moment given in the inertial frame
'''''''''''''''''''''''''''''''''''''''''''

If the external moment denoted by :math:`m_{ext} (t)` is expressed in
inertial frame, we have

.. math::

   \label{eq:115}
     M_{ext}(q,t) = R^T m_{ext}(t)= \Phi(p) m_{ext}(t)

In that case, :math:`  M_{ext}(q,t)` appears as a function :math:`q` and
we need to compute its Jacobian w.r.t :math:`q`. This computation needs
the computation of

.. math::

   \label{eq:116}
     \nabla_{p} M_{ext}(q,t) = \nabla_{p} \Phi(p) m_{ext}(t)

Let us compute first

.. math::

   \label{eq:117}
     \Phi(p) m_{ext}(t)  =
     \begin{bmatrix}
       (1-2 p_2^2- 2 p_3^2)m_{ext,1} + 2(p_1p_2-p_3p_0)m_{ext,2} + 2(p_1p_3+p_2p_0)m_{ext,3}\\
       2(p_1p_2+p_3p_0)m_{ext,1}  +(1-2 p_1^2- 2 p_3^2)m_{ext,2} + 2(p_2p_3-p_1p_0)m_{ext,3}\\
       2(p_1p_3-p_2p_0)m_{ext,1}  + 2(p_2p_3+p_1p_0)m_{ext,2}  + (1-2 p_1^2- 2 p_2^2)m_{ext,3}\\
     \end{bmatrix}

Then we get

.. math::

   \label{eq:118}
     \begin{array}{l}
     \nabla_{p} \Phi(p) m_{ext}(t)  =\\
     \begin{bmatrix}
       -2 p_3 m_{ext,2} + 2 p_2 m_{ext,3} & 2p_2 m_{ext,2}+2 p_3 m_{ext,3}  & -4 p_2 m_{ext,1} +2p_1 m_{ext,2}+2 p_0 m_{ext,3} & -3 p_3 m_{ext,1} -2p_0 m_{ext,2} +2 p_1m_{ext,3}  \\
       2p_3 m_{ext,1} -2p_1m_{ext,3}  & 2p_2m_{ext,1} -4p_1 m_{ext,2} -2p_1 m_{ext,3} & & &  \\
     \end{bmatrix}
     \end{array}

Siconos implementation
~~~~~~~~~~~~~~~~~~~~~~

The
expression: \ :math:`\mathcal R_{free}(v^\tau_{k+1}) = M(v-v_k) - h  \theta F(t_{k+1}, q(v^\tau_{k+1}),v^\tau_{k+1}) - h (1- \theta) F(t_{k}, q_{k},v_{k})`
is computed in MoreauJeanOSI::computeResidu() and saved in
ds->workspace(DynamicalSystem::freeresidu)

The
expression: \ :math:`\mathcal R(v^\tau_{k+1}) =\mathcal R_{free}(v^\tau_{k+1}) - h (1- \theta) F(t_{k}, q_{k},v_{k}) - H^\top(\bar q_k) Q_{k+1} - G^\top(\bar q_k) P_{k+1}  `
is computed in MoreauJeanOSI::computeResidu() and saved in
ds->workspace(DynamicalSystem::free).

really a bad name for the buffer ds->workspace(DynamicalSystem::free).
Why we are chosing this name ? to save some memory ?

The
expression: \ :math:`v_{free}  = v^{\tau}_{k+1} - W  \mathcal R_{free}(v^\tau_{k+1})`
is compute in MoreauJeanOSI::computeFreeState() and saved in
d->workspace(DynamicalSystem::free).

| The computation: 
:math:`v^{\tau+1}_{k+1} = v_{free} + W   H^\top(\bar q_k) Q^{\tau+1}_{k+1} + W G^\top(\bar q_k) P^{\tau+1}_{k+1}`
is done in MoreauJeanOSI::updateState and stored in d->twist().

NewtonEulerR: computation of :math:`\nabla _qH`
===============================================

Gradient computation, case of NewtonEuler with quaternion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the section, :math:`q` is the quaternion of the dynamical system.

[figCase]

The normal vector :math:`N` is view as a constant.

.. math:: ~\tilde h(q)=P_c(\frac{q}{\|q\|})

.. math:: ^t \nabla h(q)(\delta q) = \lim _{e \to 0}\frac{(\tilde h (q+e\delta q)-\tilde h (q)).N}{e}

:math:`\nabla _q h` consist in computing
:math:`P_c(\frac{q+\delta q}{\|q+\delta q\|})-P_c(q)`.

.. math:: GP(q)=qG_0P_0~^cq

.. math:: GP(\frac{q+\delta q}{\|q+\delta q\|})=(q+\delta q)G_0P_0~^c(q+\delta q)\frac{1}{\|q+\delta q\|^2}

.. math:: =(q+\delta q)~^cqGP(q)q~^c(q+\delta q)\frac{1}{\|q+\delta q\|^2}

.. math:: =((1,0,0,0)+\delta q~^cq)GP(q)((1,0,0,0)+q~^c\delta q)\frac{1}{\|q+\delta q\|^2}

.. math:: =GP(q)+\delta q~^cqGP(q) + GP(q)q~^c\delta q+0(\delta q)^2\frac{1}{\|q+\delta q\|^2}

 So, because G is independant of :math:`q`:

.. math:: P(\frac{q+\delta q}{\|q+\delta q\|})-P(q)=qGP(\frac{q+\delta q}{\|q+\delta q\|})-GP(q)=\delta q~^cqGP(q) + GP(q)q~^c\delta q+0(\delta q)^2 + GP(q)\frac{1}{\|q+\delta q\|^2}

 For the directional derivation, we chose
:math:`\delta q = \epsilon * (1,0,0,0)`. using a equivalent to
:math:`\frac{1}{1+\epsilon}`

.. math:: \lim_{\epsilon \to 0}\frac{P(\frac{q+\delta q}{\|q+\delta q\|})-P(q)}{\epsilon}=~^cqGP(q) + GP(q)q-2q_iGP(q)

 For the directional derivation, we chose
:math:`\delta q = \epsilon * (0,1,0,0)=\epsilon * e_i`

.. math:: \lim_{\epsilon \to 0}\frac{P(\frac{q+\delta q}{\|q+\delta q\|})-P(q)}{\epsilon}=e_i~^cqGP(q) - GP(q)qe_i-2q_iGP(q)

 Application to the NewtonEulerRImpact:

.. math:: H:\mathbb{R}^7 \to \mathbb{R}

.. math:: \nabla _q H \in \mathcal{M}^{1,7}

.. math::

   \nabla _q H =\left(\begin{array}{c} N_x\\N_y\\N_z\\
   (~^cqGP(q) + GP(q)q-2q_0GP(q)).N\\
   (e_2~^cqGP(q) - GP(q)qe_2-2q_1GP(q)).N\\
   (e_3~^cqGP(q) - GP(q)qe_3-2q_2GP(q)).N\\
   (e_4~^cqGP(q) - GP(q)qe_4-2q_3GP(q)).N\\
   \end{array}\right)

Ball case
~~~~~~~~~

It is the case where :math:`GP=-N`: for :math:`e2`:

.. math:: (0,1,0,0).(q_0,-\underline p).(0,-N)=

.. math:: \left(\left(\begin{array}{c}1\\0\\0\end{array}\right).\underline p,\left(\begin{array}{c}q_0\\0\\0\end{array}\right) -\left(\begin{array}{c}1\\0\\0\end{array}\right)*\underline p \right).(0,-N)=

.. math:: \left(?, -\underline p_x~N-\left(\left(\begin{array}{c}q_0\\0\\0\end{array}\right)- \left(\begin{array}{c}1\\0\\0\end{array}\right)*\underline p \right)*N\right)=

 and:

.. math:: (0,-N).(q_0,\underline p).(0,1,0,0)=

.. math:: (N.\underline p,-q_0N-N*\underline p).(0,1,0,0)=

.. math:: \left(?,(N.\underline p)\left(\begin{array}{c}1\\0\\0\end{array}\right) + \left(\begin{array}{c}1\\0\\0\end{array}\right)*(q_0N+N*\underline p)\right)=

.. math:: \left(?,(N.\underline p)\left(\begin{array}{c}1\\0\\0\end{array}\right)+q_0 \left(\begin{array}{c}1\\0\\0\end{array}\right)*N+\left(\begin{array}{c}1\\0\\0\end{array}\right)*(N*\underline p)\right)

 sub then and get the resulting vector.N:

.. math:: \left[ -\underline p_x~N -N.\underline p~\left(\begin{array}{c}1\\0\\0\end{array}\right)+()*N-\left(\begin{array}{c}1\\0\\0\end{array}\right)*(N*\underline p)\right].N=

.. math:: -\underline p_x-N_xN.\underline p+0-(\left(\begin{array}{c}1\\0\\0\end{array}\right)*(N*\underline p)).N=

 using :math:`a*(b*c)=b(a.c)-c(a.b)` leads to

.. math:: -q_1-N_xN.\underline p-(q_1~N-N_x~\underline p).N=

.. math:: -q_1-N_xN.\underline p-q_1+N_xN.\underline p=-2q_1

 for :math:`e1=(1,0,0,0)`:

.. math:: (q_0,-\underline p).(0,-N)=(?,-q_0N+\underline p*N)

.. math:: (0,-N).(q_0,\underline p)=(?,-q_0N-\underline p*N)

 So

.. math::

   \nabla _q H =\left(\begin{array}{c} N_x\\N_y\\N_z\\
   0\\
   0\\
   0\\
   0\\
   \end{array}\right)

Case FC3D: using the local frame and momentum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math:: \left(\begin{array}{c}m \dot V\\I \dot \Omega + \Omega I \Omega \end{array}\right)= \left(\begin{array}{c}Fect+R\\Mext _{R_{obj}} + (R*PG) _{R_{obj}} \end{array}\right)

with \* vectoriel product, :math:`R` reaction in the globla frame.
:math:`P` the point of contact. :math:`r` is the reaction in the local
frame. :math:`M_{R_{obj}toR_{abs}}=M_{R_{abs}toR_{obj}}^t r=R` with:

.. math:: M_{R_{C}toR_{abs}}=\left(\begin{array}{ccc} nx&t_1x&t_2x \\ny&t_1y&t_2y\\nz&t_1z&t_2z \end{array}\right)

 we have :

.. math:: \left(\begin{array}{c}R\\(R*PG) _{R_{obj}}\end{array}\right)=\left(\begin{array}{c} I_3\\M_{R_{abs}toR_{obj}}N_{PG}\end{array}\right).R

.. math:: =\left(\begin{array}{c} I_3\\M_{R_{abs}toR_{obj}}N_{PG}\end{array}\right).M_{R_{obj}toR_{abs}}r

.. math:: N_{PG}=\left(\begin{array}{ccc} 0&PG_z&-PG_y\\-PG_z&0&PG_x\\PG_y&-PG_X&0\end{array}\right)

 that is:

.. math::

   \left(\begin{array}{c}m \dot V\\I \dot \Omega + \Omega I \Omega \end{array}\right)=
   \left(\begin{array}{c} M_{R_{C}toR_{abs}} \\
     M_{R_{abs}toR_{obj}}N_{PG}M_{R_{C}toR_{abs}}
   \end{array}\right) r

So :math:`jachqt=MN`

Case FC3D: using the local frame local velocities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[figCase]

We are looking for an operator named :math:`CT` such that:

.. math:: V_C=\left(\begin{array}{c} V_N \\ V_T \\ V_S \end{array}\right)_{R_{C}}=CT \left(\begin{array}{c} V_{G1}~_{R_{abs}} \\ \Omega_1~_{R_{obj1}} \\ V_{G2}~_{R_{abs}}\\ \Omega_2~_{R_{obj2}} \end{array}\right)

.. math:: V_c=V_{G1}~_{R_{abs}} + w_1 * G_1P~_{R_{abs}} -(V_{G2}~_{R_{abs}} + w_2 * G_1P~_{R_{abs}})

where :math:`w_1` and :math:`w_2` are given in :math:`R_{abs}`. We note
:math:`M_{R_{obj1}toR_{abs}}` the matrice converting the object 1
coordinate to the absolute coordinate. We note :math:`N_{GP}` the
matrice such that :math:`w_1*G_1P~_{R_{abs}} = N_{GC} w_1`. Endly, we
note :math:`M_{R_{abs}toR_C}` converting the absolute coordinate to the
:math:`R_C` frame. we get:

.. math:: CT= M_{R_{abs}toR_C}   \left(\begin{array}{cccc} I_3 & N_{G_1C}M_{R_{obj1}toR_{abs}} & -I_3 & -N_{G_2C}M_{R_{obj2}toR_{abs}} \end{array}\right)

Expression of :math:`M_{R_{obj1}toR_{abs}}`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using quaternion, we get :

.. math::

   \label{eq:newton_Mobjtoabs}
   M_{R_{obj1}toR_{abs}} = \left(\begin{array}{ccc} q \left(\begin{array}{c}1\\0\\0 \end{array}\right)~^cq & q \left(\begin{array}{c}  0\\1\\0 \end{array}\right)~ ^cq & q \left(\begin{array}{c}  0\\0\\1 \end{array}\right)~ ^cq  \end{array}\right)

Expression of :math:`N_1`
^^^^^^^^^^^^^^^^^^^^^^^^^

.. math:: N_{GC}=\left(\begin{array}{ccc} 0&G_1C_z&-G_1C_y\\-G_1C_z&0&G_1C_x\\G_1C_y&-G_1C_X&0\end{array}\right)

Projection On constraints
=========================

Velocity formulation
~~~~~~~~~~~~~~~~~~~~

The first step consists in doing a velocity formulation of the system:

.. math::

   \label{NE_Dyn1}
   \begin{array}{l}
     M \dot v = F_{ext}+B \lambda \\
     \dot q = T v\\
     y=h(q) \\
     NSLAW(y,\lambda,...)\\
   \end{array}

The constraint :math:`\dot q = T v` is suffisiant to keep a normal
quaternion. Because of the speed formulation, :math:`h(q)` could violate
the NSLAW. A solution coul be to add a formulation in position. We must
underline that the constraints :math:`\mid Q \mid = 1` is implicit in
this system. Endeed, the direction :math:`\dot q = T v` is tangential to
the sphere.

Posion formulation
~~~~~~~~~~~~~~~~~~

It consists in writting a position formulation on the system:

.. math::

   \label{NE_Dyn1}
   \begin{array}{l}
   h(q) = \left(\begin{array}{l}
     HI(q)\\HE(q)
   \end{array}\right)
   \end{array}

Approach using q
^^^^^^^^^^^^^^^^

| We are looking for :math:`q_1` from :math:`q_0`:

.. math:: q_1=q_0+\nabla HI \Lambda _I + \nabla HE \Lambda _E

| Assume that :math:`h(q_0)` doesn’t satisfy the constraints, ie
:math:` HI(q_0) \ngeq 0` or :math:` HE(q_0) \neq 0)`. Linearize
:math:`h` leads to:

.. math:: 0 \leq HI(q_0) + \nabla ^t HI (\nabla HI  \Lambda _I + \nabla HE \Lambda _E) \bot \Lambda _I \geq 0

.. math:: 0=HE(q_0) + \nabla ^t HE (\nabla HI  \Lambda _I + \nabla HE \Lambda _E)

The getting system could be written has a MLCP:

.. math:: C \ni h(q_0)+ \nabla ^t h(\nabla h\Lambda) , \Lambda \in C^*

In the case of a quaternion :math:`Q` for the rotation representation,
it is noteworthy that this system doesn’t deal with the constraints
:math:`\mid Q \mid = 1`. Thus, the direction :math:`(q_1,q_0)` can be
normal to this constraint, in that case this approach doesn’t work. (It
happens in practice) The solution that consists in normaliaed q after
this formulation is not convenient because, it could be incompatible
with :math:`\mid Q \mid = 1`. A better approach is to add this
constraint.

The constraint :math:`\mid Q \mid = 1` in the system HE:

.. math::

   \tilde {HE}(q)= \left(\begin{array}{l}
       HE(q)\\
       \mid Q \mid -1
   \end{array}\right)

The formulation described above can be done.

Approach using V
^^^^^^^^^^^^^^^^

It consists in building the OSNSP using :math:`CT` instead of :math:`C`.

.. math::

   \label{NE_projV}
   h(q_1) = h(q_0)+\nabla ^t H \delta q

ie:

.. math::

   \label{NE_projV}
   h(q_1) = h(q_0)+\nabla ^t H T V

We are looking for :math:`q_1` such that:

.. math:: q_1-q_0 = \nabla H \Lambda

We have

.. math:: \delta q=TV, \qquad ^tT\delta q=^tTTV, \qquad(^tTT)^{-1}~^tT\delta q=V

ie

.. math:: h(q_1)=h(q_0)+^t\nabla _q hT(^tTT)^{-1}~^tT\nabla _q h\Lambda

With :math:`C=^t\nabla _q h` leading to the prolem:

.. math:: K \ni h(q_0)+ CT~ (^tTT)^{-1}~ ^t(CT)   \Lambda \in K^*

Simulation of a Cam Follower System
===================================

| **M**\ ain Contributors: *Mario di Bernardo, Gustavo Osorio, Stefania
Santini*
| *University of Naples Federico II, Italy.*

The free body dynamics can be described by a linear second order system.
An external input is considered acting directly on the follower. This
input is a non linear forcing component coming from the valve. The
follower motion is constrained to a phase space region bounded by the
cam position. The non conservative Newton restitution law is used for
the computation of the post impact velocity. The cam is assumed to be
massive therefore only rotational displacement is allowed. Under these
assumptions, the free body dynamics of the follower can be described by

.. math::

   \begin{aligned}
     \label{eq:sols}
     \mu\frac{d^2u(t)}{dt^2}+\zeta\frac{du(t)}{dt}+\kappa
     u(t)=f_{v}(t),  \; \text{\hspace{5mm} \text{if} \hspace{3mm}$u(t) > c(t)$.}\end{aligned}

where :math:`\mu`, :math:`\zeta` and :math:`\kappa` are constant
parameters for the follower mass, friction viscous damping and spring
stiffness respectively. The state of the follower is given by the
position :math:`u(t)` and velocity :math:`v(t)={\frac{du}{dt}}`. The
external forcing is given by :math:`f_v(t)`. The cam angular position
determines :math:`c(t)` that defines the holonomic (i.e. constraint only
on the position) rheonomic (i.e. time varying) constraint. The dynamic
behavior when impacts occurs (i.e. :math:`u(t) = c(t)`) is modelled via
Newton’s impact law that in this case is given by

.. math::

   \begin{aligned}
     \label{eq:il}
     v(t^+)=
     \frac{dc}{dt}-r\left(v(t^-)-\frac{dc}{dt}\right)=(1+r)\frac{dc}{dt}-rv(t^-), \; \text{ \text{if}\hspace{3mm}$
   u(t)=c(t)$.}\end{aligned}

where :math:`v(t^+)` and :math:`v(t^-)` are the post and pre impact
velocities respectively, :math:`\frac{dc}{dt}` is the velocity vector of
the cam at the contact point with the follower, and :math:`r \in [0,1]`
is the restitution coefficient to model from plastic to elastic impacts.
In Figure [Fig:cam-shaft] is presented the schematic diagram of the
physical cam-follower system. In Figure [Fig:cam-shaft].a for
:math:`t=0`, [Fig:cam-shaft].b for :math:`t=\beta`, and
[Fig:cam-shaft].c the profile of the constraint position
:math:`\delta c(t)`, velocity :math:`\frac{dc}{dt}(t)` and acceleration
:math:`\frac{d^2c}{dt^2}(t)`. It is possible to visualize the follower
displacement as a function of the cam position. It is also important to
notice that different types of cams and followers profiles are used in
practical applications.

(90,80)(0,0) (0,0) (0,52) (2,34.5) (26,46) (15,46) (9,62) (32,62)
(38,40) (66,28) (2.5,6) (52.4,6) (18.5,-1) (69.4,-1)

(90,80)(-3,0) (0,0) (-3,75) (-3,49) (-3,25) (30,60) (58,60) (10,60)
(32,-1)

[Fig:cam-shaft]

The cam-follower as a Lagrangian NSDS.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to completely describe the cam-follower system as a
driven impact oscillator into the framework of *Lagrangian NSDS* using a
translation in space. Setting :math:`\hat u(t)=u(t)-c(t)` and
:math:`\hat
v(t)= v(t)-dc/dt`, then equations ([eq:sols]) and ([eq:il]) can be
expressed as (the argument :math:`t` will not be explicitly written)

.. math::

   \begin{aligned}
     \label{eq:trans}
     \mu\frac{d^2\hat u}{dt^2}+\zeta\frac{d\hat u}{dt}+\kappa
     \hat u=f_{v}-\left(\mu\frac{d^2c}{dt^2}+\zeta\frac{dc}{dt}+\kappa
     c\right)&\equiv &\hat f,  \; \text{\hspace{6.5mm} \text{if} \hspace{3mm}$\hat u >
    0$.}\\
   \hat v^+&=&-r \hat v^- , \; \text{ \text{if}\hspace{3mm}$\hat
   u=0$.}\end{aligned}

Using the framework presented in [2] we have that the equation of motion
of a Lagrangian system may be stated as follows :

.. math::

   \begin{aligned}
     \label{eq:lag1}
     M(q)\ddot q + Q(q,\dot q) + F(\dot q, q , t) = F_{ext}(t) + R\end{aligned}

From the ([eq:trans]) we can derive all of the terms which define a
Lagrangian NSDS. In our case the model is completely linear:

.. math::

   \begin{aligned}
     \nonumber
     q&=& \left[\begin{array}{c}  \hat u  \end{array}\right]    \\
     \nonumber
     M(q)&=&  \left[\begin{array}{c} \mu  \end{array}\right] \\
     \label{eq:lag2}
     Q(q,\dot q )& = &\left[\begin{array}{c} 0  \end{array}\right]  \\
     \nonumber
     F(q, \dot q ) &=&  \left[\begin{array}{c} \zeta \end{array}\right] \dot q +  \left[\begin{array}{c} \kappa  \end{array}\right] q\\
     \nonumber
     F_{ext}& = & \left[\begin{array}{c} \hat f \end{array}\right]\end{aligned}

The unilateral constraint requires that:

.. math::

   \begin{aligned}
   \label{eq:constr} \nonumber
    \hat u \geq 0\end{aligned}

so we can obtain

.. math::

   \begin{aligned}
   y &= & H^T q + b \\
   \nonumber H^T &=&\left[\begin{array}{c} 1 \end{array}\right]\\
   \nonumber b&=&0\end{aligned}

In the same way, the reaction force due to the constraint is written as
follows:

.. math::

   \begin{aligned}
   \nonumber R=H \lambda, \hspace{1cm}  \text{with }
   H=\left[\begin{array}{c} 1
   \end{array}\right]\end{aligned}

The unilataral contact law may be formulated as follow:

.. math::

   \begin{aligned}
     \label{eq:119}
     0 \leq y \perp \lambda\geq 0\end{aligned}

and the Newton’s impact law:

.. math::

   \begin{aligned}
     \label{eq:120}
   \text{If } y=0, \dot{y}^+ =-r\dot{y}^-\end{aligned}

Implementation in the platform
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the simulation of the cam follower system follow the steps

#. Move to the working directory ``sample/CamFollower``

   ``$cd sample/CamFollower ``

#. Clean the directory form binary files using the ``siconos`` command

   ``$siconos -c ``

#. Compile the file ``CamFollowerNoXml.cpp`` in the sample folder (*See*
   the code at the end of the section)

   ``$siconos CamFollowerNoXml.cpp``

#. Change the simulation parameters (*i.e.* Follower initial position
   and velocity, cam initial angle, simulations time, cam rotational
   speed in rpm, etc.) in the file ``CamFollowerNoXml.cpp``.

Next we present the sample code for the ``CamFollowerNoXml.cpp`` file:

=

=

=

=

| 
|  int main(int argc, char\* argv[]) **{**
|  **{**

=

=

=

=

=

| 
|  *// ======== Creation of the model =============
 *// User-defined main parameters
**

| double rpm=358;
|  double phi\_0=0;

unsigned int dsNumber = 1; *// the Follower and the ground
 unsigned int nDof = 1; *// degrees of freedom for the Follower
 double t0 = 0; *// initial computation time
 double T = 5; *// final computation time
 double h = 0.0001; *// time step
 int Kplot;
 Kplot=(int)(Tplot/h);
*****

double position\_init = 0.4; *// initial position for lowest bead.
 double velocity\_init = 0.4; *// initial velocity for lowest bead.
 *// ======= Dynamical systems =========
 vector<DynamicalSystem \*> vectorDS; // the list of DS
 vectorDS.resize(dsNumber,NULL);
 SiconosMatrix \*Mass, \*K, \*C; // mass/rigidity/viscosity
 Mass = new SiconosMatrix(nDof,nDof);
 (\*Mass)(0,0) = 1.221;
 K = new SiconosMatrix(nDof,nDof);
 (\*K)(0,0) = 1430.8;
 C = new SiconosMatrix(nDof,nDof);
 (\*C)(0,0) = 0;
 // Initial positions and velocities
 vector<SimpleVector \*> position\_0;
 vector<SimpleVector \*> velocity\_0;
 position\_0.resize(dsNumber,NULL);
 velocity\_0.resize(dsNumber,NULL);
 position\_0[0] = new SimpleVector(nDof);
 velocity\_0[0] = new SimpleVector(nDof);
 (\*(position\_0[0]))(0) = position\_init;
 (\*(velocity\_0[0]))(0) = velocity\_init;
 vectorDS[0] =
 new
LagrangianLinearTIDS(0,nDof,\*(position\_0[0]),\*(velocity\_0[0]),\*Mass,\*K,\*C);
 static\_cast<LagrangianDS\*>(vectorDS[0])
->setComputeFExtFunction(“FollowerPlugin.so”, “FollowerFExt”);
 // Example to set a list of parameters in FExt function.
 // 1 - Create a simple vector that contains the required parameters.
 // Here we set two parameters, the DS number.
 SimpleVector \* param = new SimpleVector(2);
 (\*param)(0)=rpm;
 (\*param)(1)=phi\_0;
 // 2 - Assign this param to the function FExt

static\_cast<LagrangianDS\*>(vectorDS[0])->setParametersListPtr(param,2);
 // 2 corresponds to the position of FExt in the stl vector of possible
parameters.
 // 0 is mass, 1 FInt.
 // Now the cam rotational velocity in rpms will be available in FExt
plugin.
 // ===== Interactions =====
 vector<Interaction\*> interactionVector;
 interactionVector.resize(1,NULL);
 vector<DynamicalSystem\*> \*dsConcerned =
 new vector<DynamicalSystem\*>(dsNumber);
 // ===== Non Smooth Law =====
 double e = 0.8;
***

| // Interaction Follower-floor
|  SiconosMatrix \*H = new SiconosMatrix(1,nDof);
|  (\*H)(0,0) = 1.0;
|  NonSmoothLaw \* nslaw = new NewtonImpactLawNSL(e);
|  Relation \* relation = new LagrangianLinearR(\*H);
|  (\*dsConcerned)[0] = vectorDS[0];

| interactionVector[0] = new Interaction(“Follower-Ground”,0,1,
dsConcerned);
|  interactionVector[0]->setRelationPtr(relation);
|  interactionVector[0]->setNonSmoothLawPtr(nslaw);

| // ===== Interactions =====
|  // ===== NonSmoothDynamicalSystem =====

| bool isBVP =0;
|  NonSmoothDynamicalSystem \* nsds =
|  new NonSmoothDynamicalSystem(isBVP);
| // Set DS of this NonSmoothDynamicalSystem
|  nsds->setDynamicalSystems(vectorDS);
|  // Set interactions of the NonSmoothDynamicalSystem
|  nsds->setInteractions(interactionVector);
|  // ===== Model =====
|  Model \* Follower = new Model(t0,T);
|  // set NonSmoothDynamicalSystem of this model
|  Follower->setNonSmoothDynamicalSystemPtr(nsds);
|  // ====== Strategy ======
|  double theta = 0.5; // theta for Moreau integrator
|  string solverName = “QP” ;
|  Strategy\* S = new TimeStepping(Follower);
|  // – Time discretisation –
|  TimeDiscretisation \* t = new TimeDiscretisation(h,S);
|  // – OneStepIntegrators –
|  vector<OneStepIntegrator \*> vOSI;
|  vOSI.resize(dsNumber,NULL);
|  vOSI[0] = new Moreau(t,vectorDS[0],theta);
|  S->setOneStepIntegrators(vOSI);
|  // – OneStepNsProblem –
|  OneStepNSProblem \* osnspb = new LCP(S,solverName,101,
0.0001,“max”,0.6);
|  S->setOneStepNSProblemPtr(osnspb); // set OneStepNSProblem of the
strategy
|  cout << “=== End of model loading === ” << endl;
|  // ==== End of model definition======
|  // ========= Computation============
|  // — Strategy initialization —
|  S->initialize();
|  cout <<“End of strategy initialisation” << endl;

| int k = t->getK(); // Current step
|  int N = t->getNSteps(); // Number of time steps
|  // — Get the values to be plotted —
|  // -> saved in a matrix dataPlot
|  unsigned int outputSize = 8;
|  SiconosMatrix DataPlot(Kplot+1,outputSize );
|  // For the initial time step:
|  // time
|  DataPlot(k,0) = k\*t->getH();
|  DataPlot(k,1) = static\_cast<LagrangianDS\*>(vectorDS[0])->getQ()(0);
|  DataPlot(k,2) =
static\_cast<LagrangianDS\*>(vectorDS[0])->getVelocity()(0);
|  DataPlot(k,3) = (Follower->getNonSmoothDynamicalSystemPtr()->
|  getInteractionPtr(0)->getLambda(1))(0);
|  DataPlot(k,4) =
static\_cast<LagrangianDS\*>(vectorDS[0])->getFExt()(0);
|  // State of the Cam
|  double CamEqForce,CamPosition,CamVelocity,CamAcceleration;

| CamEqForce=
|  CamState(k\*t->getH(),rpm,CamPosition,CamVelocity,CamAcceleration);
|  // Position of the Cam
|  DataPlot(k, 5) = CamPosition;
|  // Velocity of the Cam
|  DataPlot(k, 6) = CamVelocity;
|  // Acceleration of the Cam
|  DataPlot(k, 7) =
|  CamPosition+static\_cast<LagrangianDS\*>(vectorDS[0])->getQ()(0);
|  // — Time loop —
|  cout << “Start computation ... ” << endl;
|  while(k < N)
|  **{ **

=

=

=

=

=

| 
|  // — Get values to be plotted —
|  DataPlot(k,0) = k\*t->getH();
|  DataPlot(k,1) =
|  static\_cast<LagrangianDS\*>(vectorDS[0])->getQ()(0);
|  DataPlot(k,2) =
|  static\_cast<LagrangianDS\*>(vectorDS[0])->getVelocity()(0);
|  DataPlot(k,3) =
|  (Follower->getNonSmoothDynamicalSystemPtr()->
|  getInteractionPtr(0)->getLambda(1))(0);
|  DataPlot(k,4) =
static\_cast<LagrangianDS\*>(vectorDS[0])->getFExt()(0);

| CamEqForce=
|  CamState(k\*t->getH(),rpm,CamPosition,CamVelocity,CamAcceleration);
|  DataPlot(k, 5) = CamPosition;
|  DataPlot(k, 6) = CamVelocity;
|  DataPlot(k, 7) = CamPosition+
|  static\_cast<LagrangianDS\*>(vectorDS[0])->getQ()(0);

| // transfer of state i+1 into state i and time incrementation
|  S->nextStep();

| // get current time step
|  k = t->getK();
|  // solve ...
|  S->computeFreeState();
|  S->computeOneStepNSProblem();
|  // update
|  S->update();
|  **} **
|  // — Output files —
|  DataPlot.rawWrite(“result.dat”, “ascii”);

| // — Free memory —
|  delete osnspb;
|  delete vOSI[0];
|  delete t;
|  delete S;
|  delete Follower;
|  delete nsds;
|  delete interactionVector[0];
|  delete relation;
|  delete nslaw;
|  delete H;
|  delete dsConcerned;
|  delete vectorDS[0];
|  delete position\_0[0];
|  delete velocity\_0[0];
|  delete C;
|  delete K;
|  delete Mass;

| 
|  **}**

Simulation
~~~~~~~~~~

We have perform the simulation of the cam follower system for different
values of the cam rotational speed with the SICONOS software package
using a time-stepping numerical scheme with step size
(:math:`h=1e^{-4}`) and an event-driven scheme with minimum step size
(:math:`h_{min}=1e^{-12}`). Fig. [Fig:time:sub:`c`\ omparison] and
[Fig:state:sub:`c`\ omparison] show the time simulations for different
values of the cam rotational speed and Fig.
[Fig:attractor:sub:`c`\ omparison] show the chaotic attractor at
:math:`rpm=660` for impact and stroboscopic Poincarè sections.

(60,60)(0,-7) (0,0) (35,-4)

(60,60)(15,-7) (0,0) (35,-4)

(60,60)(-40,-2) (0,0) (35,-4)

[Fig:time:sub:`c`\ omparison]

(60,60)(0,-7) (0,0) (35,-4)

(60,60)(15,-7) (0,0) (35,-4)

(60,60)(0,-2) (0,0) (35,-4)

(60,60)(-17,-2) (0,0) (35,-4)

[Fig:state:sub:`c`\ omparison]

(60,60)(0,-7) (0,0) (35,-4)

(60,60)(15,-7) (0,0) (35,-4)

(60,60)(0,-2) (0,0) (35,-4)

(60,60)(-17,-2) (0,0) (35,-4)

[Fig:attractor:sub:`c`\ omparison]

Quartic Formulation
===================

Slidding ?
~~~~~~~~~~

It consists in finding :math:`\alpha >0` and
:math:`R \in \partial K_{\mu}` such that
:math:`-\alpha \left(\begin{array}{l} 0\\ R_T\end{array}\right)=MR+q`.
That is :

.. math::

   \label{eq_quartic1}
   \left[\begin{array}{c}
   M+ \left(\begin{array}{ccc} 0&0&0\\ 0&\alpha&0 \\ 0&0&\alpha \end{array}\right)
   \end{array}\right]R+q=0

:math:`R_T` is on a conic
^^^^^^^^^^^^^^^^^^^^^^^^^

The first line of the system [eq\ :sub:`q`\ uartic1] and the
:math:`R \in \partial K_{\mu}` is the intersection between a plan and a
cone in :math:`\mathbb{R}^3`, endeed:

.. math::

   \label{eq_quartic2}
   \begin{array}{l}
    \mu R_N =  \parallel R_T \parallel  \\
   \frac{M_{11}}{\mu} \parallel R_T \parallel = -q_1-M_{12}R_{T1}-M_{13}R_{T2}
   \end{array}

That is:

.. math::

   \label{eq_quartic2}
   \begin{array}{l}
   \mu^2 R_N^2 =  (R_{T1}^2 +R_{T1}^2)  \\
   \frac{M_{11}^2}{\mu^2} (R_{T1}^2 +R_{T1}^2)=(-q_1-M_{12}R_{T1}-M_{13}R_{T2})^2
   \end{array}

That means that :math:`R_T` is contained in a conic, focus and
directrice are:

.. math::

   \label{eq_quartic3}
   \begin{array}{l}
   \mathcal{D} : q_1+M_{12}R_{T1}+M_{13}R_{T2} =0  \\
   focus : \mathcal{O}\\
   \frac{M_{11}^2}{\mu^2}  Dist(\mathcal{O}, R_T) ^2=Dist(\mathcal{D},R_T)^2 (M_{12}^2+M_{13}^2)\\
   \frac{Dist(\mathcal{O}, R_T)}{Dist(\mathcal{D},R_T)}=\frac{\mu\sqrt{(M_{12}^2+M_{13}^2)}}{M_{11} }=e
   \end{array}

The parametric equation is:

.. math::

   \label{eq_quartic4}
   \begin{array}{l}
   R_{T1}=r cos(\theta )\\
   R_{T2}=r sin(\theta )\\
   r=\frac{p}{1+ecos(\theta - \phi)}
   \end{array}

With :math:`p` an simple expression of :math:`M_{11},M_{12},M_{13}`, and
:math:`\phi` a constant angle between :math:`\mathcal{D}` and
:math:`(O,R_{T1})`

The two last line of the system [eq\ :sub:`q`\ uartic1]
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::

   \label{eq_quartic5}
   \frac{\parallel R_T \parallel}{\mu} \tilde M_{1.} +\left(\tilde M+\left(\begin{array}{cc} \alpha&0 \\ 0&\alpha \end{array}\right)\right)R_T+\tilde q=0

:math:`\tilde M` is symetric, so it exists a unitary matrix :math:`V`
such that
:math:`V \tilde M V^T = \left(\begin{array}{cc} d_1&0 \\ 0&d_2 \end{array}\right)`.
One can get:

.. math::

   \label{eq_quartic6}
   \frac{\parallel R_T \parallel}{\mu} V \tilde M_{1.} +V \left(\tilde M+\left(\begin{array}{cc} \alpha&0 \\ 0&\alpha \end{array}\right)\right)V^TVR_T+V\tilde q=0

Rename:

.. math::

   \label{eq_quartic7}
     \frac{\parallel \bar R_T \parallel}{\mu} \bar M_{1.} +\left(\begin{array}{cc} d_1+\alpha&0 \\ 0&d_2+\alpha \end{array}\right)\overline R_T+\bar q=0

In the plan, either :math:`V` is a rotation or a symetrie. So
:math:` \bar R_T=VR_T` is a conic with the same focus and a rotated
directrice, it means that it exists :math:`\phi_1` such that :

.. math::

   \label{eq_quartic8}
   \begin{array}{l}
   \bar R_{T1}=r cos(\theta )\\
   \bar R_{T2}=r sin(\theta )\\
   r=\frac{p}{1+ecos(\theta - \phi_1)}
   \end{array}

The equation [eq\ :sub:`q`\ uartic7] is :

.. math::

   \label{eq_quartic9}
   \begin{array}{l}
     (d_1+\alpha)\bar R_{T1}=-\bar q_1+a_1 \parallel R_T \parallel\\
   (d_2+\alpha)\bar R_{T2}=-\bar q_2+a_2 \parallel R_T \parallel
   \end{array}

The case (:math:`\bar R_{T1} = 0` or :math:`\bar R_{T2} = 0`) has to be
examine. We try to eliminate :math:`alpha`:

.. math::

   \label{eq_quartic10}
     \begin{array}{l}
       d_1 \bar R_{T1} \bar R_{T2}+\alpha \bar R_{T1} \bar R_{T2} =-\bar q_1\bar R_{T2}+a_1 \bar R_{T2} \parallel R_T \parallel\\
   d_2 \bar R_{T1} \bar R_{T2}+\alpha \bar R_{T1} \bar R_{T2} =-\bar q_2\bar R_{T1}+a_2 \bar R_{T1} \parallel R_T \parallel
   \end{array}

that leads to:

.. math::

   \label{eq_quartic10}
     (d_1-d_2) \bar R_{T1} \bar R_{T2}=-\bar q_1\bar R_{T2}+\bar q_2\bar R_{T1}+(a_1 \bar R_{T2}-a_2 \bar R_{T1}) \parallel R_T \parallel\\

The parametric expression of :math:`\bar R_T` leads to:

.. math::

   \label{eq_quartic11}
   \begin{array}{l}
     (d_1-d_2)r^2cos(\theta )sin(\theta )=-\bar q_1rsin(\theta )+\bar q_2rcos(\theta )+r(a_1 rsin(\theta )-a_2 rcos(\theta )) \\
     \textrm{ie:}(d_1-d_2)rcos(\theta )sin(\theta )=-\bar q_1sin(\theta )+\bar q_2cos(\theta )+r(a_1 sin(\theta )-a_2 cos(\theta ))\\
     \end{array}

with the expression of r:

.. math::

   \label{eq_quartic12}
   \begin{array}{l}
   (d_1-d_2)\frac{p}{1+ecos(\theta - \phi_1)}cos(\theta )sin(\theta )=\\-\bar q_1sin(\theta )+\bar q_2cos(\theta )+\frac{p}{1+ecos(\theta - \phi_1)}(a_1  sin(\theta )-a_2 cos(\theta ))\\\\
   \textrm{ie:}(d_1-d_2)pcos(\theta )sin(\theta )=\\(1+ecos(\theta - \phi_1))(-\bar q_1sin(\theta )+\bar q_2cos(\theta ))+p(a_1  sin(\theta )-a_2 cos(\theta ))\\\\
   \textrm{ie:}(d_1-d_2)pcos(\theta )sin(\theta )=\\(1+e(cos(\theta)cos(\phi_1)+sin(\theta)sin(\phi_1)))(-\bar q_1sin(\theta )+\bar q_2cos(\theta ))+p(a_1  sin(\theta )-a_2 cos(\theta ))\\\\
   \textrm{ie:}(d_1-d_2)pcos(\theta )sin(\theta )+\\(1+ecos(\theta)cos(\phi_1)+esin(\theta)sin(\phi_1))(\bar q_1sin(\theta )-\bar q_2cos(\theta ))+p(-a_1  sin(\theta )+a_2 cos(\theta ))=0
    \end{array}

rename :

.. math::

   \label{eq_quartic13}
   \begin{array}{l}
   Acos(\theta )^2+Bsin(\theta)^2+Csin(\theta )cos(\theta )+Dsin(\theta )+Ecos(\theta )=0
    \end{array}

with

.. math::

   \label{eq_quartic12}
   \begin{array}{l}
   A=- e\bar q_2cos(\phi_1)\\
   B=e \bar q_1sin(\phi_1)\\
   C=(d_1-d_2)p+ecos(\phi_1)\bar q_1-esin(\phi_1)\bar q_2\\
   D=\bar q_1-pa_1\\
   E=-\bar q_2+pa_2\\
   \end{array}

rename : Using the following set of unknown :

.. math::

   \label{eq_quartic14}
   \begin{array}{l}
   t=tan(\theta /2)\\
   sin(\theta )=\frac{2t}{1+t^2}\\
   cos(\theta )=\frac{1-t^2}{1+t^2}
    \end{array}

leads to:

.. math::

   \label{eq_quartic13}
   \begin{array}{l}
     A\frac{(1-t^2)^2}{1+t^2} +B\frac{4t^2}{1+t^2}+ C\frac{2t(1-t^2)}{1+t^2}+D2t+E(1-t^2)=0\\
   \textrm{ie:}A(1-t^2)^2 + 4Bt^2+C2t(1-t^2)+2Dt(1+t^2)+E(1-t^2)(1+t^2)=0\\\\
   \textrm{ie:}P_4=A-E\qquad P_3=-2C+2D \qquad P_2=4B-2A \qquad P_1=2C+2D \qquad P_0=A+E
    \end{array}

Finally, we get 4 possible values for :math:`R_T`, checking the sign of
:math:`\alpha` and :math:`R_N` selects the solutions.

case :math:`R_{T12}=0`
^^^^^^^^^^^^^^^^^^^^^^

From [eq\ :sub:`q`\ uartic9], :math:`R_{T1}` leads to:

.. math::

   \label{eq_quartic14}
   \begin{array}{l}
     \parallel R_T \parallel=|\bar R_{T2}|=\frac{\bar q_1}{a_1}\\\\
     \bar R_T=\left(\begin{array}{c} 0 \\ \pm \frac{\bar q_1}{a_1} \end{array}\right)
    \end{array}

From [eq\ :sub:`q`\ uartic9], :math:`R_{T2}` leads to:

.. math::

   \label{eq_quartic14}
   \begin{array}{l}
     \parallel R_T \parallel=|\bar R_{T1}|=\frac{\bar q_2}{a_2}\\\\
     \bar R_T=\left(\begin{array}{c}  \pm \frac{\bar q_2}{a_2} \\ 0 \end{array}\right)
    \end{array}

From :math:`\bar R_T`, we have to check the coherence with the
equation [eq\ :sub:`q`\ uartic8]. If it is on the conic, we compute R,
and the sign condition of the equation [eq\ :sub:`q`\ uartic1] must be
check.

Alart–Curnier Formulation
=========================

Reduced formulation to local variables.
---------------------------------------

Formulation
~~~~~~~~~~~

Let us start with

.. math::

   \label{eq:AC-L7}
     \begin{array}{l}
     \varPhi_1(U,P) =  - U_{k+1}  + \widehat W P_{k+1}  + U_{\mathrm{free}}\\ \\
     \varPhi_2(U,P) =  P_{\n} - \operatorname{proj}_{\nbR^{a}_+} (P_{\n} - \rho_{\n}\circ (U_{\n} +e \circ  U_{\n,k}) ) \\ \\
     \varPhi_3(U,P) =  P_{\t} - \operatorname{proj}_{\widehat {\bf D}(P_{\n},U_{\n})} (P_{{\t}} - \rho_{\t}\circ \,U_{\t} )
   \end{array}

where the modified friction disk for a contact :math:`\alpha` is

.. math::

   \label{eq:AC-L3}
     \widehat {\bf D}^\alpha(P^\alpha_{\n,k+1},U_{\n,k+1}^{\alpha}) = {\bf D}(\mu(\operatorname{proj}_{\nbR_+} (P^\alpha_{\n,k+1} - \rho^\alpha_{\n}\,(U_{\n,k+1}^{\alpha}+e^\alpha U_{\n,k}^{\alpha}) )).

Structure of the Jacobians
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us denote the one element of the generalized Jacobian by
:math:` H(U,P) \in \partial \Phi(U,P)` which has the structure

.. math::

   \label{eq:AC-L6}
      H(U,P) = 
      \left[\begin{array}{cccc}
          - I & 0 &  \widehat W_{\n\n} & \widehat W_{\n\t} \\ \\
          0  & -I  &  \widehat W_{\t\n} & \widehat W_{\t\t} \\ \\
          \partial_{U_{\n}} \Phi_2(U,P) & 0 &   \partial_{P_{\n}} \Phi_2(U,P) & 0 \\ \\
          \partial_{U_{\n}} \Phi_3(U,P) &  \partial_{U_{\t}} \Phi_3(U,P) &  \partial_{P_{\n}} \Phi_3(U,P)  & \partial_{P_{\t}} \Phi_3(U,P)
      \end{array}\right]

Computation of the gradients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us consider the single contact case.

Computation of the gradients of :math:`\Phi_2`
''''''''''''''''''''''''''''''''''''''''''''''

.. math::

   \label{eq:AC-T1}
     \begin{array}{l}
     \varPhi_2(U,P) =  P_{\n} - \operatorname{proj}_{\nbR^{a}_+} (P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) ) \\ \\
   \end{array}

-  **If** :math:`P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) \geq 0 `, we
   get

   .. math::

      \label{eq:AC-T2}
          \begin{array}{l}
            \varPhi_2(U,P) =  + \rho_{\n} (U_{\n} +e  U_{\n,k})
          \end{array}

   and

   .. math::

      \label{eq:AC-T3}
          \begin{array}{l}
           \partial_{U_{\n}} \varPhi_2(U,P) =  + \rho_{\n} \\ \\
           \partial_{P_{\n}} \varPhi_2(U,P) =  0 \\ \\ 
          \end{array}

-  **If** :math:`P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k})  < 0 `, we get

   .. math::

      \label{eq:AC-T4}
          \begin{array}{l}
            \varPhi_2(U,P) =  P_{\n}
          \end{array}

   and

   .. math::

      \label{eq:AC-T5}
          \begin{array}{l}
           \partial_{U_{\n}} \varPhi_2(U,P) =  0 \\ \\
           \partial_{P_{\n}} \varPhi_2(U,P) =  1 \\ \\ 
          \end{array}

Computation of the gradients of :math:`\Phi_3`
''''''''''''''''''''''''''''''''''''''''''''''

.. math::

   \label{eq:AC-TT1}
     \begin{array}{l}
     \varPhi_3(U,P) =  P_{\t} - \operatorname{proj}_{\widehat {\bf D}(P_{\n},U_{\n})} (P_{\t} - \rho_{\t} U_{\t} ) \\ \\
   \end{array}

-  **If**
   :math:`\|P_{\t} - \rho_{\t} U_{\t}\| \leq \mu \max (0 ,P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) ) `
   , we get

   .. math::

      \label{eq:AC-TT2}
        \begin{array}{l}
        \varPhi_3(U,P) =  + \rho_{\t} U_{\t} 
      \end{array}

   and

   .. math::

      \label{eq:AC-TT3}
          \begin{array}{l}
           \partial_{U_{\n}} \varPhi_3(U,P) =  0 \\ \\
           \partial_{P_{\n}} \varPhi_3(U,P) =  0 \\ \\ 
           \partial_{U_{\t}} \varPhi_3(U,P) =  + \rho_{\t} \\ \\
           \partial_{P_{\t}} \varPhi_3(U,P) =  0 \\ \\ 
          \end{array}

-  **If**
   :math:`\|P_{\t} - \rho_{\t} U_{\t}\| > \mu \max (0 ,P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) ) `
   , we get

   .. math::

      \label{eq:AC-TT4}
        \begin{array}{l}
        \varPhi_3(U,P) =  P_{\t} - \mu \max(0,P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) )  {\displaystyle \frac{P_{\t} - \rho_{\t} U_{\t} }{ \| P_{\t} - \rho_{\t} U_{\t}\| }}
      \end{array}

   -  **If** :math:`P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) \leq 0`, we
      get

      .. math::

         \label{eq:AC-TT5}
           \begin{array}{l}
           \varPhi_3(U,P) =   P_{\t}
         \end{array}

      and

      .. math::

         \label{eq:AC-TT6}
            \begin{array}{l}
              \partial_{U_{\n}} \varPhi_3(U,P) =  0 \\ \\
              \partial_{P_{\n}} \varPhi_3(U,P) =  0 \\ \\ 
              \partial_{U_{\t}} \varPhi_3(U,P) =  0 \\ \\
              \partial_{P_{\t}} \varPhi_3(U,P) =  I_2 \\ \\ 
            \end{array}

   -  **If** :math:`P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) > 0`, we
      get

      .. math::

         \label{eq:AC-TT7}
           \begin{array}{l}
           \varPhi_3(U,P) =  P_{\t} - \mu (P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) )  {\displaystyle \frac{P_{\t} - \rho_{\t} U_{\t} }{ \| P_{\t} - \rho_{\t} U_{\t}\| }}
         \end{array}

      and

      .. math::

         \label{eq:AC-TT8}
            \begin{array}{l}
              \partial_{U_{\n}} \varPhi_3(U,P) =  \mu \rho_{\n}  {\displaystyle \frac{P_{\t} - \rho_{\t} U_{\t} }{ \| P_{\t} - \rho_{\t} U_{\t}\| }}\text{{\bf WARNING} case was not taken into account}\\ \\
              \partial_{P_{\n}} \varPhi_3(U,P) =  -\mu  {\displaystyle \frac{P_{\t} - \rho_{\t} U_{\t} }{ \| P_{\t} - \rho_{\t} U_{\t}\| }} \\ \\ 
              \partial_{U_{\t}} \varPhi_3(U,P) =  \mu\rho_{\t}(P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) ) \Gamma(P_{\t} - \rho_{\t} U_{\t})  \\ \\
              \partial_{P_{\t}} \varPhi_3(U,P) =  I_2-\mu(P_{\n} - \rho_{\n} (U_{\n} +e  U_{\n,k}) ) \Gamma(P_{\t} - \rho_{\t} U_{\t})  \\ \\ 
            \end{array}

Rearranging the cases
~~~~~~~~~~~~~~~~~~~~~

**TO BE COMPLETED**

Formulation with global variables.
----------------------------------

Formulation
~~~~~~~~~~~

Let us start with

.. math::

   \label{eq:GAC-L1}
     \begin{array}{l}
     \Psi_{1}^{a}(v,U,P) =  - \widehat M v_{k+1}  +  H P_{k+1}  + q \\ \\
     \Psi_{1}^{b}(v,U,P) =  - U_{k+1}  + H^\top v _{k+1}  + b \\ \\
     \Psi_2(v,U,P) =  P_{\n} - \operatorname{proj}_{\nbR^{a}_+} (P_{\n} - \rho_{\n}\circ (U_{\n} +e \circ  U_{\n,k}) ) \\ \\
     \Psi_3(v,U,P) =  P_{\t} - \operatorname{proj}_{\widehat {\bf D}(P_{\n},U_{\n})} (P_{{\t}} - \rho_{\t}\circ \,U_{\t} )
   \end{array}

where the modified friction disk for a contact :math:`\alpha` is

.. math::

   \label{eq:GAC-L2}
     \widehat {\bf D}^\alpha(P^\alpha_{\n,k+1},U_{\n,k+1}^{\alpha}) = {\bf D}(\mu(\operatorname{proj}_{\nbR_+} (P^\alpha_{\n,k+1} - \rho^\alpha_{\n}\,(U_{\n,k+1}^{\alpha}+e^\alpha U_{\n,k}^{\alpha}) )).

Structure of the Jacobians
~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us denote the one element of the generalized Jacobian by
:math:` H(v,U,P) \in \partial \Psi(s,U,P)` which has the structure

.. math::

   \label{eq:GAC-L3}
      H(v,U,P) = 
      \left[\begin{array}{ccccc}
          - \widehat M & 0 & 0 & H_{\n} & H_{\t} \\ \\
           H_{\n}^\top &  - I & 0 & 0 &0 \\ \\
           H_{\t}^\top &  0  & -I & 0 &0 \\ \\
           0 & \partial_{U_{\n}} \Psi_2(v,U,P) & 0 &   \partial_{P_{\n}} \Psi_2(v,U,P) & 0 \\ \\
           0 & \partial_{U_{\n}} \Psi_3(v,U,P) &  \partial_{U_{\t}} \Psi_3(v,U,P) &  \partial_{P_{\n}} \Psi_3(v,U,P)  & \partial_{P_{\t}} \Psi_3(v,U,P)
      \end{array}\right]

We clearly have

.. math::

   \label{eq:equivalentJacobian}
     \begin{array}{lcl}
        \partial_{U} \Psi_2(v,U,P) &=& \partial_{U} \Phi_2(U,P) \\ 
        \partial_{P} \Psi_2(v,U,P) &=& \partial_{P} \Phi_2(U,P) \\     
        \partial_{U} \Psi_3(v,U,P) &=& \partial_{U} \Phi_3(U,P) \\ 
        \partial_{P} \Psi_3(v,U,P) &=& \partial_{P} \Phi_3(U,P) \\
     \end{array}

and we get

.. math::

   \label{eq:GAC-L4}
      H(v,U,P) = 
      \left[\begin{array}{ccccc}
          - \widehat M & 0 & 0 & H_{\n} & H_{\t} \\ \\
           H_{\n}^\top &  - I & 0 & 0 &0 \\ \\
           H_{\t}^\top &  0  & -I & 0 &0 \\ \\
           0 & \partial_{U_{\n}} \Phi_2(U,P) & 0 &   \partial_{P_{\n}} \Phi_2(U,P) & 0 \\ \\
           0 & \partial_{U_{\n}} \Phi_3(U,P) &  \partial_{U_{\t}} \Phi_3(U,P) &  \partial_{P_{\n}} \Phi_3(U,P)  & \partial_{P_{\t}} \Phi_3(U,P)
      \end{array}\right]

Simplification ?
~~~~~~~~~~~~~~~~

Since the second line :math:`\Psi_1^b` is linear, we should be able to
derive a reduced Jacobian using the chain rule. Let us define
:math:`\widetilde \Psi`

.. math::

   \label{eq:chainrule}
     \widetilde \Psi(v,P)  = \Psi(v,H^\top v +b,P)

.. math::

   \label{eq:GAC-L5}
     \begin{array}{l}
     \widetilde \Psi_{1}(v,P) =  - \widehat M v_{k+1}  +  H P_{k+1}  + q \\ \\
     \widetilde \Psi_2(v,P) =  P_{\n} - \operatorname{proj}_{\nbR^{a}_+} (P_{\n} - \rho_{\n}\circ (H^\top_{\n}v+b_{\n} +e \circ  U_{\n,k}) ) \\ \\
     \widetilde \Psi_3(v,P) =  P_{\t} - \operatorname{proj}_{\widehat {\bf D}(P_{\n},U_{\n})} (P_{{\t}} - \rho_{\t}\circ \,(H^\top_\t v + b_\t) )
   \end{array}

Chain rule
''''''''''

.. math::

   \label{eq:chainrule1}
     \begin{array}{lcl}
     \partial_v \widetilde \Psi_{2,3}(v,P) &=&  \partial_v \Psi_{2,3}(v,H^\top v +b,P)  \\ \\
     &=& H_{\n}^\top \partial_{U_\n} \Phi_{2,3}(H^\top v + b,P) + H_{\t}^\top \partial_{U_\t} \Phi_{2,3}(H^\top v + b,P)  
   \end{array}

.. math::

   \label{eq:GAC-L6}
      H(v,P) = 
      \left[\begin{array}{ccc}
          - \widehat M &   H_{\n} & H_{\t} \\ \\
          H_{\n}^\top \partial_{U_\n} \Phi_{2}(H^\top v + b,P) &   \partial_{P_{\n}} \Phi_2(H^\top v + b,P) & 0 \\ \\
          \begin{array}{c}
            H_{\n}^\top \partial_{U_\n} \Phi_{3}(H^\top v + b,P) \\
            \quad \quad + H_{\t}^\top \partial_{U_\t} \Phi_{3}(H^\top v + b,P)\\
        \end{array}
        &  \partial_{P_{\n}} \Phi_3(H^\top v + b,P)  & \partial_{P_{\t}} \Phi_3(H^\top v + b,P)
      \end{array}\right]

discussion
''''''''''

-  Formulae has to be checked carefully

-  I do not known if there an interest in the simplification. With
   sparse matrices, it is perhaps easier to deal with ([eq:GAC-L4])

.. |image| image:: ./DSClassDiagram.pdf
