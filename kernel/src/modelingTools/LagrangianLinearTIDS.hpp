/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*! \file LagrangianLinearTIDS.hpp

 */
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.hpp"


/** Lagrangian Linear Systems with time invariant coefficients - \f$M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$

    The class LagrangianLinearTIDS  allows to define  and compute a generic ndof-dimensional
    Lagrangian Linear Time Invariant Dynamical System of the form :

    \f$
    M \ddot q + C \dot q + K q =  F_{ext}(t,z) + p,
    \f$

    where
    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,
    the variable p contains the impulse and not the force.
    -  \f$ M \in  R^{ndof \times ndof} \f$ is the Mass matrix (access : mass() method).
    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix (access : K() method).
    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix (access : C() method).
    -  \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.

   The equation of motion is also shortly denoted as:
    \f$
    M(q,z) \dot v = F(v, q, t, z) + p
    \f$

    where
    - \f$F(v, q, t, z) \in R^{ndof} \f$ collects the total forces
    acting on the system, that is
    \f$ F(v, q, t, z) =  F_{ext}(t, z) -  Cv - Kq \f$.

    This vector is saved and may be accessed using forces() method.

    If required (e.g. for Event-Driven like simulation), reformulation as a first-order system is also available, and writes:

    - \f$ n= 2 ndof \f$
    - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
    - rhs given by:

    \rst
    .. math::

        rhs(x,t,z) = \left[\begin{array}{c}
        \\dot q  \\
        \ddot q = M^{-1}\left[F_{ext}(t, z) - C \\dot q - K q  + p \right]\\
        \end{array}\right]

    \endrst
    Its jacobian is:

    \rst
    .. math::

        \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
        0   & I \\
        -M^{-1}K & -M^{-1}C \\
        \end{array}\right]

   \endrst

   with the input due to the non smooth law:
  \rst

   .. math::

      r = \left[\begin{array}{c}0 \\ p \end{array}\right]

  \endrst
*/
class LagrangianLinearTIDS : public LagrangianDS
{

protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** stiffness matrix */
  SP::SiconosMatrix _K;

  /** damping matrix */
  SP::SiconosMatrix _C;

  /** default constructor */
  LagrangianLinearTIDS():LagrangianDS() {};

public:

  /*! @name public constructors */
  //@{

  /** constructor from initial state and all matrix operators.
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   *  \param K stiffness matrix
   *  \param C damping matrix
   */
  LagrangianLinearTIDS(SP::SiconosVector q0, SP::SiconosVector v0,
                       SP::SiconosMatrix M, SP::SiconosMatrix K, SP::SiconosMatrix C);

  /** constructor from initial state and mass matrix only. Leads to \f$ M\dot v = F_{ext}(t,z) + p\f$.
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param M mass matrix
   */
  LagrangianLinearTIDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix M):
    LagrangianDS(q0, v0, M){};

  /** destructor */
  ~LagrangianLinearTIDS(){};

  ///@}

  /*! @name Right-hand side computation */
  //@{

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param t time of initialization
   */
  void initRhs(double t) ;

  /** Compute \f$F(v,q,t,z)\f$
   *  \param time the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  void computeForces(double time, SP::SiconosVector q, SP::SiconosVector velocity);

  ///@}

  /*! @name Attributes access
    @{ */

  /** get a copy of the stiffness matrix
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getK() const
  {
    return *_K;
  }

  /** get stiffness matrix (pointer link)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix K() const
  {
    return _K;
  }

  /** set (copy) the value of the stiffness matrix
   *  \param K new stiffness matrix
   */
  void setK(const SiconosMatrix& K);

  /** set stiffness matrix (pointer link)
   *  \param newPtr pointer to the new Stiffness matrix
   */
  void setKPtr(SP::SiconosMatrix newPtr);

  /** get a copy of the damping matrix
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getC() const { return *_C; }

  /** get damping matrix (pointer link)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix C() const
  {
    return _C;
  }

  /** set (copy) the value of the damping matrix
   *  \param C new damping matrix
   */
  void setC(const SiconosMatrix& C);

  /** set damping matrix (pointer link)
   * \param newPtr pointer to the new damping matrix
   */
  void setCPtr(SP::SiconosMatrix newPtr) ;

  /** get \f$ \nabla_qF(v,q,t,z)\f$ (pointer  link)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqForces() const
  {
    return _K;
  }

  /** get \f$ \nabla_{\dot q}F(v,q,t,z)\f$ (pointer  link)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqvForces() const
  {
    return _C;
  }

  ///@}

  /*! @name Miscellaneous public methods */
  //@{

  /**\return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }

  /** print the data onto the screen
   */
  void display(bool brief =true) const;

  ///@}

  ACCEPT_STD_VISITORS();

};
#endif // LAGRANGIANTIDS_H
