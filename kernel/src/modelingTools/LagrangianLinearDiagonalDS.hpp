/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file LagrangianLinearDiagonalDS.hpp */
#ifndef LAGRANGIANLINEARDIAGONALDS_H
#define LAGRANGIANLINEARDIAGONALDS_H

#include "LagrangianDS.hpp"

/** Lagrangian Linear Systems with time invariant and diagonal coefficients - \f$M\dot v + Cv + Kq = F_{ext}(t,z) + p \f$

    \author SICONOS Development Team - copyright INRIA
    \date (Creation) Apr 29, 2004

    The class LagrangianLinearDiagonalDS  allows to define  and compute a generic ndof-dimensional
    Lagrangian Linear Time Invariant Dynamical System of the form :
    where
    \f[
    M \ddot q + C \dot q + K q =  F_{ext}(t,z) + p,
    \f]
    where
    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,
    the variable p contains the impulse and not the force.
    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass().
    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K().
    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C().
    -  \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
    - all matrix coefficients are diagonal

    If required (e.g. for Event-Driven like simulation), reformulation as a first-order system (DynamicalSystem)
    is possible, with:

    - \f$ n= 2 ndof \f$
    - \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
    - rhs given by:
    \f[
    rhs(x,t,z) = \left[\begin{array}{c}
    \dot q  \\
    \ddot q = M^{-1}\left[F_{ext}(t, z) - C \dot q - K q  + p \right]\\
    \end{array}\right]
    \f]
    Its jacobian is:
    \f[
    \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
    0  & I \\
    -M^{-1}K & -M^{-1}C \\
    \end{array}\right]
    \f]
    The input due to the non smooth law is:
    \f[
    r = \left[\begin{array}{c}0 \\ p \end{array}\right]
    \f]
*/
class LagrangianLinearDiagonalDS : public LagrangianDS
{

protected:

  /* serialization hooks */
  ACCEPT_SERIALIZATION(LagrangianLinearDiagonalDS);

  /** stiffness matrix */
  SP::SiconosVector _stiffness;

  /** damping matrix */
  SP::SiconosVector _damping;

  /** mass density*/
  double _mu;

  /** default constructor */
  LagrangianLinearDiagonalDS():LagrangianDS() {};

public:

  /** constructor from initial state and all operators.
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param stiffness diagonal of the stiffness matrix
   *  \param damping diagonal of the damping matrix
   *  \param mass, diagonal of the mass matrix
   */
  LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosVector stiffness, SP::SiconosVector damping, SP::SiconosVector mass);

  /** constructor for complete system with identity mass matrix
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param mass, diagonal of the mass matrix
   *  \param stiffness diagonal of the stiffness matrix
   */
  LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosVector stiffness, SP::SiconosVector damping);

  /** constructor for undamped system and identity mass matrix
   *  \param q0 initial coordinates
   *  \param v0 initial velocity
   *  \param stiffness diagonal of the stiffness matrix
   */
  LagrangianLinearDiagonalDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosVector stiffness);

  /** destructor */
  ~LagrangianLinearDiagonalDS(){};

  /*! @name Attributes access
    @{ */

  /** get a copy of the stiffness matrix (diagonal only)
   *  \return SiconosVector
   */
  inline const SiconosVector get_stiffness() const
  {
    return *_stiffness;
  }

  /** get stiffness matrix (diagonal only, pointer link)
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector stiffness() const
  {
    return _stiffness;
  }

  /** get a copy of the damping matrix (diagonal only)
   *  \return SiconosVector
   */
  inline const SiconosVector get_damping() const
  {
    return *_damping;
  }

  /** get damping matrix (diagonal only, pointer link)
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector damping() const
  {
    return _damping;
  }

  ///@}

  /*! @name Right-hand side computation */
  //@{

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param time of initialization
   */
  void initRhs(double t) ;

  /** Compute \f$F(v,q,t,z)\f$
   *  \param time the current time
   *  \param q SP::SiconosVector: pointers on q
   *  \param velocity SP::SiconosVector: pointers on velocity
   */
  void computeForces(double time, SP::SiconosVector q, SP::SiconosVector velocity);

  ///@}

  /*! @name Miscellaneous public methods */
  //@{

  /**\return true if the Dynamical system is linear. */
  bool isLinear()
  {
    return true;
  }

  /** print the data of the dynamical system on the standard output
   */
  void display() const;

  ///@}
  ACCEPT_STD_VISITORS();

};
#endif // LAGRANGIANLINEARDIAGONALDS_H
