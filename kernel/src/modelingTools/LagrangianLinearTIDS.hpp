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

/*! \file LagrangianLinearTIDS.hpp

*/
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.hpp"


/** Lagrangian Linear Systems with time invariant coefficients - Derived from LagrangianDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianLinearTIDS  allows to define  and compute a generic ndof-dimensional
 * Lagrangian Linear Time Invariant Dynamical System of the form :
 * where
 * \f[
 * M \ddot q + C \dot q + K q =  F_{Ext}(t,z) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,
 *       the variable p contains the impulse and not the force.
 *    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass.
 *    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K.
 *    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C.
 *    -  \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *
 *
 *
 * Links with first order DynamicalSystem top-class are:
 *
 * \f$ n= 2 ndof \f$
 * \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
 *
 * The rhs is given by:
 * \f[
 * rhs(x,t,z) = \left[\begin{array}{c}
 *  \dot q  \\
 * \ddot q = M^{-1}(q)\left[F_{Ext}( q , t, z) - C \dot q - K q   \right]\\
 * \end{array}\right]
 * \f]
 * Its jacobian is:
 * \f[
 * \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
 *  0  & I \\
 * -M^{-1}K & -M^{-1}C) \\
 * \end{array}\right]
 * \f]
 *  The input due to the non smooth law is:
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 */
class LagrangianLinearTIDS : public LagrangianDS
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);

  /** specific matrix for a LagrangianLinearTIDS */
  SP::SiconosMatrix _K;

  /** specific matrix for a LagrangianLinearTIDS */
  SP::SiconosMatrix _C;

  /** default constructor */
  LagrangianLinearTIDS():LagrangianDS() {};

public:

  /** constructor from a set of data
   *  \param q0 initial coordinates of this DynamicalSystem
   *  \param v0 initial velocity of this DynamicalSystem
   *  \param M mass matrix of the DynamicalSystem
   *  \param K stiffness matrix of the DynamicalSystem
   *  \param C damping matrix of the DynamicalSystem
   */
  LagrangianLinearTIDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix M, SP::SiconosMatrix K, SP::SiconosMatrix C);

  /** constructor from a set of data
   *  \param q0 initial coordinates of this DynamicalSystem
   *  \param v0 initial velocity of this DynamicalSystem
   *  \param M mass matrix of this DynamicalSystem
   */
  LagrangianLinearTIDS(SP::SiconosVector q0, SP::SiconosVector v0, SP::SiconosMatrix M):
    LagrangianDS(q0, v0, M){};

  /** destructor */
  ~LagrangianLinearTIDS(){};

  /**
   * \return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }

  /** Initialization function for the rhs and its jacobian.
   *  \param t time of initialization
   */
  void initRhs(double t) ;

  // --- GETTERS AND SETTERS ---

  // -- K --
  /** get the value of K
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getK() const
  {
    return *_K;
  }

  /** get K
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix K() const
  {
    return _K;
  }

  /** set the value of K to newValue
   *  \param K new stiffness matrix
   */
  void setK(const SiconosMatrix& K);

  /** set K to pointer newPtr
   *  \param newPtr pointer to the new Stiffness matrix
   */
  void setKPtr(SP::SiconosMatrix newPtr);

  // -- C --
  /** get the value of C
   *  \return SimpleMatrix

  inline const SimpleMatrix getC() const { return *C; }
  */
  /** get C
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix C() const
  {
    return _C;
  }

  /** set C to pointer newPtr
   * \param newPtr pointer to the new damping matrix
   */
  void setCPtr(SP::SiconosMatrix newPtr) ;

  /** Default function to the right-hand side term
   *  \param t current time
   *  \param isDup flag to avoid recomputation of operators
   *
   */
  void computeRhs(double t, bool isDup = false);

  /** function to compute forces with some specific values for q and velocity (ie not those of the current state).
   *  \param t the current time
   *  \param q pointer to positions vector
   *  \param v pointer to velocities vector
   */
  void computeForces(double t, SP::SiconosVector q, SP::SiconosVector v);

  /** Default function to compute forces
   *  \param t the current time
   */
  void computeForces(double t);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param t current time
   *  \param isDup flag to avoid recomputation of operators
   *
   */
  void computeJacobianRhsx(double t, bool isDup = false);

  // --- Miscellaneous ---

  /** print the data onto the screen
   */
  void display() const;

  /** overload LagrangianDS corresponding function
   * \return a double, always zero.
   */
  double dsConvergenceIndicator()
  {
    return 0.0;
  }

  ACCEPT_STD_VISITORS();

};
#endif // LAGRANGIANTIDS_H
