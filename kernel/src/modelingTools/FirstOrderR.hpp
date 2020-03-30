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

/*! \file FirstOrderR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.hpp"
#include "Interaction.hpp"

/** FirstOrder Relation
 *
 * This is an abstract class for all relation operating on first order systems.
 * It should not be used. Rather, the following classes should be used:
 *
 * - FirstOrderNonlinearR: for fully nonlinear relations: \f$ y = h(t, X, \lambda, Z)\f$, \f$ R = g(t, X, \lambda, Z)\f$.
 * - FirstOrderType2R: specialization with \f$ y = h(t, X, \lambda, Z)\f$, \f$ R = g(t, \lambda, Z)\f$.
 * - FirstOrderType1R: further specialization with \f$ y = h(t, X, Z)\f$, \f$ R = g(t, \lambda, Z)\f$.
 * - FirstOrderLinearR: linear case: \f$ y = C(t)x + D(t)\lambda + F(t) z + e\f$, \f$ R = B(t)\lambda\f$.
 * - FirstOrderLinearR: time-invariant linear case: \f$ y = Cx + D\lambda + F z + e\f$, \f$ R = B\lambda\f$.
 *
 * If the relation involves only one DynamicalSystem, then \f$R = r\f$, \f$X = x\f$, and \f$Z = z\f$.
 * With two, then \f$R = [r_1, r_2] \f$, \f$X = [x_1 x_2] \f$, and \f$Z = [z_1 z_2]\f$.
 *
 * Remember that $y$ and $\lambda$ are relation from the Interaction, and have the same size.
 *
 */
class FirstOrderR : public Relation
{
public:
  enum FirstOrderRDS  {x,z,r, DSlinkSize};
  enum FirstOrderRVec {e,  relationVectorsSize};
  enum FirstOrderRMat {mat_C, mat_D, mat_F, mat_B, mat_K, relationMatricesSize};

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderR);

  /** basic constructor
  *  \param newType the type of the relation
  */
  FirstOrderR(RELATION::SUBTYPES newType): Relation(RELATION::FirstOrder, newType) {}

  /** The following matrices are used if the relation is linear w.r.t to some variables.
   * If the matricesa are, the computation of the Jacobian is not done.
   */

  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t X */
  SP::SimpleMatrix _C;
  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t lambda */
  SP::SimpleMatrix _D;
  /** A matrix to store the constant Jacobian of h(t, X, lambda, Z) w.r.t Z */
  SP::SimpleMatrix _F;

  /** A matrix to store the constant Jacobian of g(t, X, lambda, Z) w.r.t lambda */
  SP::SimpleMatrix _B;
  /** A matrix to store the constant Jacobian of g(t, X, lambda, Z) w.r.t X */
  SP::SimpleMatrix _K;

  /** Continuous memory vector of size of x to call plugin */
  SP::SiconosVector _vec_x;
  /** Continuous memory vector of size of z to call plugin */
  SP::SiconosVector _vec_z;
  /** Continuous memory vector of size of r to call plugin */
  SP::SiconosVector _vec_r;



public:

  /** destructor
  */
  virtual ~FirstOrderR() {};


  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   */
  virtual void initialize(Interaction& inter);

  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction& inter) = 0;

  /** set C to pointer newC
  *  \param newC the C matrix
  */
  inline void setCPtr(SP::SimpleMatrix newC)
  {
    _C = newC;
  }

  /** set B to pointer newB
  *  \param newB the B matrix
  */
  inline void setBPtr(SP::SimpleMatrix newB)
  {
    _B = newB;
  }

  /** set D to pointer newPtr
  *  \param newD the D matrix
  */
  inline void setDPtr(SP::SimpleMatrix newD)
  {
    _D = newD;
  }

  /** set F to pointer newPtr
  *  \param newF the F matrix
  */
  inline void setFPtr(SP::SimpleMatrix newF)
  {
    _F = newF;
  }

  /** get C
  *  \return C matrix
  */
  inline SP::SimpleMatrix C() const
  {
    return _C;
  }

  /** get D
  *  \return D matrix
  */
  inline SP::SimpleMatrix D() const
  {
    return _D;
  }

  /** get F
  *  \return F matrix
  */
  inline SP::SimpleMatrix F() const
  {
    return _F;
  }

  /** get B
  *  \return B matrix
  */
  inline SP::SimpleMatrix B() const
  {
    return _B;
  }

  /** get K
  *  \return K matrix
  */
  inline SP::SimpleMatrix K() const
  {
    return _K;
  }

};
#endif
