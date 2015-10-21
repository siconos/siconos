/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file FirstOrderR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderR_H
#define FirstOrderR_H

#include "Relation.hpp"
#include "Interaction.hpp"

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(X,t,\lambda,Z)\\
 * R &=& g(X,t,\lambda,Z)
 * \f}
 *  X, Z, R corresponds to DynamicalSystem variables.
 *  If DS1 and DS2 are involved in the linked Interaction, then X =[x1 x2], Z=[z1 z2] ...
 *
 *  \f$ y \ and \ \lambda \f$ are specific variables of the interaction (see this class for more details).
 *  h and g are plugged on external functions, via plug-in mechanism (see SiconosSharedLibrary).
 *
 * h <=> output
 *
 * g <=> input
 *
 * Operators (and their corresponding plug-in):
- h: saved in Interaction as y (plug-in: output[0])
- \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
- \f$ \nabla_\lambda h \f$: jacobianH[1] ( output[2] )
- g: saved in DS as r ( input[0])
- \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )


Note: we use a vector for jacobianG while there is only one jacobian. Just for future changes and to allow easy new implementations if some other
variables are required in g.

 *
 */
// namespace FirstOrderR {enum {xfree, z, x, r, deltax, xPartialNS, DSlinkSize};}
// namespace FirstOrderRVec {enum {xfree, z, x, r, e, g_alpha, residuR, workVecSize};}
// namespace FirstOrderRMat {enum {C, D, F, B, K, Ktilde, Khat, workMatSize};}


class FirstOrderR : public Relation
{
public:
  enum FirstOrderRDS {xfree, z, x, r, deltax, xPartialNS, DSlinkSize};
  enum FirstOrderRVec {vec_xfree, vec_z, vec_x, vec_r, e, g_alpha, vec_residuR, workVecSize};
  enum FirstOrderRMat  {mat_C, mat_D, mat_F, mat_B, mat_K, mat_Ktilde, mat_Khat, mat_workMatSize};


protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderR);

  /** basic constructor
  *  \param newType the type of the relation
  */
  FirstOrderR(RELATION::SUBTYPES newType): Relation(RELATION::FirstOrder, newType) {}

  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM) = 0;

  SP::SimpleMatrix _C;
  SP::SimpleMatrix _D;
  SP::SimpleMatrix _F;

  SP::SimpleMatrix _B;
  SP::SimpleMatrix _K;

public:

  /** destructor
  */
  virtual ~FirstOrderR() {};

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
   * \param DSlink
   * \param workV
   * \param workM
   */
  virtual void initialize(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

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
