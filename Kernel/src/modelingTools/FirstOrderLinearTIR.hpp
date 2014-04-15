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
/*! \file FirstOrderLinearTIR.hpp

 */
#ifndef FirstOrderLinearTIR_H
#define FirstOrderLinearTIR_H

#include "FirstOrderR.hpp"

/** Linear Time Invariant Relation, derived from class FirstOrderR

\author SICONOS Development Team - copyright INRIA
\version 3.0.0.
\date Apr 15, 2007

Linear Relation for First Order Dynamical Systems:

\f{eqnarray}
y &=& Cx(t) + Fz + D\lambda + e\\

R &=& B\lambda
\f}

 */
class FirstOrderLinearTIR : public FirstOrderR
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIR);

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
  */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  SP::SiconosVector _e;

public:

  /** default constructor, protected
  */
  FirstOrderLinearTIR();

  /** create the Relation from a set of data
  *  \param C the matrix C
  *  \param B the matrix B
  */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix B);

  /** create the Relation from a set of data
  *  \param C the C matrix
  *  \param D the D matrix
  *  \param F the F matrix
  *  \param e the e matrix
  *  \param B the B matrix
  */
  FirstOrderLinearTIR(SP::SimpleMatrix C, SP::SimpleMatrix D, SP::SimpleMatrix F, SP::SiconosVector e, SP::SimpleMatrix B);

  /** destructor
  */
  virtual ~FirstOrderLinearTIR() {};

  // GETTERS/SETTERS

  /** default function to compute h
  *  \param time current time
  *  \param xXXX
  *  \param zXXX
  *  \param y value of h
  */
  void computeh(BlockVector& x, SiconosVector& lambda, BlockVector& z, SiconosVector& y);

  /** default function to compute g
  *  \param lambdaXXX
  *  \param r non-smooth input
  */
  void computeg(SiconosVector& lambda, BlockVector& r);

  /** default function to compute y
  *  \param time current time
  *  \param inter Interaction using this Relation
  *  \param DSlink
  *  \param workV
  *  \param workM
  *  \param level not used
  */
  virtual void computeOutput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level = 0);

  /** default function to compute r
  *  \param time current time
  *  \param inter Interaction using this Relation
  *  \param DSlink
  *  \param workM
  *  \param level not used
  */
  virtual void computeInput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level = 0);

  /** print the data to the screen
  */
  void display() const;

  /** compute the jacobian of h: nothing to be done here
   */
  virtual void computeJach(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM) {};

  /** compute the jacobian of g: nothing to be done here
   */
  virtual void computeJacg(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM) {};

  /** set e
  *  \param  newe the new value of e
  */
  inline void setePtr(SP::SiconosVector newe)
  {
    _e = newe;
  }

  /** get e
  *  \return e matrix
  */
  inline SP::SiconosVector e() const
  {
    return _e;
  }

  /**
  * return true if the relation is linear.
  */
  virtual bool isLinear()
  {
    return true;
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIR)

#endif
