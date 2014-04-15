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

/*! \file FirstOrderType1R.hpp
\brief non linear relations, with y depending on dynamical systems state and r on lambda.
 */

#ifndef FirstOrderType1R_H
#define FirstOrderType1R_H

#include "FirstOrderR.hpp"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*Type1Ptr)(unsigned int, double*, unsigned int, double*, unsigned int, double*);


/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27, 2004
 *
 * Derived from FirstOrderR - See this class for more comments.
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(x,z)\\
 * r &=& g(\lambda,z)
 * \f}
 *
 * Operators (and their corresponding plug-in):
- h: saved in Interaction as y (plug-in: output[0])
- \f$ \nabla_x h \f$: jacobianH[0] ( output[1] )
- g: saved in DS as r ( input[0])
- \f$ \nabla_\lambda g \f$: jacobianG[0] ( input[1] )
 *
 */
class FirstOrderType1R : public FirstOrderR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderType1R);
  

public:

  /** default constructor */
  FirstOrderType1R() : FirstOrderR(RELATION::Type1R) {};

  /** data constructor
  *  \param pluginh the plugin to compute h
  *  \param pluging the plugin to compute g
  */
  FirstOrderType1R(const std::string& pluginh, const std::string& pluging);

  /** data constructor
  *  \param pluginh the plugin to compute h
  *  \param pluging the plugin to compute g
  *  \param pluginJachx the plugin to compute \f$\nabla_x h\f$
  *  \param pluginJacglambda the plugin to compute \f$\nabla_{\lambda} g\f$
  */
  FirstOrderType1R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJachx, const std::string& pluginJacglambda);

  /** destructor
  */
  ~FirstOrderType1R() {};

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
   * \param DSlink link to DS variable
   * \param workV work vectors to initialize
   * \param workM work matrices to initialize
  */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /** default function to compute y = h(x, z, t)
  * \param time the current time
  * \param x the state vector
  * \param z the external input
  * \param y the "output" vector
  */
void computeh(double time, SiconosVector& x, SiconosVector& z, SiconosVector& y);

  /** default function to compute g
  * \param time the current time
  * \param lambda the lambda vector
  * \param z the external input
  * \param r the nonsmooth "input" vector
  */
void computeg(double time, SiconosVector& lambda, SiconosVector& z, SiconosVector& r);

  /** default function to compute \f$\nabla_x h\f$
  *  \param time current time (not used)
  *  \param x the state used to evaluate the jacobian
  *  \param z the extra input used to evaluate the jacobian
  *  \param C the matrix used to store the jacobian
  */
void computeJachx(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& C);

  /** default function to compute \f$\nabla_z h\f$
  *  \param time current time (not used)
  *  \param x the state used to evaluate the jacobian
  *  \param z the extra input used to evaluate the jacobian
  *  \param D the matrix used to store the jacobian
  */
void computeJachz(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& D);

  /** default function to compute jacobianG according to lambda
  *  \param time current time (not used)
  *  \param lambda the nonsmooth input used to evaluate the jacobian
  *  \param z the extra input used to evaluate the jacobian
  *  \param B the matrix used to store the jacobian
  */
void computeJacglambda(double time, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B);

  /** default function to compute y, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param DSlink
  *  \param workV
  *  \param workM
  *  \param level not used
  */
virtual void computeOutput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level = 0);

  /** default function to compute r, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param DSlink
  *  \param workV
  *  \param workM
  *  \param level not used
  */
virtual void computeInput(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM, SiconosMatrix& osnsM, unsigned int level = 0);

    virtual void computeJach(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

    virtual void computeJacg(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

    ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderType1R)

#endif
