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

/*! \file FirstOrderNonLinearR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderNonLinearR_H
#define FirstOrderNonLinearR_H

#include "FirstOrderR.hpp"
#include "Interaction.hpp"

/** Pointer to function for plug-in for operators related to output and its gradients.*/
typedef void (*OutPtr)(unsigned int, double*, double, unsigned int, double*, double*, unsigned int, double*);

/** Pointer to function for plug-in for operators related to input and its gradients.*/
typedef void (*InPtr)(unsigned int, double*, double, unsigned int, double*, unsigned int, double*);

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

class FirstOrderNonLinearR : public FirstOrderR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearR);


  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

public:

  /** basic constructor
  */
  FirstOrderNonLinearR(): FirstOrderR(RELATION::NonLinearR) {}

  /** destructor
  */
  virtual ~FirstOrderNonLinearR() {};

  /** default function to compute h
   * \param time : current time
   * \param x
   * \param lambda
   * \param z
   * \param y
   */
  virtual void computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& y);

  /** default function to compute g
   * \param time : current time
   * \param x
   * \param lambda
   * \param z
   * \param r
  */
  virtual void computeg(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SiconosVector& r);

  /** default function to compute jacobianH
   * \param time : current time
   * \param x
   * \param lambda
   * \param z
   * \param C
   */
  virtual void computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& C);
  virtual void computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& D);
  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  /** default function to compute jacobianG according to lambda
  * \param time current time
  * \param x
  * \param lambda
  * \param z
  * \param B
  */
  virtual void computeJacglambda(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B);
  /** default function to compute jacobianG according to x
   * \param time  double : current time
   * \param x
   * \param lambda
   * \param z
   * \param K
   */
  virtual void computeJacgx(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& K);

  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp);
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);
  virtual void prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp);

  /** return true if the relation requires the computation of residu
   * \return true if residu are required, false otherwise
   */
  virtual bool requireResidu()
  {
    return true;
  }
 
};
TYPEDEF_SPTR(FirstOrderNonLinearR)

#endif
