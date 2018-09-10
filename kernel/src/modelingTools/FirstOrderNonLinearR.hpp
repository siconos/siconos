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

/*! \file FirstOrderNonLinearR.hpp
\brief General interface for relations.
 */

#ifndef FirstOrderNonLinearR_H
#define FirstOrderNonLinearR_H

#include "FirstOrderR.hpp"
#include "Interaction.hpp"

/** Pointer to function for plug-in for operators related to input and its gradients.*/

/** FirstOrder Non Linear Relation.
 *
 *  Relation for First Order Dynamical Systems, with:
 * \rststar
 *
 * .. math::
 *
 *     y &=& h(X,t,\lambda,Z) \\
 *     R &=& g(X,t,\lambda,Z)
 *
 * \endrststar
 *
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




public:

  /** basic constructor
  */
  FirstOrderNonLinearR(): FirstOrderR(RELATION::NonLinearR) {}

  /** destructor
  */
  virtual ~FirstOrderNonLinearR() {};

  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction using this relation
   */
  void initialize(Interaction& inter);

  /** check sizes of the relation specific operators.
   * \param inter an Interaction using this relation
   */
  virtual void checkSize(Interaction& inter);

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

  virtual void computeJach(double time, Interaction& inter);

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

  virtual void computeJacg(double time, Interaction& inter);

  /** default function to compute y, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param level not used
  */
  virtual void computeOutput(double time, Interaction& inter, unsigned int level = 0);

  /** default function to compute r, using the data from the Interaction and DS
   *  \param time current time (not used)
   *  \param inter Interaction using this Relation
   *  \param level not used
   */
  virtual void computeInput(double time, Interaction& inter, unsigned int level = 0);

  
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
