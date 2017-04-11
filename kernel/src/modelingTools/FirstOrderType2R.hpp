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

/*! \file FirstOrderType2R.hpp
  \brief non linear relations: \f$y=h(x,\lambda,z)\quadr=g(\lambda,z)\f$
 */

#ifndef FirstOrderType2R_H
#define FirstOrderType2R_H

#include "FirstOrderR.hpp"

typedef void (*Type2PtrH)(unsigned int, double*, unsigned int, double*, unsigned int, double*);
typedef void (*Type2PtrG)(unsigned int, double*,  unsigned int, double*);

/** FirstOrder Non Linear Relation.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 27
 *
 * Derived from FirstOrderR - See this class for more comments.
 *
 *  Relation for First Order Dynamical Systems, with:
 * \f{eqnarray}
 * y &=& h(x,\lambda,z)\\
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
class FirstOrderType2R : public FirstOrderR
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderType2R);

public:
  
  /** Basic contructor */
  FirstOrderType2R();
  
  /** data constructor
   *  \param pluginh name of the plugin to compute h
   *  \param pluging name of the plugin to compute g
   */
  FirstOrderType2R(const std::string& pluginh, const std::string& pluging);

  /** data constructor
  *  \param pluginh name of the plugin to compute h
  *  \param pluging name of the plugin to compute g
  *  \param pluginJacobianhx name of the plugin to compute the Jacobian of h according to x \f$\nabla_x h\f$
  *  \param pluginJacobianglambda name of the plugin to compute the jacobian of g according to lambda
  */
  FirstOrderType2R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJacobianhx, const std::string& pluginJacobianglambda);

  /** destructor
  */
  ~FirstOrderType2R() {};


  /** initialize the relation (check sizes, memory allocation ...)
   * \param inter the interaction that owns this relation
   * \param DSlink link to DS variable
   * \param workV work vectors to initialize
   * \param workM work matrices to initialize
  */
  virtual void initComponents(Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);

  /** default function to compute y = h(x, lambda, t)
  * \param time the current time
  * \param x the state vector
  * \param lambda XXX
  * \param y the "output" vector
  */
  virtual void computeh(double time, SiconosVector& x, SiconosVector& lambda, SiconosVector& y);

  /** default function to compute g
  * \param time the current time
  * \param lambda XXX
  * \param r the nonsmooth "input" vector
  */
  virtual void computeg(double time, SiconosVector& lambda, SiconosVector& r);

  /** default function to compute \f$\nabla_x h\f$
  *  \param time current time (not used)
  *  \param x the state used to evaluate the jacobian
  *  \param lambda
  *  \param C the matrix used to store the jacobian
  */
  virtual void computeJachx(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& C);
//  virtual void computeJachx(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& C);


//  virtual void computeJachz(double time, Interaction& inter, VectorOfBlockVectors& DSlink, VectorOfVectors& workV, VectorOfSMatrices& workM);
//  virtual void computeJachz(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& D);

  /** default function to compute jacobianG according to lambda
  *  \param time current time (not used)
  *  \param lambda the nonsmooth input used to evaluate the jacobian
  *  \param B the matrix used to store the jacobian
  */
  virtual void computeJacglambda(double time, SiconosVector& lambda, SimpleMatrix& B);
//  virtual void computeJacglambda(double time, SiconosVector& lambda, SiconosVector& z, SimpleMatrix& B);

  /** default function to compute y, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param interProp
  *  \param level not used
  */
  virtual void computeOutput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** default function to compute r, using the data from the Interaction and DS
  *  \param time current time (not used)
  *  \param inter Interaction using this Relation
  *  \param interProp
  *  \param level not used
  */
  virtual void computeInput(double time, Interaction& inter, InteractionProperties& interProp, unsigned int level = 0);

  /** return true if the relation requires the computation of residu
      \return true if residu are required, false otherwise
   */
  virtual bool requireResidu()
  {
    return true;
  }

  virtual void prepareNewtonIteration(Interaction& inter, InteractionProperties& interProp);

  virtual void computeJachlambda(double time, SiconosVector& x, SiconosVector& lambda, SimpleMatrix& D);

  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp);

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderType2R)

#endif
