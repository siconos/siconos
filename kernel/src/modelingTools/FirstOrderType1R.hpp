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

  /** build from plugin for \f$h(x,z)\f$ and \f$g(\lambda, z)\f$
  *  \param pluginh the plugin to compute h
  *  \param pluging the plugin to compute g
  */
  FirstOrderType1R(const std::string& pluginh, const std::string& pluging);

  /** build from plugin for \f$h(x,z)\f$,\f$g(\lambda, z) \f$ and their gradients
  *  \param pluginh the plugin to compute h
  *  \param pluging the plugin to compute g
  *  \param pluginJachx the plugin to compute \f$\nabla_x h\f$
  *  \param pluginJacglambda the plugin to compute \f$\nabla_{\lambda} g\f$
  */
  FirstOrderType1R(const std::string& pluginh, const std::string& pluging, const std::string& pluginJachx, const std::string& pluginJacglambda);

  /** destructor */
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
  *  \param F the matrix used to store the jacobian
  */
void computeJachz(double time, SiconosVector& x, SiconosVector& z, SimpleMatrix& F);

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

  virtual void computeJach(double time, Interaction& inter, InteractionProperties& interProp);

  virtual void computeJacg(double time, Interaction& inter, InteractionProperties& interProp);

  /** return true if the relation requires the computation of residu
      \return true if residu are required, false otherwise
   */
  virtual bool requireResidu()
  {
    return true;
  }

    ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderType1R)

#endif
