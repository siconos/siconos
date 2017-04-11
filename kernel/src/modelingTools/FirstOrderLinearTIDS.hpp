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
/*! \file FirstOrderLinearTIDS.hpp

*/
#ifndef LINEARTIDS_H
#define LINEARTIDS_H

#include "FirstOrderLinearDS.hpp"


/** First order linear and time-invariant coeff systems - \f$M \dot x = Ax(t)+ b + r, x(t_0)=x_0\f$.
 
   \author SICONOS Development Team - copyright INRIA
   \date (Creation) Apr 29, 2004
 
 
   This class represents first order linear systems of the form:
 
   \f[
   M\dot x(t) = A x(t) + b + r,
   x(t_0)=x_0
   \f]

   where
   - \f$x \in R^{n} \f$ is the state,
   - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
   - \f$M \in R^{n\times n} \f$ is a constant invertible matrix
   - \f$A \in R^{n\times n}\f$
   - \f$b \in R^{n} \f$

   No plugged operators for this class.
   
**/

class FirstOrderLinearTIDS : public FirstOrderLinearDS
{
private:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIDS);

  /** default constructor
   */
  FirstOrderLinearTIDS() {};

public:

  /** initial state and constant A matrix
   *  \param x0 the initial state vector
   *  \param A the A matrix
   */
  FirstOrderLinearTIDS(SP::SiconosVector x0, SP::SiconosMatrix A):
    FirstOrderLinearDS(x0, A){};

  /** initial state, constant A matrix, constant b vector
   *  \param x0 the initial state vector
   *  \param A matrix
   *  \param b vector
   */
  FirstOrderLinearTIDS(SP::SiconosVector x0, SP::SiconosMatrix A, SP::SiconosVector b):
    FirstOrderLinearDS(x0, A, b){};

  /** Copy constructor
   * \param FOLTIDS the FirstOrderLinearTIDS to copy
   */
  FirstOrderLinearTIDS(const FirstOrderLinearTIDS & FOLTIDS): FirstOrderLinearDS(FOLTIDS) {};

  /** destructor */
  ~FirstOrderLinearTIDS() {};

  /*! @name Right-hand side computation */
  
  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization.
   */
  void initRhs(double time);

  /** Default function to the right-hand side term
   *  \param time current time
   *  \param isDSup flag to avoid recomputation of operators
   *
   */
  void computeRhs(double time, bool isDSup = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param time current time
   *  \param isDSup flag to avoid recomputation of operators
   *
   */
  void computeJacobianRhsx(double time, bool isDSup = false);

  ///@}

  /*! @name Miscellaneous public methods */

  /** data display on screen
   */
  void display() const;

  ///@}

  /** Dumb function, there is no plugin here
   * \param time unused
   */
  virtual void updatePlugins(double time)
  {};

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIDS)

#endif // LINEARTIDS_H
