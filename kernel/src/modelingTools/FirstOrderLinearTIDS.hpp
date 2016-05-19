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


/** First order linear systems - Inherits from DynamicalSystems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * \dot x(t) = A x(t) + b + r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *
 *  The  rhs is specialized by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$b \in R^{n} \f$
 *
 *
 * This class inherits from FirstOrderLinearDS one. The difference is that here A and b do not depend on time.
 *
 *
 * A and b are constant matrix or vector, and thus can not be set using a plug-in.
 *
 * To build and use a linearDS, you first need to call a constructor, with A as a required data.
 * Then, this system has to be initialized -> compute members value at time t0. This is usually done during call to simulation->initialize.
 * Finally, the state of the DS can be obtained by calling "compute" functions. In FirstOrderLinearTIDS case, since A and b are fixed, you can
 * only call computeRhs(time), to compute rhs = Ax+b.
 *
 **/
class FirstOrderLinearTIDS : public FirstOrderLinearDS
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderLinearTIDS);


  /** default constructor
   */
  FirstOrderLinearTIDS() {};

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** constructor from a set of data
   *  \param x0 the initial state of this DynamicalSystem
   *  \param A the A matrix
   */
  FirstOrderLinearTIDS(SP::SiconosVector x0, SP::SiconosMatrix A);

  /** constructor from a set of data
   *  \param x0 the initial state of this DynamicalSystem
   *  \param A the A matrix
   *  \param b the b vector
   */
  FirstOrderLinearTIDS(SP::SiconosVector x0, SP::SiconosMatrix A, SP::SiconosVector b);

  /** Copy constructor
   * \param FOLTIDS the FirstOrderLinearTIDS to copy
   */
  FirstOrderLinearTIDS(const FirstOrderLinearTIDS & FOLTIDS): FirstOrderLinearDS(FOLTIDS) {};

  /** destructor */
  ~FirstOrderLinearTIDS() {};

  /** indicate that the DS is linear
   * \return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }



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

  /** data display on screen
   */
  void display() const;

  /** Dumb function, there is no plugin here
   * \param time unused
   */
  virtual void updatePlugins(double time)
  {
    ;
  };

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearTIDS)

#endif // LINEARTIDS_H
