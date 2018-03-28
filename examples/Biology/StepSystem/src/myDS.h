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

/*! \file MyDSDS.h
 The considered dynamical system is a first order non linear system,
 of size 2. See details in siconos examples manual.

*/

#ifndef MYDSDS_H
#define MYDSDS_H

#include "SiconosKernel.hpp"

/** This class inherits from first Order linear dynamical system
 */
class MyDS : public FirstOrderNonLinearDS
{

public:
  /** default constructor
   * \param the type of the system
   */
  MyDS(SP::SiconosVector x0);


  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~MyDS() {};

  using FirstOrderNonLinearDS::computef;


  /** Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   */
  virtual void computef(double, SP::SiconosVector x);

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  //virtual void computef(double, SP::SiconosVector);

  using FirstOrderNonLinearDS::computeJacobianfx;


  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  //virtual void computeJacobianfx(double, bool  = false);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$ with x different from current saved state.
   *  \param double time : current time
   *  \param SP::SiconosVector
   */
  virtual void computeJacobianfx (double, SP::SiconosVector);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  /* virtual void computeRhs(double, bool  = false); */

  using FirstOrderNonLinearDS::resetNonSmoothPart;

  virtual void resetNonSmoothPart();

};

TYPEDEF_SPTR(MyDS);

#endif


