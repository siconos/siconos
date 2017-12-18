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
/*! \file NewtonImpactNSL.hpp

*/
#ifndef NEWTONIMPACTNSL_H
#define NEWTONIMPACTNSL_H

#include "NonSmoothLaw.hpp"

/** Newton impact Non Smooth Law
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) June 29, 2004
 *
 * This class formalizes the Newton Impact law together with a complementarity condition. i.e.
 * \f[
 * \left\{\begin{array}{l}
 * y \geq 0, \lambda \geq 0, y^{T} \lambda=0\\
 *  if y \leq 0 \quad \mbox{then} \quad \dot y(t^{+}) - e \dot y(t^{-}) \geq 0, \quad  \lambda \geq 0, (\dot y(t^{+}) - e \dot y(t^{-}))^{T} \lambda=0
 * \end{array}\right.
 * \f]
 *
 * nsLawSize is equal to 1.
 *
 */
class NewtonImpactNSL : public NonSmoothLaw
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonImpactNSL);

  /** The Newton normal coefficient of restitution  */
  double _e;

public:

  /** default constructor
   */
  NewtonImpactNSL();

  /** constructor with the value of the NewtonImpactNSL attributes
  *  \param e the value of the coefficient of restitution
  */
  NewtonImpactNSL(double e);

  /** Apply multiple-axis impact */
  NewtonImpactNSL(unsigned int size, double e);

  /** destructor
   */
  ~NewtonImpactNSL();

  /** check the ns law to see if it is verified
  *  \return a boolean value whioch determines if the NS Law is verified
  */
  bool isVerified() const;

  /** getter of e
  *  \return the value of e
  */
  inline double e() const
  {
    return _e;
  };

  /** setter of e
  *  \param newVal a double to set e
  */
  inline void setE(double newVal)
  {
    _e = newVal;
  };

  /** print the data to the screen
  */
  void display() const;

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();
};

DEFINE_SPTR(NewtonImpactNSL)


#endif // NewtonImpactNSL_H
