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
/*! \file NewtonImpactRollingFrictionNSL.hpp
  Newton-Impact Non-Smooth Law
*/

#ifndef NEWTONIMPACTROLLINGFRICTIONNSLAW_H
#define NEWTONIMPACTROLLINGFRICTIONNSLAW_H

#include "NonSmoothLaw.hpp"

/** Newton Impact-Friction Non Smooth Law
 *
 */
class NewtonImpactRollingFrictionNSL : public NonSmoothLaw
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonImpactRollingFrictionNSL);

  /** The Newton coefficient of restitution
   */
  double _en;
  double _et;
  /** friction coefficient */
  double _mu;

  /** friction coefficient */
  double _muR;

  /** default constructor
   */
  NewtonImpactRollingFrictionNSL();

public:

  /** basic constructor
   *  \param size size of the ns law
   */
  explicit NewtonImpactRollingFrictionNSL(unsigned int size);

  /** constructor with the value of the NewtonImpactRollingFrictionNSL attributes
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param muR double : rolling friction coefficient
   *  \param size unsigned int: size of the ns law
   */
  NewtonImpactRollingFrictionNSL(double en, double et, double mu, double muR, unsigned int size);

  /** Destructor */
  ~NewtonImpactRollingFrictionNSL();

  /** check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  // GETTERS/SETTERS

  /** getter of en
   *  \return the value of en
   */
  inline double en() const
  {
    return _en;
  };

  /** setter of en
   *  \param newVal a double to set en
   */
  inline void setEn(double newVal)
  {
    _en = newVal;
  };

  /** getter of et
   *  \return the value of et
   */
  inline double et() const
  {
    return _et;
  };

  /** setter of et
   * \param newVal a double to set et
   */
  inline void setEt(double newVal)
  {
    _et = newVal;
  };

  /** getter of mu
   * \return the value of mu
   */
  inline double mu() const
  {
    return _mu;
  };

  /** getter of mu
   * \return the value of mu
   */
  inline double muR() const
  {
    return _muR;
  };

  /** setter of mu
   * \param newVal a double to set mu
   */
  inline void setMu(double newVal)
  {
    _mu = newVal;
  };

  /** setter of muR
   * \param newVal a double to set muR
   */
  inline void setMuR(double newVal)
  {
    _muR = newVal;
  };

  // OTHER FUNCTIONS

  /** print the data to the screen
   */
  void display() const;

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();

};
DEFINE_SPTR(NewtonImpactRollingFrictionNSL)
#endif // NewtonImpactRollingFrictionNSL_H
