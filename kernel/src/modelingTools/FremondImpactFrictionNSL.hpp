/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
/*! \file NewtonImpactFrictionNSL.hpp
  Fremond-Impact Non-Smooth Law
*/

#ifndef FREMONDIMPACTFRICTIONNSLAW_H
#define FREMONDIMPACTFRICTIONNSLAW_H

#include "NonSmoothLaw.hpp"

/** Fremond Impact-Friction Non Smooth Law
 *
 */
class FremondImpactFrictionNSL : public NonSmoothLaw {

private:
  ACCEPT_SERIALIZATION(FremondImpactFrictionNSL);

  /** The Fremond coefficient of restitution
   */
  double _en;
  double _et;
  /** friction coefficient */
  double _mu;

  /** default constructor
   */
  FremondImpactFrictionNSL();

public:
  /** basic constructor
   *
   *  \param size size of the ns law
   */
  FremondImpactFrictionNSL(unsigned int size);

  /** constructor with the value of the FremondImpactFrictionNSL attributes
   *
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param size unsigned int: size of the ns law
   */
  FremondImpactFrictionNSL(double en, double et, double mu, unsigned int size);

  /** Destructor */
  ~FremondImpactFrictionNSL();

  /** check the ns law to see if it is verified
   *
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const override;

  // GETTERS/SETTERS

  /** \return the value of en
   */
  inline double en() const { return _en; };

  /** setter of en
   *
   *  \param newVal a double to set en
   */
  inline void setEn(double newVal) { _en = newVal; };

  /** \return the value of et
   */
  inline double et() const { return _et; };

  /** setter of et
   *
   *  \param newVal a double to set et
   */
  inline void setEt(double newVal) { _et = newVal; };

  /** \return the value of mu
   */
  inline double mu() const { return _mu; };

  /** setter of mu
   *
   *  \param newVal a double to set mu
   */
  inline void setMu(double newVal) { _mu = newVal; };

  // OTHER FUNCTIONS

  /** print the data to the screen
   */
  void display() const override;

  ACCEPT_STD_VISITORS();
};
DEFINE_SPTR(FremondImpactFrictionNSL)
#endif // FremondImpactFrictionNSL_H
