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
/*! \file MohrCoulombPlasticityNSL.hpp
  Newton-Impact Non-Smooth Law
*/

#ifndef MOHRCOULOMBPLASTICITYNSLAW_H
#define MOHRCOULOMBPLASTICITYNSLAW_H

#include "NonSmoothLaw.hpp"

namespace siconos::modeling {

/** Newton Impact-Friction Non Smooth Law
 *
 */
class MohrCoulombPlasticityNSL : public NonSmoothLaw {
 private:
  ACCEPT_SERIALIZATION(NewtonImpactFrictionNSL);

  /** The Newton coefficient of restitution
   */
  double _c{0.}; // Cohesion
  double _phi{0.}; // Friction Angle

  /** default constructor
   */
  MohrCoulombPlasticityNSL() = default;

 public:
  /** basic constructor
   *
   *  \param size size of the ns law
   */
  MohrCoulombPlasticityNSL(unsigned int size) : NonSmoothLaw(size){};

  /** constructor with the value of the NewtonImpactFrictionNSL attributes
   *
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param size unsigned int: size of the ns law
   */
  MohrCoulombPlasticityNSL(double c, double phi, unsigned int size)
      : NonSmoothLaw(size), _c(c), _phi(phi){};

  /** Destructor */
  ~MohrCoulombPlasticityNSL() noexcept = default;

  /** check the ns law to see if it is verified
   *
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const override;

  // GETTERS/SETTERS

  /** \return the value of en
   */
  inline double c() const { return _c; };

  /** setter of en
   *
   *  \param newVal a double to set en
   */
  inline void setCohesion(double newVal) { _c = newVal; };

  /** \return the value of et
   */
  inline double phi() const { return _phi; };

  /** setter of et
   *
   *  \param newVal a double to set et
   */
  inline void setFrictionAngle(double newVal) { _phi = newVal; };

  // OTHER FUNCTIONS

  /** print the data to the screen
   */
  void display() const override;

  // visitors hook
  void accept(siconos::internal::SiconosVisitor& tourist) const override {
    tourist.visit(*this);
  }
  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};
}  // namespace siconos::modeling
#endif  // NewtonImpactFrictionNSL_H
