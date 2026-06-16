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
/*! \file CohesiveZoneModelNIFNSL.hpp
  Base class for cohesive zone models based on NewtonImpactFrictionNSL

  A cohesive zone model describes the progressive damage and failure of an interface
  between two materials. When the interface is intact, it can sustain traction up to
  a critical value. Beyond this threshold, the interface progressively degrades until
  complete failure.
*/

#ifndef COHESIVEZONEMODELNIFNSL_H
#define COHESIVEZONEMODELNIFNSL_H

#include "NewtonImpactFrictionNSL.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {

// Forward declaration
class Interaction;

/** CohesiveZoneModelNIFNSL
 * Base class for cohesive zone models based on NewtonImpactFrictionNSL
 * 
 * This class extends the Newton impact-friction law with cohesive behavior,
 * where the interface can sustain traction up to a critical value before
 * progressive damage occurs.
 */
class CohesiveZoneModelNIFNSL : public NewtonImpactFrictionNSL {
 private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(CohesiveZoneModelNIFNSL);

  /** fallback law when the interface is broken */
  std::shared_ptr<NonSmoothLaw> _nslaw_broken;

 protected:
  /** default constructor
   */
  CohesiveZoneModelNIFNSL() = delete;  // Base class has private default constructor

 public:
  /** basic constructor
   *  \param size size of the ns law
   */
  explicit CohesiveZoneModelNIFNSL(siconos::algebra::Index size);

  /** constructor with the value of the CohesiveZoneModelNIFNSL attributes
   *  \param en double : normal restitution coefficient
   *  \param et double : tangent restitution coefficient
   *  \param mu double : friction coefficient
   *  \param size unsigned int: size of the ns law
   */
  CohesiveZoneModelNIFNSL(double en, double et, double mu, siconos::algebra::Index size);

  /** Destructor */
  ~CohesiveZoneModelNIFNSL() noexcept override = default;

  /** getter of r_cohesion - pure virtual
   * \param inter the Interaction containing internal variables
   * \return pointer to the cohesion force vector data
   */
  virtual double* r_cohesion(Interaction& inter) const = 0;

  /** getter for the broken law
   * \return the fallback non-smooth law when interface is broken
   */
  std::shared_ptr<NonSmoothLaw> nslawBroken() const { return _nslaw_broken; };

  /** Check if the NS law is active at a given level
   * \param inter the Interaction
   * \param level the level to check
   * \return true if active
   */
  /** Check if the NS law is active at a given level
   *  \param inter the Interaction
   *  \param level the level to check
   *  \return true if active
   */
  bool isActiveAtLevel(Interaction& inter, unsigned int level) const override;

  /** print the data to the screen
   */
  void display() const override;

  /** Visitors hook
   */
  void accept(nonsmooth_laws::Visitor& tourist) const override {
    tourist.visit(*this);
  }

  Type acceptType(types::FindType& ft) const override { return ft.visit(*this); }
};

}  // namespace siconos::modeling

#endif  // COHESIVEZONEMODELNIFNSL_H
