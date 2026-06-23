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

/**
 * \file CohesiveZoneModelNIFNSL.cpp
 * \brief Implementation of the CohesiveZoneModelNIFNSL base class
 *
 * This file implements the common functionality for cohesive zone models,
 * including construction and default implementations of virtual methods.
 *
 * The cohesive zone model extends NewtonImpactFrictionNSL to support
 * interfaces that can sustain traction before contact occurs. This is
 * essential for modeling material failure, delamination, and adhesive
 * contact.
 *
 * Key implementation aspects:
 * - The fallback law (_nslaw_broken) is created during construction with
 *   the same restitution and friction parameters as the cohesive law
 * - isActiveAtLevel() returns true only for level 1 (predictor step),
 *   ensuring cohesive forces are computed before contact detection
 * - Derived classes must implement the cohesive force computation and
 *   internal variable management
 */

#include "CohesiveZoneModelNIFNSL.hpp"

#include <iostream>

#include "Interaction.hpp"
#include "NewtonImpactFrictionNSL.hpp"

namespace siconos::modeling {

CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(siconos::algebra::Index size)
    : NewtonImpactFrictionNSL(size) {
  // Create fallback law for when the interface is fully broken.
  // This uses zero restitution to ensure energy dissipation after failure.
  _nslaw_broken = std::make_shared<NewtonImpactFrictionNSL>(0.0, 0.0, 0.0, size);
}

CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(double en, double et, double mu,
                                                   siconos::algebra::Index size)
    : NewtonImpactFrictionNSL(en, et, mu, size) {
  // Create fallback law with the same parameters as the cohesive law.
  // This ensures consistent behavior after the interface breaks.
  _nslaw_broken = std::make_shared<NewtonImpactFrictionNSL>(en, et, mu, size);
}

bool CohesiveZoneModelNIFNSL::isActiveAtLevel(Interaction& inter, unsigned int level) const {
  // Cohesive forces are computed at level 1 (predictor step) to influence
  // the contact detection and reaction computation at level 0 (corrector).
  // This allows cohesive attraction to bring bodies into contact.
  return (level == 1);
}

void CohesiveZoneModelNIFNSL::display() const {
  // Display base class parameters (restitution, friction)
  NewtonImpactFrictionNSL::display();
  
  // Display cohesive-specific header
  std::cout << "=== CohesiveZoneModelNIFNSL data display ==============================="
            << std::endl;
  std::cout << "(Abstract base class - concrete parameters in derived class)" << std::endl;
  std::cout << "==================================================================" << std::endl;
}

}  // namespace siconos::modeling
