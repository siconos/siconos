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

#include "CohesiveZoneModelNIFNSL.hpp"

#include <iostream>

#include "Interaction.hpp"
#include "NewtonImpactFrictionNSL.hpp"

namespace siconos::modeling {

CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(siconos::algebra::Index size)
    : NewtonImpactFrictionNSL(size) {
  _nslaw_broken = std::make_shared<NewtonImpactFrictionNSL>(0.0, 0.0, 0.0, size);
}

CohesiveZoneModelNIFNSL::CohesiveZoneModelNIFNSL(double en, double et, double mu,
                                                   siconos::algebra::Index size)
    : NewtonImpactFrictionNSL(en, et, mu, size) {
  _nslaw_broken = std::make_shared<NewtonImpactFrictionNSL>(en, et, mu, size);
}

bool CohesiveZoneModelNIFNSL::isActiveAtLevel(Interaction& inter, unsigned int level) const {
  // Default implementation: check if level is 1
  // Derived classes may override this based on internal variables
  return (level == 1);
}

void CohesiveZoneModelNIFNSL::display() const {
  NewtonImpactFrictionNSL::display();
  std::cout << "=== CohesiveZoneModelNIFNSL data display ==============================="
            << std::endl;
  std::cout << "==================================================================" << std::endl;
}

}  // namespace siconos::modeling
