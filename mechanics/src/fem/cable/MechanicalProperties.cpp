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

#include "MechanicalProperties.h"

#include <iostream>

void siconos::fem::cable::MechanicalProperties::set_T(double a_T) {
  if (a_T > 0.0) initialTension_ = a_T;
}

void siconos::fem::cable::MechanicalProperties::set_rho(double a_rho) {
  if (a_rho > 0.0) linearDensity_ = a_rho;
}

void siconos::fem::cable::MechanicalProperties::display() const {
  std::cout << "--- Mechanical properties: \n";
  std::cout << "- Rigidity: " << crossSectionRigidity_;
  std::cout << " - Linear density: " << linearDensity_;
  std::cout << " - Initial tension: " << initialTension_;
  std::cout << "\n --------------\n\n";
}