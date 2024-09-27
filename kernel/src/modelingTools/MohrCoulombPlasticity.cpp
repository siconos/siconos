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
#include "MohrCoulombPlasticityNSL.hpp"
#include "SiconosException.hpp"
#include <iostream>

bool siconos::modeling::MohrCoulombPlasticityNSL::isVerified() const
{
  bool res = false;
  // to do
  THROW_EXCEPTION("MohrCoulombPlasticityNSL:: isVerified, not yet implemented!");
  return res;
}

void siconos::modeling::MohrCoulombPlasticityNSL::display() const
{
  std::cout << "=== Mohr Coulomb Plasticity non-smooth law data display ===" << std::endl;
  std::cout << " Cohesion: " << _c << std::endl;
  std::cout << " Friction Angle: " << _phi << std::endl;
  std::cout << "==========================================================" << std::endl;
}
