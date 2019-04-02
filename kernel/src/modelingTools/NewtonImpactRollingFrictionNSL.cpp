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
#include "NewtonImpactRollingFrictionNSL.hpp"

#include <iostream>

// Default (private)
NewtonImpactRollingFrictionNSL::NewtonImpactRollingFrictionNSL():
  NonSmoothLaw(), _en(0.0), _et(0.0), _mu(0.0), _muR(0.0)
{}

NewtonImpactRollingFrictionNSL::NewtonImpactRollingFrictionNSL(unsigned int size):
  NonSmoothLaw(size), _en(0.0), _et(0.0), _mu(0.0), _muR(0.0)
{}

NewtonImpactRollingFrictionNSL::NewtonImpactRollingFrictionNSL(double newEn, double newEt, double newMu, double newMuR, unsigned int newSize):
  NonSmoothLaw(newSize), _en(newEn), _et(newEt), _mu(newMu), _muR(newMuR)
{}

NewtonImpactRollingFrictionNSL::~NewtonImpactRollingFrictionNSL()
{}

bool NewtonImpactRollingFrictionNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactRollingFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactRollingFrictionNSL::display() const
{
  std::cout << "=== Newton impact-friction non-smooth law data display ===" <<std::endl;
  std::cout << " Normal Newton coefficient of restitution: " << _en <<std::endl;
  std::cout << " Tangential Newton coefficient of restitution: " << _et <<std::endl;
  std::cout << "Friction coefficient: " << _mu <<std::endl;
  std::cout << "Rolling friction coefficient: " << _muR <<std::endl;
  std::cout << "==========================================================" <<std::endl;
}
