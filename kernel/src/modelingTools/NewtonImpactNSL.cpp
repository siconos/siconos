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
#include "NewtonImpactNSL.hpp"

#include <iostream>

NewtonImpactNSL::NewtonImpactNSL(): NonSmoothLaw(1), _e(0.0)
{}

NewtonImpactNSL::NewtonImpactNSL(double e):
  NonSmoothLaw(1), _e(e)
{}

NewtonImpactNSL::NewtonImpactNSL(unsigned int size, double e):
  NonSmoothLaw(size), _e(e)
{}

NewtonImpactNSL::~NewtonImpactNSL()
{}

bool NewtonImpactNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactFrictionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactNSL::display() const
{
  std::cout << "===============================================================================" <<std::endl;
  std::cout << "=== Newton impact (frictionless) non-smooth law coefficient of restitution: " << _e <<std::endl;
  std::cout << "===============================================================================" <<std::endl;
}
