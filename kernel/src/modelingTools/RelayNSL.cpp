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

#include "RelayNSL.hpp"

#include <iostream>

// Default (private)
RelayNSL::RelayNSL(): NonSmoothLaw(), _lb(-1.0), _ub(1.0)
{}

RelayNSL::RelayNSL(unsigned int size, double lb, double ub):
  NonSmoothLaw(size), _lb(lb), _ub(ub)
{}

RelayNSL::~RelayNSL()
{}

bool RelayNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void RelayNSL::display() const
{
  std::cout << "------------------------------------" <<std::endl;
  std::cout << "____ data of the RelayNSL" <<std::endl;
  std::cout << "| nSLawSize : " << _size <<std::endl;
  std::cout << "| lb : " << _lb <<std::endl;
  std::cout << "| ub : " << _ub <<std::endl;
  std::cout << "____________________________" <<std::endl;
  std::cout << "------------------------------------" <<std::endl;
}
