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

#include <cassert>
#include <iostream>
#include "NormalConeNSL.hpp"

#include "SimpleMatrix.hpp"
#include "SiconosVector.hpp"


// Default (private)
NormalConeNSL::NormalConeNSL(): NonSmoothLaw(), _H(std11::shared_ptr<SimpleMatrix>()), _K(std11::shared_ptr<SiconosVector>())
{}

NormalConeNSL::NormalConeNSL(unsigned size, SP::SimpleMatrix H, SP::SiconosVector K):
  NonSmoothLaw(size), _H(H), _K(K)
{
assert(H->size(1) == size &&
      "NormalConeNSL::NormalConeNSL - the number of columns in H and the declared size are not equal, check your code !");
}

NormalConeNSL::~NormalConeNSL()
{}

bool NormalConeNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void NormalConeNSL::display() const
{
  std::cout << "------------------------------------" <<std::endl;
  std::cout << "____ data of the NormalConeNSL" <<std::endl;
  std::cout << "| nSLawSize : " << _size <<std::endl;
  std::cout << "| H : " << std::endl;
  _H->display();
  std::cout << "| K : " << std::endl;
  _K->display();
  std::cout << "____________________________" <<std::endl;
  std::cout << "------------------------------------" <<std::endl;
}
