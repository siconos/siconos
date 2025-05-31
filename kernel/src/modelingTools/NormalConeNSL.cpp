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

#include "NormalConeNSL.hpp"

#include <cassert>
#include <iostream>

#include "SiconosException.hpp"
#include "SiconosVector.hpp"
#include "SiconosMatrix.hpp"

siconos::modeling::NormalConeNSL::NormalConeNSL(
    unsigned size, std::shared_ptr<siconos::algebra::SiconosMatrix> H,
    std::shared_ptr<siconos::algebra::SiconosVector> K)
    : NonSmoothLaw(size), _H(H), _K(K)
{
  assert(H->cols() == size &&
         "NormalConeNSL::NormalConeNSL - the number of columns in H and the declared size are "
         "not equal, check your code !");
}

bool siconos::modeling::NormalConeNSL::isVerified(void) const
{
  bool res = false;
  // to do
  return res;
}

void siconos::modeling::NormalConeNSL::display() const
{
  std::cout << "------------------------------------" << std::endl;
  std::cout << "____ data of the NormalConeNSL" << std::endl;
  std::cout << "| nSLawSize : " << _size << std::endl;
  std::cout << "| H : " << std::endl;
  siconos::algebra::print(*_H);
  std::cout << "| K : " << std::endl;
  siconos::algebra::print(*_K);
  std::cout << "____________________________" << std::endl;
  std::cout << "------------------------------------" << std::endl;
}
