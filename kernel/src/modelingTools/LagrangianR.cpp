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

// \todo : create a work vector for all tmp vectors used in computeg, computeh ...

#include "LagrangianR.hpp"

#include <iostream>

#include "SiconosMatrix.hpp"
#include "Tools.hpp"

void siconos::modeling::LagrangianR::display() const {
  std::cout << "=====> Relation of type " << tools::enum_to_string(_relationType)
            << " and subtype " << tools::enum_to_string(_subType) << "\n";
  if (jacobianhOver_q_view_)
    std::cout << " jacobianhOver_q_ :\n" << *jacobianhOver_q_view_ << "\n";
}
