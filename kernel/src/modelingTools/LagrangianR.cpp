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
#include "SiconosVisitor.hpp"

void siconos::modeling::LagrangianR::display() const {
  std::cout << "=====> Relation of type "
            << static_cast<std::underlying_type<RelationType>::type>(_relationType)
            << " and subtype "
            << static_cast<std::underlying_type<RelationSubType>::type>(_subType) << "\n";
  std::cout << " jacobianhOver_q_ :\n" << jacobianhOver_q_view_ << "\n";
}

void siconos::modeling::LagrangianR::accept(
    std::shared_ptr<siconos::internal::SiconosVisitor> tourist) const {
  tourist->visit(*this);
}
