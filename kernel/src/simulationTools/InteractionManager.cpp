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

#include "InteractionManager.hpp"

#include <algorithm>

#include "NonSmoothLaw.hpp"
#include "siconos_debug.h"

void siconos::simulation::InteractionManager::insertNonSmoothLaw(
    std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw, long unsigned int group1,
    long unsigned int group2)
{
  // ublas::matrix size type is not the same on 32 bits and 64 bits
  auto maxgroup = std::max((siconos::modeling::NSLawMatrix::size_type)group1,
                           (siconos::modeling::NSLawMatrix::size_type)group2);
  _nslaws.resize(std::max(_nslaws.size1(), maxgroup + 1));
  _nslaws(group1, group2) = nslaw;
}

std::shared_ptr<siconos::modeling::NonSmoothLaw>
siconos::simulation::InteractionManager::nonSmoothLaw(long unsigned int group1,
                                                      long unsigned int group2)
{
  if (group1 < _nslaws.size1() && group2 < _nslaws.size2())
    return _nslaws(group1, group2);
  else
    return std::shared_ptr<siconos::modeling::NonSmoothLaw>();  // ??
}
