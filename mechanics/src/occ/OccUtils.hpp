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
#ifndef OCC_UTILS
#define OCC_UTILS

#include <Standard_TypeDef.hxx>  // From Opencascade
#include <memory>

#include "ContactShapeDistance.hpp"

// OpenCascade forward declarations
class TopoDS_Shape;

namespace siconos::mechanics::occ {

struct OccContactFace;
struct OccContactEdge;

void occ_move(TopoDS_Shape& shape, const std::array<double, 7>& pos);

auto occ_distanceFaceFace(std::shared_ptr<OccContactFace> csh1,
                          std::shared_ptr<OccContactFace> csh2) -> ContactShapeDistance;

auto occ_distanceFaceEdge(std::shared_ptr<OccContactFace> csh1,
                          std::shared_ptr<OccContactEdge> csh2) -> ContactShapeDistance;

}  // namespace siconos::mechanics::occ
#endif
