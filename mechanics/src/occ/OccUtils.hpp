/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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
