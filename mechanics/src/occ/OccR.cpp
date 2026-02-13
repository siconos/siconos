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
#include "OccR.hpp"

#include "ContactPoint.hpp"
#include "Geometer.hpp"
#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "OccContactShape.hpp"
#include "SiconosVector.hpp"

// #define  DEBUG_MESSAGES
#include "siconos_debug.h"

siconos::mechanics::occ::OccR::OccR(const ContactPoint& contact1, const ContactPoint& contact2,
                                    const DistanceCalculator& distance_calculator)
    : NewtonEuler3DR(), _contact1{contact1}, _contact2{contact2}, dt{distance_calculator} {}

void siconos::mechanics::occ::OccR::computeh(
    const Eigen::Ref<const siconos::algebra::SiconosVector7>&,
    const std::optional<Eigen::Ref<const siconos::algebra::SiconosVector>>&,
    Eigen::Ref<siconos::algebra::SiconosVector> y) {
  DEBUG_BEGIN("siconos::mechanics::occ::OccR::computeh()\n");
  // std::shared_ptr<ContactShapeDistance> distance{nullptr};
  ContactShapeDistance distance{};

  if (std::get_if<OccDistanceType>(&dt)) {
    std::cout << "OCC Case \n";
    distance = std::visit(Geometer<OccDistanceType>{}, _contact1.contactShape,
                          _contact2.contactShape);

  } else if (std::get_if<CadmbtbDistanceType>(&dt)) {
    //(std::is_same<decltype(distance_calculator), CadmbtbDistanceType>::value)
    std::cout << "CAM Case \n";
    distance = std::visit(Geometer<CadmbtbDistanceType>{}, _contact1.contactShape,
                          _contact2.contactShape);
  }

  DEBUG_PRINTF("---->%g P1=(%g, %g, %g) P2=(%g,%g,%g) N=(%g, %g, %g)\n", distance.value,
               distance.point1.X(), distance.point1.Y(), distance.point1.Z(),
               distance.point2.X(), distance.point2.Y(), distance.point2.Z(),
               distance.normal.X(), distance.normal.Y(), distance.normal.Z());

  (*_Pc1)(0) = distance.point1.X() + _offset1 * distance.normal.X();
  (*_Pc1)(1) = distance.point1.Y() + _offset1 * distance.normal.Y();
  (*_Pc1)(2) = distance.point1.Z() + _offset1 * distance.normal.Z();
  (*_Pc2)(0) = distance.point2.X() - _offset2 * distance.normal.X();
  (*_Pc2)(1) = distance.point2.Y() - _offset2 * distance.normal.Y();
  (*_Pc2)(2) = distance.point2.Z() - _offset2 * distance.normal.Z();

  (*_Nc)(0) = distance.normal.X();
  (*_Nc)(1) = distance.normal.Y();
  (*_Nc)(2) = distance.normal.Z();

  distance.value -= (_offset1 + _offset2);

  y(0) = distance.value;

  DEBUG_EXPR(siconos::algebra::print(y););
  DEBUG_EXPR(siconos::algebra::print(*_Nc););
  DEBUG_END("siconos::mechanics::occ::OccR::computeh()\n");
}
