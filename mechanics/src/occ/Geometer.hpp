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
#ifndef GEOMETER_HPP
#define GEOMETER_HPP

// OpenCASCADE include
#include <Standard_TypeDef.hxx>  // For Standard_Real ...
#include <gp_Dir.hxx>
#include <iostream>  //TMP
#include <limits>
#include <memory>
#include <variant>

#include "ContactShapeDistance.hpp"
#include "OccContactShape.hpp"
#include "OccUtils.hpp"
#include "cadmbtb.hpp"

namespace siconos::mechanics::occ {

struct OccDistanceType {};
struct CadmbtbDistanceType {};
using DistanceCalculator = std::variant<OccDistanceType, CadmbtbDistanceType>;

struct OccContactFace;
struct OccContactEdge;

template <typename DistType>
ContactShapeDistance distanceFaceFace(std::shared_ptr<OccContactFace> csh1, std::shared_ptr<OccContactFace> csh2) {
}

template <typename DistType>
ContactShapeDistance distanceFaceEdge(std::shared_ptr<OccContactFace> csh1, std::shared_ptr<OccContactEdge> csh2) {
}

template <typename DistType>
ContactShapeDistance distanceEdgeEdge(std::shared_ptr<OccContactEdge> csh1, std::shared_ptr<OccContactEdge> csh2) {
  throw "Geometer: Edge-Edge distance unimplemented";
}

template <>
ContactShapeDistance distanceFaceFace<CadmbtbDistanceType>(std::shared_ptr<OccContactFace> csh1,
                                                           std::shared_ptr<OccContactFace> csh2) {
  return cadmbtb::distanceFaceFace(csh1, csh2);
}

template <>
ContactShapeDistance distanceFaceEdge<CadmbtbDistanceType>(std::shared_ptr<OccContactFace> csh1,
                                                           std::shared_ptr<OccContactEdge> csh2) {
  return cadmbtb::distanceFaceEdge(csh1, csh2);
}

template <>
ContactShapeDistance distanceFaceFace<OccDistanceType>(std::shared_ptr<OccContactFace> csh1,
                                                       std::shared_ptr<OccContactFace> csh2) {
  return occ_distanceFaceFace(csh1, csh2);
}

template <>
ContactShapeDistance distanceFaceEdge<OccDistanceType>(std::shared_ptr<OccContactFace> csh1,
                                                       std::shared_ptr<OccContactEdge> csh2) {
  return occ_distanceFaceEdge(csh1, csh2);
}

template <typename DistType>
struct Geometer {
  // std::shared_ptr<ContactShapeDistance> dist{nullptr};

  Geometer() {
    std::cout << "BUILD A MG\n";
    // dist = std::make_shared<ContactShapeDistance>();
    // dist->value = std::numeric_limits<double>::infinity();
  };

  ContactShapeDistance operator()(const OccContactFace& face1, const OccContactFace& face2) {
    std::cout << "ok ........\n";
    return distanceFaceFace<DistType>(face1, face2);
  };

  ContactShapeDistance operator()(const OccContactFace& face1, const OccContactEdge& edge2) {
    auto dist2 = distanceFaceEdge<DistType>(face1, edge2);
    dist2.normal.Reverse();
    return dist2;
  };

  ContactShapeDistance operator()(const OccContactEdge& edge1, const OccContactFace& face2) {
    return distanceFaceEdge<DistType>(face2, edge1);
  };

  ContactShapeDistance operator()(const OccContactEdge& edge1, const OccContactEdge& edge2) {
    return distanceEdgeEdge<DistType>(edge1, edge2);
  };

  ContactShapeDistance operator()(std::shared_ptr<OccContactFace> face1,
                                  std::shared_ptr<OccContactFace> face2) {
    std::cout << "ok ........\n";
    return distanceFaceFace<DistType>(face1, face2);
  };

  ContactShapeDistance operator()(std::shared_ptr<OccContactFace> face1,
                                  std::shared_ptr<OccContactEdge> edge2) {
    auto dist2 = distanceFaceEdge<DistType>(face1, edge2);
    dist2.normal.Reverse();
    return dist2;
  };

  ContactShapeDistance operator()(std::shared_ptr<OccContactEdge> edge1,
                                  std::shared_ptr<OccContactFace> face2) {
    return distanceFaceEdge<DistType>(face2, edge1);
  };

  ContactShapeDistance operator()(std::shared_ptr<OccContactEdge> edge1,
                                  std::shared_ptr<OccContactEdge> edge2) {
    return distanceEdgeEdge<DistType>(edge1, edge2);
  };
};

}  // namespace siconos::mechanics::occ
#endif
