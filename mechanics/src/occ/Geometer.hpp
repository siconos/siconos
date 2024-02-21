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
#include <iostream>              //TMP
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
void distanceFaceFace(const OccContactFace& csh1, const OccContactFace& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist) {}

template <typename DistType>
void distanceFaceEdge(const OccContactFace& csh1, const OccContactEdge& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist) {}

template <typename DistType>
void distanceEdgeEdge(const OccContactEdge& csh1, const OccContactEdge& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist) {
  throw "Geometer: Edge-Edge distance unimplemented";
}

template <>
void distanceFaceFace<CadmbtbDistanceType>(const OccContactFace& csh1,
                                           const OccContactFace& csh2, Standard_Real& X1,
                                           Standard_Real& Y1, Standard_Real& Z1,
                                           Standard_Real& X2, Standard_Real& Y2,
                                           Standard_Real& Z2, Standard_Real& nX,
                                           Standard_Real& nY, Standard_Real& nZ,
                                           Standard_Real& MinDist) {
  cadmbtb::distanceFaceFace(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ, MinDist);
}

template <>
void distanceFaceEdge<CadmbtbDistanceType>(const OccContactFace& csh1,
                                           const OccContactEdge& csh2, Standard_Real& X1,
                                           Standard_Real& Y1, Standard_Real& Z1,
                                           Standard_Real& X2, Standard_Real& Y2,
                                           Standard_Real& Z2, Standard_Real& nX,
                                           Standard_Real& nY, Standard_Real& nZ,
                                           Standard_Real& MinDist) {
  cadmbtb::distanceFaceEdge(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ, MinDist);
}

template <>
void distanceFaceFace<OccDistanceType>(const OccContactFace& csh1, const OccContactFace& csh2,
                                       Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                       Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                       Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                       Standard_Real& MinDist) {
  occ_distanceFaceFace(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ, MinDist);
}

template <>
void distanceFaceEdge<OccDistanceType>(const OccContactFace& csh1, const OccContactEdge& csh2,
                                       Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                       Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                       Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                       Standard_Real& MinDist) {
  occ_distanceFaceEdge(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ, MinDist);
}

template <typename DistType>
struct Geometer {
  std::shared_ptr<ContactShapeDistance> dist{nullptr};

  Geometer() {
    std::cout << "BUILD A MG\n";
    dist = std::make_shared<ContactShapeDistance>();
    dist->value = std::numeric_limits<double>::infinity();
  };

  std::shared_ptr<ContactShapeDistance> operator()(const OccContactFace& face1,
                                                   const OccContactFace& face2) {
    std::cout << "ok ........\n";
    distanceFaceFace<DistType>(face1, face2, dist->x1, dist->y1, dist->z1, dist->x2, dist->y2,
                               dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(const OccContactFace& face1,
                                                   const OccContactEdge& edge2) {
    distanceFaceEdge<DistType>(face1, edge2, dist->x1, dist->y1, dist->z1, dist->x2, dist->y2,
                               dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    dist->nx = -dist->nx;
    dist->ny = -dist->ny;
    dist->nz = -dist->nz;
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(const OccContactEdge& edge1,
                                                   const OccContactFace& face2) {
    distanceFaceEdge<DistType>(face2, edge1, dist->x1, dist->y1, dist->z1, dist->x2, dist->y2,
                               dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    return dist;
  };
  std::shared_ptr<ContactShapeDistance> operator()(const OccContactEdge& edge1,
                                                   const OccContactEdge& edge2) {
    distanceEdgeEdge<DistType>(edge1, edge2, dist->x1, dist->y1, dist->z1, dist->x2, dist->y2,
                               dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(std::shared_ptr<OccContactFace> face1,
                                                   std::shared_ptr<OccContactFace> face2) {
    std::cout << "ok ........\n";
    std::cout << "yyyayyazyzyayzyazyzyz " << dist->value << "\n";
    distanceFaceFace<DistType>(*face1, *face2, dist->x1, dist->y1, dist->z1, dist->x2,
                               dist->y2, dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    std::cout << "yyyayyazyzyayzyazyzyz " << dist->value << "\n";
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(std::shared_ptr<OccContactFace> face1,
                                                   std::shared_ptr<OccContactEdge> edge2) {
    distanceFaceEdge<DistType>(*face1, *edge2, dist->x1, dist->y1, dist->z1, dist->x2,
                               dist->y2, dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    dist->nx = -dist->nx;
    dist->ny = -dist->ny;
    dist->nz = -dist->nz;
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(std::shared_ptr<OccContactEdge> edge1,
                                                   std::shared_ptr<OccContactFace> face2) {
    distanceFaceEdge<DistType>(*face2, *edge1, dist->x1, dist->y1, dist->z1, dist->x2,
                               dist->y2, dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    return dist;
  };

  std::shared_ptr<ContactShapeDistance> operator()(std::shared_ptr<OccContactEdge> edge1,
                                                   std::shared_ptr<OccContactEdge> edge2) {
    distanceEdgeEdge<DistType>(*edge1, *edge2, dist->x1, dist->y1, dist->z1, dist->x2,
                               dist->y2, dist->z2, dist->nx, dist->ny, dist->nz, dist->value);
    return dist;
  };
};

}  // namespace siconos::mechanics::occ
#endif
