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
#include "MBTB_Body.hpp"

#include "StorageTools.hpp"

siconos::mechanisms::MBTB_Body::MBTB_Body(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0, double mass,
    Eigen::Ref<siconos::algebra::SiconosMatrix> I,
    Eigen::Ref<siconos::algebra::SiconosVector> centerOfMass, const std::string& BodyName,
    const std::string& CADFile)
    : NewtonEulerDS{q0, v0, mass, I, siconos::algebra::alias_t},
      _mBodyName{BodyName},
      _cadFileName{CADFile} {
  centerOfMass_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(centerOfMass.data(), ndof_);
}
