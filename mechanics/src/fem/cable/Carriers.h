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

/*! \file Carriers.h
  The set carriers (cabins)

*/
#pragma once
#include <nlohmann/json.hpp>

namespace siconos::fem::cable {

/** Vehicles description (all of them in one single object)
 *
 * Build from json (in TransportCableModel) and read-only use afterwards
 */
struct Carriers {
  Carriers() = delete;
  Carriers(Carriers &&) = delete;
  Carriers(const Carriers &) = delete;
  Carriers &operator=(const Carriers &) = delete;
  Carriers &operator=(Carriers &&) = delete;
  virtual ~Carriers() noexcept = default;

  explicit Carriers(const nlohmann::json &j)
      : number(j.value("number", 0)),
        mass(j.value("mass", 0.0)),
        distanceBetweenCarriers(j.value("distanceBetweenCarriers", 0.0)),
        linearDensity(j.value("linearDensity", 0.0)),
        firstCarrierInitialPosition(j.value("firstCarrierInitialPosition", 0.0)) {}

  /** number of carriers*/
  const int number;

  /** mass of an empty vehicle (kg)*/
  const double mass;

  /** distance between two vehicles (m)*/
  const double distanceBetweenCarriers;  //

  /** fictitious linear density due to carriers */
  const double linearDensity;

  // /** mass (max) of a loaded carrier (kg) */
  // const double loaded_mass_{0.};

  /** Distance between starting pylon and the first carrier */
  const double firstCarrierInitialPosition;
};
}  // namespace siconos::fem::cable

// Serialization
namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::Carriers> {
  static siconos::fem::cable::Carriers from_json(const json &j) {
    return siconos::fem::cable::Carriers(j);
  }

  static void to_json(ordered_json &j, const siconos::fem::cable::Carriers &c) {
    j = json{{"number", c.number},
             {"mass", c.mass},
             {"distanceBetweenCarriers", c.distanceBetweenCarriers},
             {"linearDensity", c.linearDensity},
             {"firstCarrierInitialPosition", c.firstCarrierInitialPosition}};
  }
};
}  // namespace nlohmann
