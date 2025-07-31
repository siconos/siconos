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

/*! \file MechanicalProperties.h
  \brief Some mechanical Properties of a cable
*/

#pragma once
#include <nlohmann/json.hpp>

namespace siconos::fem::cable {

/** Class to handle mechanical properties of a cable */
class MechanicalProperties {
 private:
  /** rigidity (N) */
  double crossSectionRigidity_{1.};

  /** linear density	(kg/m) */
  double linearDensity_{0.};

  /** initial tension (kN)*/
  double initialTension_{1.};

 public:
  /** Default and only constructor */
  MechanicalProperties() = default;

  explicit MechanicalProperties(const nlohmann::json &j)
      : crossSectionRigidity_(j.value("EA", 1.0)),
        linearDensity_(j.value("linearDensity", 0.0)),
        initialTension_(j.value("T0", 1.0)) {}

  MechanicalProperties &operator=(const MechanicalProperties &) = default;
  ~MechanicalProperties() noexcept = default;

  inline auto crossSectionRigidity() const { return crossSectionRigidity_; }
  inline auto linearDensity() const { return linearDensity_; }
  inline auto initialTension() const { return initialTension_; }
  inline auto get_alpha() const { return 9.81 * linearDensity_ / initialTension_; }
  inline auto get_beta() const { return initialTension_ / crossSectionRigidity_; }

  void set_T(double a_T);
  void set_rho(double a_rho);

  /** Display parameters to screen */
  void display() const;
};
}  // namespace siconos::fem::cable

// Serialization
namespace nlohmann {
template <>
struct adl_serializer<siconos::fem::cable::MechanicalProperties> {
  static siconos::fem::cable::MechanicalProperties from_json(const json &j) {
    return siconos::fem::cable::MechanicalProperties(j);
  }

  static void to_json(json &j, const siconos::fem::cable::MechanicalProperties &input) {
    j = json{{"EA", input.crossSectionRigidity()},
             {"linearDensity", input.linearDensity()},
             {"T0", input.initialTension()}};
  }
};
}  // namespace nlohmann