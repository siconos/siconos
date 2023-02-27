/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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
  Description of the vehicles

*/
#pragma once
#include "BaseModel.h"

namespace siconos::fem::cable {
// punctual masses
class Carriers : public BaseModel {
 public:
  Carriers() = default;
  Carriers(const Carriers &) = default;

  virtual ~Carriers() noexcept = default;
  const double &get_rho() const;
  const double &get_d_inter_vehicules() const;

 private:
  Carriers(Carriers &&) = delete;
  Carriers &operator=(const Carriers &) = delete;
  Carriers &operator=(Carriers &&) = delete;

  void from_json(const json &j);

  /** number of vehicules*/
  int m_n{0};  //
  /** mass of one vehicle (kg) */
  double m_mass{0.};
  /** distance between two vehicles (m)*/
  double m_d{0.};  //

  double m_rho{0.};
};
}  // namespace siconos::fem::cable
