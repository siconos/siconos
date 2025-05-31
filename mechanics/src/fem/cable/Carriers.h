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
  Class Carriers, description of the vehicles

*/
#pragma once
#include "BaseModel.h"

namespace siconos::fem::cable {

/** Vehicles description (all of them in one single object)
 *
 * Build from json (in TransportCableModel) and read-only use afterwards
 */
class Carriers : public BaseModel {
 public:
  Carriers() = default;
  Carriers(const Carriers &) = default;

  virtual ~Carriers() noexcept = default;

  /** \return the number of vehicles  */
  auto get_number_of_vehicles() const { return m_n; }

  /** \return density  */
  auto get_rho() const { return m_rho; }

  /** \return distance between two vehicles  */
  auto get_d_inter_vehicules() const { return m_d; }

  /** \return distance between first pile and the first vehicule */
  auto get_d_start() const { return m_d_start; }

  /** \return  up load  */
  auto up_load() const { return m_up_load; }

  /** \return down load  */
  auto down_load() const { return m_down_load; }

  /** \return mass 100% of one vehicule (kg)  */
  auto loaded_mass() const { return m_loaded_mass; }

 private:
  Carriers(Carriers &&) = delete;
  Carriers &operator=(const Carriers &) = delete;
  Carriers &operator=(Carriers &&) = delete;

  void from_json(const json &j);

  /** number of vehicules*/
  int m_n{0};

  /** mass of one vehicle (kg) */
  double m_mass{0.};

  /** distance between two vehicles (m)*/
  double m_d{0.};  //

  double m_rho{0.};

  /** mass 100% of one vehicule (kg) */
  double m_loaded_mass{0.};

  double m_up_load{1.};    // % up load
  double m_down_load{0.};  // % down load
  double m_d_start{0.};    // % (m_d), distance of the first vehicule
};
}  // namespace siconos::fem::cable
