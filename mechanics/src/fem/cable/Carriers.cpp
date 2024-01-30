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
#include "Carriers.h"

const double &siconos::fem::cable::Carriers::get_rho() const { return m_rho; }

const double &siconos::fem::cable::Carriers::get_d_inter_vehicules() const { return m_d; }

const double &siconos::fem::cable::Carriers::get_d_start() const { return m_d_start; }

void siconos::fem::cable::Carriers::from_json(const json &j) {
  j.at("rho").get_to(m_rho);
  j.at("mass").get_to(m_mass);
  j.at("d").get_to(m_d);

  if (j.contains("d_start")) {
    j.at("d_start").get_to(m_d_start);
  }
}
