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

#include "MechanicalProperties.h"

void siconos::fem::cable::MechanicalProperties::set_T(double a_T) {
  if (a_T > 0.0) m_T0 = a_T;
}

void siconos::fem::cable::MechanicalProperties::set_rho(double a_rho) {
  if (a_rho > 0.0) m_rho = a_rho;
}

void siconos::fem::cable::MechanicalProperties::from_json(const json &j) {
  j.at("EA").get_to(m_EA);
  j.at("rho").get_to(m_rho);
  j.at("T0").get_to(m_T0);
}
