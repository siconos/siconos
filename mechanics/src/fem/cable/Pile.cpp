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

#include "Pile.h"

siconos::fem::cable::Pile::Pile(const Pile &a_pile, bool a_isStation) : Point(a_pile) {
  m_radius = a_pile.get_radius();
  m_h = 0;
  m_dDown = 0;
  m_dUp = 0;
  m_isStation = a_isStation;
}

siconos::fem::cable::Pile::~Pile() {}

const double &siconos::fem::cable::Pile::get_radius() const { return m_radius; }

void siconos::fem::cable::Pile::from_json(const json &j) {
  Point::from_json(j);
  j.at("R").get_to(m_radius);
  try {
    j.at("dUp").get_to(m_dUp);
    j.at("dDown").get_to(m_dDown);
    j.at("h").get_to(m_h);
  } catch (const std::exception &) {
  }
}

bool siconos::fem::cable::Pile::isStation() const { return m_isStation; }

void siconos::fem::cable::Pile::transform(bool a_Up) {
  if (!a_Up) {
    if (m_isStation) {
      y += 2.0 * m_radius;
    } else {
      y += m_dUp + m_dDown;
    }
  }
}
