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

/*! \file Pylon.h
  Pylon class
*/
#pragma once

#include "Point.h"

namespace siconos::fem::cable {
/**   Description of a pylon (position, geometry ...)
 */
class Pylon : public Point {
 public:
  Pylon() = default;
  Pylon(const Pylon &) = default;
  Pylon(Pylon &&) = default;
  Pylon &operator=(const Pylon &) = delete;
  Pylon &operator=(Pylon &&) = default;

  Pylon(const Pylon &a_pile, bool a_isStation);

  virtual ~Pylon();
  inline auto get_radius() const { return m_radius; }
  // const double &get_dUp() const;
  // const double &get_dDown() const;

  void from_json(const json &j);
  bool isStation() const;
  void transform(bool a_Up);

 private:
  /** True if it's a station pylon else false */
  bool m_isStation{false};

  /** radius (of the pulley for a station, of curvature in other cases) */
  double m_radius{0.};

  /** distance between the pylon and the up-rope */
  double m_dUp{0.};

  /* distance between the pylon and the down-rope  */
  double m_dDown{0.};

  /** pylon height */
  double m_h{0.};
};
}  // namespace siconos::fem::cable
