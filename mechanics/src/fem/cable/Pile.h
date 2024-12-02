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

/*! \file Pile.h

  Description of a pylon (position, geometry ...)

*/
#pragma once

#include "Point.h"

namespace siconos::fem::cable {
class Pile : public Point {
 public:
  Pile() = default;
  Pile(const Pile &) = default;
  Pile(Pile &&) = default;
  Pile &operator=(const Pile &) = delete;
  Pile &operator=(Pile &&) = default;

  Pile(const Pile &a_pile, bool a_isStation);

  virtual ~Pile();
  const double &get_radius() const;
  // const double &get_dUp() const;
  // const double &get_dDown() const;

  void from_json(const json &j);
  bool isStation() const;
  void transform(bool a_Up);

 private:
  bool m_isStation{false};
  /** rayon de la poulie dans le cas station, rayon courbure sinon */
  double m_radius{0.};
  /** distance poteau - cable montant */
  double m_dUp{0.};
  /* distance poteau - cable descendant */
  double m_dDown{0.};
  /** hauteur du poteau */
  double m_h{0.};
};
}  // namespace siconos::fem::cable
