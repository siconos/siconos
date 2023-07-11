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

/*! \file Cable.h
  \brief Material properties of a cable
*/

#pragma once
#include "BaseModel.h"

namespace siconos::fem::cable {

class Cable : public BaseModel {
 public:
  Cable() = default;
  Cable &operator=(const Cable &) = default;
  virtual ~Cable() noexcept = default;

  const double &get_EA() const;
  const double &get_rho() const;
  const double &get_T0() const;

  double get_alpha() const;
  double get_beta() const;

  void set_T(double a_T);
  void set_rho(double a_rho);

 private:
  void from_json(const json &j);

  /** rigidity (N) */
  double m_EA{1.};

  /** linear density	(kg/m) */
  double m_rho{0.};

  /** initial tension (kN)*/
  double m_T0{1.};
};
}  // namespace siconos::fem::cable
