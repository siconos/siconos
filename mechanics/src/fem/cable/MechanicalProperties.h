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
#include "BaseModel.h"

namespace siconos::fem::cable {

/** Class to handle mechanical properties of a cable */
class MechanicalProperties : public BaseModel {
 public:
  /** Default and only constructor */
  MechanicalProperties() = default;

  MechanicalProperties &operator=(const MechanicalProperties &) = default;
  virtual ~MechanicalProperties() noexcept = default;

  inline auto get_EA() const { return m_EA; }
  inline auto get_rho() const { return m_rho; }
  inline auto get_T0() const { return m_T0; }
  inline auto get_alpha() const { return 9.81 * m_rho / m_T0; }
  inline auto get_beta() const { return m_T0 / m_EA; }

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
