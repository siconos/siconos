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

/*! \file Pulley.h

  Support (pile + contact description) corresponding to the last and first piles (station)

*/

#pragma once

#include "Support.h"

namespace siconos::fem::cable {
class Pulley : public Support {
 public:
  Pulley(const Pile &a_pile);
  virtual ~Pulley() noexcept = default;
  virtual double get_radius() const override { return m_radiusP; };
  const Point &get_center() const override;
  double get_L(const class Ropeway &a_rope) const;
  double get_tension() const { return m_TR; }

  //------------ statique -------------
  virtual void prepare(const Rope &a_rope) override {};
  virtual void prepare(const Pile &a_start, const Pile &a_end, double T) override;

  int compute(int nb, std::vector<Point> &a_q, int q_offset = 0) const;
  virtual void compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c) override;
  //------------ dynamique -------------
  virtual bool isContact(const std::shared_ptr<siconos::algebra::SiconosVector> &a_p,
                         const double &a_tol) override;

 private:
  double m_radiusP{0.};

  /** Tension */
  double m_TR{0.};

  double dy{0.};
};

//------ Export ----------
void to_json(ojson &j, const Pulley &p);

}  // namespace siconos::fem::cable
