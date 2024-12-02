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

/*! \file Support.h

  A pile and its contact point description (nonsmooth law ...)
*/
#pragma once

#include "Point.h"

namespace siconos::algebra {
class SiconosVector;
}
namespace siconos::modeling {

class NonSmoothLaw;
}

namespace siconos::fem::cable {

class Pile;
class Point;
class Rope;

class Support {
 protected:
  const Pile &r_pile;  // reference to the cable model
  Point m_p;

  std::shared_ptr<siconos::modeling::NonSmoothLaw> m_nslaw{nullptr};

  std::shared_ptr<siconos::algebra::SiconosVector> m_pc2{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> m_normal{nullptr};
  std::shared_ptr<siconos::algebra::SiconosVector> m_tangent{nullptr};

  // Rule of 5
  Support() = delete;
  Support(const Support &) = delete;
  Support(Support &&) = delete;
  Support &operator=(const Support &) = delete;
  Support &operator=(Support &&) = delete;

 public:
  Support(const Pile &a_pile);
  virtual ~Support() noexcept = default;
  virtual double get_radius() const;
  virtual const Point &get_center() const { return m_p; };

  //------------ statique -------------
  virtual void prepare(const Rope &a_rope);
  virtual void prepare(const Pile &a_start, const Pile &a_end, double T);

  virtual void compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c);

  //------------ dynamique -------------
  virtual bool isContact(const std::shared_ptr<siconos::algebra::SiconosVector> &a_p,
                         const double &a_tol);

  virtual void InitFriction(double a_mu);
  std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw();

  inline std::shared_ptr<siconos::algebra::SiconosVector> pc2() const { return m_pc2; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> normal() const { return m_normal; }
  inline std::shared_ptr<siconos::algebra::SiconosVector> tangent() const { return m_tangent; }

  bool isContact(const double &a_tol, const double &dx, const double &dy, const double &dz,
                 double &g,

                 double &nx, double &ny, double &nz,

                 double &tx, double &ty, double &tz);
};

//------ Export ----------
void to_json(ojson &j, const Support &p);

}  // namespace siconos::fem::cable
