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
Support class

*/
#pragma once

#include "Point.h"
#include "SiconosVector.hpp"

namespace siconos::modeling {

class NonSmoothLaw;
}

namespace siconos::fem::cable {

class Pylon;
class Point;
class Rope;

/** A pylon, its contact point description (coordinates, normal ...) and a nonsmooth law
 *
 */
class Support {
 protected:
  /** The pylon (geometrical descr) */
  const Pylon &pylon_;

  /**  Pylon position */
  Point m_p;

  /** Non smooth law (for the contact pylon-rope) */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> m_nslaw{nullptr};

  /** Contact point coordinates  */
  std::shared_ptr<siconos::algebra::SiconosVector3> m_pc2{nullptr};

  /** Normal at contact */
  std::shared_ptr<siconos::algebra::SiconosVector3> m_normal{nullptr};

  /** Tangent at contact */
  std::shared_ptr<siconos::algebra::SiconosVector3> m_tangent{nullptr};

  // Rule of 5
  Support() = delete;
  Support(const Support &) = delete;
  Support(Support &&) = delete;
  Support &operator=(const Support &) = delete;
  Support &operator=(Support &&) = delete;

 public:
  Support(const Pylon &a_pile);
  virtual ~Support() noexcept = default;
  virtual double get_radius() const;
  virtual const Point &get_center() const { return m_p; };

  //------------ static -------------
  virtual void prepare(const Rope &a_rope);
  virtual void prepare(const Pylon &a_start, const Pylon &a_end, double T);

  virtual void compute(const Point &a_p, double a_tol, double &g, Point &G, Point &T, int &c);

  //------------ dynamic -------------

  /** \return true if contact is on */
  virtual bool isContact(const Eigen::Ref<siconos::algebra::SiconosVector3> &a_p,
                         double a_tol);

  /** \return true if contact is on */
  bool isContact(double a_tol, double dx, double dy, double dz, double &g, double &nx,
                 double &ny, double &nz, double &tx, double &ty, double &tz);

  /** Build the nonsmooth law (if required) */
  virtual void InitFriction(double a_mu);

  /** \return the nonsmooth law used for this support */
  std::shared_ptr<siconos::modeling::NonSmoothLaw> nslaw();

  /** \return the coordinates of the contact point */
  inline auto pc2() const { return m_pc2; }

  /** \return the normal at contact */
  inline auto normal() const { return m_normal; }

  /**\return the tangent at contact */
  inline auto tangent() const { return m_tangent; }
};

//------ Export ----------
void to_json(ojson &j, const Support &p);

}  // namespace siconos::fem::cable
