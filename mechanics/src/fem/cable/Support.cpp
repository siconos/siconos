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

#include "Support.h"

#include "CableTools.h"  // for sgn function
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothLaw.hpp"
#include "Pylon.h"
#include "Rope.h"
#include "SiconosVector.hpp"

siconos::fem::cable::Support::Support(const Pylon &a_pile) : pylon_(a_pile) {
  m_p = pylon_;

  m_pc2 = std::make_shared<siconos::algebra::SiconosVector3>();
  m_normal = std::make_shared<siconos::algebra::SiconosVector3>();
  m_tangent = std::make_shared<siconos::algebra::SiconosVector3>();
}

double siconos::fem::cable::Support::get_radius() const { return pylon_.get_radius(); }

void siconos::fem::cable::Support::prepare(const Rope &a_rope) {
  double radius = siconos::fem::cable::tools::sgn(a_rope.supportReaction().z) * get_radius();
  m_p.z -= radius;
}

void siconos::fem::cable::Support::prepare(const Pylon &a_start, const Pylon &a_end,
                                           double T) {
  // Does nothing ...
}

void siconos::fem::cable::Support::compute(const Point &a_p, double a_tol, double &g, Point &G,
                                           Point &T, int &c) {
  c = isContact(a_tol, a_p.x - m_p.x, 0, a_p.z - m_p.z, g, G.x, G.y, G.z, T.x, T.y, T.z) ? 1
                                                                                         : 0;
}

void siconos::fem::cable::to_json(ojson &j, const Support &s) {
  j["radius"] = s.get_radius();
  j["p"] = s.get_center();
}

bool siconos::fem::cable::Support::isContact(
    const Eigen::Ref<siconos::algebra::SiconosVector3> &a_p, double a_tol) {
  double dx = a_p.getValue(0) - m_p.x;
  double dz = a_p.getValue(2) - m_p.z;
  double d = sqrt(dx * dx + dz * dz);
  double go = d - get_radius();
  bool isCt = (go <= a_tol);
  if (isCt) {
    m_normal->setValue(0, dx / d);
    m_normal->setValue(2, dz / d);
    m_tangent->setValue(0, -m_normal->getValue(2));
    m_tangent->setValue(2, m_normal->getValue(0));
  }
  return isCt;
}

void siconos::fem::cable::Support::InitFriction(double a_mu) {
  if (m_nslaw != nullptr) {
    // m_nslaw->
  } else {
    m_nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0.0, 0.0, a_mu, 2);
  }
}

std::shared_ptr<siconos::modeling::NonSmoothLaw> siconos::fem::cable::Support::nslaw() {
  return m_nslaw;
}

bool siconos::fem::cable::Support::isContact(double a_tol, double dx, double dy, double dz,
                                             double &g, double &nx, double &ny, double &nz,
                                             double &tx, double &ty, double &tz) {
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  double go = d - get_radius();
  bool isCt = (go <= a_tol);
  if (isCt) {
    g = go;
    nx = dx / d;
    ny = dy / d;
    nz = dz / d;
    tx = -ny - nz;
    ty = (ny) ? nx : 0.0;
    tz = (nz) ? nx : 0.0;
  }
  return isCt;
}
