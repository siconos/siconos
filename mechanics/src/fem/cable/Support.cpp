// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "Support.h"

#include "CableTools.h"  // for sgn function
#include "NewtonImpactFrictionNSL.hpp"
#include "Rope.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::fem::cable::Support::Support(const siconos::algebra::SiconosVector3& coordinates,
                                      double radius)
    : center_pos_{coordinates}, radius_{radius} {
  m_pc2 = std::make_shared<siconos::algebra::SiconosVector3>();
  m_normal = std::make_shared<siconos::algebra::SiconosVector3>();
  m_tangent = std::make_shared<siconos::algebra::SiconosVector3>();
}

void siconos::fem::cable::Support::prepare(const Rope& a_rope) {
  center_pos_(2) -= siconos::fem::cable::tools::sgn(a_rope.supportReaction()(2)) * radius_;
}

void siconos::fem::cable::Support::compute(const siconos::algebra::SiconosVector3& a_p,
                                           double a_tol, double& g,
                                           Eigen::Ref<siconos::algebra::SiconosVector3> G,
                                           Eigen::Ref<siconos::algebra::SiconosVector3> T,
                                           int& c) {
  c = isContact(a_tol, a_p(0) - center_pos_(0), 0, a_p(2) - center_pos_(2), g, G(0), G(1),
                G(2), T(0), T(1), T(2))
          ? 1
          : 0;
}

bool siconos::fem::cable::Support::isContact(
    const Eigen::Ref<siconos::algebra::SiconosVector3>& a_p, double a_tol) {
  // Roller support is assumed to be in the x-z plan
  double dx = a_p(0) - center_pos_(0);
  double dz = a_p(2) - center_pos_(2);
  double d = sqrt(dx * dx + dz * dz);
  double go = d - radius_;
  bool isCt = (go <= a_tol);
  if (isCt) {
    (*m_normal)(0) = dx / d;
    (*m_normal)(2) = dz / d;
    (*m_tangent)(0) = -(*m_normal)(2);
    (*m_tangent)(2) = (*m_normal)(0);
  }
  return isCt;
}

void siconos::fem::cable::Support::InitFriction(double a_mu) {
  if (!m_nslaw)
    m_nslaw = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(0.0, 0.0, a_mu, 2);
}

bool siconos::fem::cable::Support::isContact(double a_tol, double dx, double dy, double dz,
                                             double& g, double& nx, double& ny, double& nz,
                                             double& tx, double& ty, double& tz) {
  double d = sqrt(dx * dx + dy * dy + dz * dz);
  double go = d - radius_;
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

void siconos::fem::cable::Support::display() const {
  std::cout << "--- Support: \n";
  std::cout << "- Center position: ";
  siconos::algebra::print(center_pos_.transpose());
  std::cout << "- Radius: " << center_pos_ << ", " << radius_ << "\n";
  std::cout << "  Normal and tangent vectors at contact: \n";
  if (m_normal) {
    siconos::algebra::print(m_normal->transpose());
    siconos::algebra::print(m_tangent->transpose());
  }
}
