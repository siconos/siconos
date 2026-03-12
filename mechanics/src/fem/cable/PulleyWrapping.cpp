#include "PulleyWrapping.h"

#include <numbers>

#include "Pylon.h"
#include "Rope.h"
#include "Ropeway.h"
#include "SiconosVector.hpp"

void siconos::fem::cable::PulleyWrapping::prepare(const Pylon &a_start, const Pylon &a_end,
                                                  double T) {
  const auto &coords_s = a_start.coords();
  const auto &coords_e = a_end.coords();

  radius_ = 0.5 * (coords_s - coords_e).norm();
  center_pos_ = 0.5 * (coords_s + coords_e);
  theta_start_ = std::numbers::pi * 0.5;
  if (coords_s(1) - coords_e(1) <= 0) {
    theta_start_ *= -1.;
  }

  tension_ = T;
}

int siconos::fem::cable::PulleyWrapping::computeMesh(
    int number_of_nodes_, siconos::algebra::SiconosVector &positions, int q_offset) const {
  int n = positions.size();
  for (int i = 0; i < number_of_nodes_ && 3 * (i + q_offset) < n; i++) {
    double theta = theta_start_ + (std::numbers::pi / (number_of_nodes_ - 1)) * i;
    int pos = 3 * (i + q_offset);
    positions(pos) = center_pos_(0) + radius_ * std::cos(theta);
    positions(pos + 1) = center_pos_(1) + radius_ * std::sin(theta);
    positions(pos + 2) = center_pos_(2);
  }
  return q_offset + number_of_nodes_ - 1;  // TO BE CHECKED !!
}

void siconos::fem::cable::PulleyWrapping::compute(
    const siconos::algebra::SiconosVector3 &a_p, double a_tol, double &g,
    Eigen::Ref<siconos::algebra::SiconosVector3> G,
    Eigen::Ref<siconos::algebra::SiconosVector3> T, int &c) {
  c = Support::isContact(a_tol, a_p(0) - center_pos_(0), a_p(1) - center_pos_(1), 0, g, G(0),
                         G(1), G(2), T(0), T(1), T(2))
          ? 1
          : 0;
}

double siconos::fem::cable::PulleyWrapping::length(const Ropeway &a_rope) const {
  return radius_ * std::numbers::pi /
         (1 + tension_ / a_rope.mechanicalProperties0().crossSectionRigidity());
}

bool siconos::fem::cable::PulleyWrapping::isContact(
    const Eigen::Ref<siconos::algebra::SiconosVector3> &a_p, double a_tol) {
  // PulleyWrapping support is assumed to be in the x-y plan
  double dx = a_p(0) - center_pos_(0);
  double dy = a_p(1) - center_pos_(1);
  double d = sqrt(dx * dx + dy * dy);
  double go = d - radius_;
  bool isCt = (go <= a_tol);
  if (isCt) {
    (*m_normal)(0) = dx / d;
    (*m_normal)(1) = dy / d;
    (*m_tangent)(0) = -(*m_normal)(1);
    (*m_tangent)(1) = (*m_normal)(0);
  }
  return isCt;
}

void siconos::fem::cable::PulleyWrapping::display() const {
  std::cout << " Pulley type support \n";
  std::cout << " Tension: " << tension_ << "\n";
  Support::display();
}