#include "Pulley.h"

#include <numbers>

#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothLaw.hpp"
#include "Pylon.h"
#include "Rope.h"
#include "Ropeway.h"
#include "SiconosVector.hpp"

siconos::fem::cable::Pulley::Pulley(const Pylon &a_pile) : Support(a_pile) {}

void siconos::fem::cable::Pulley::prepare(const Pylon &a_start, const Pylon &a_end, double T) {
  Point delta;
  delta.diff(a_start, a_end);
  m_radiusP = delta.norm() / 2.0;
  m_p.add(a_start, a_end);
  m_p.mult(0.5);
  dy = a_start.y - m_p.y;
  m_TR = T;
}

int siconos::fem::cable::Pulley::compute(int nb, std::vector<Point> &a_q, int q_offset) const {
  double pi2 = std::numbers::pi * 0.5;
  if (dy <= 0) {
    pi2 = -pi2;
  }
  int n = a_q.size();
  for (int i = 0; i < nb && (i + q_offset) < n; i++) {
    double theta = pi2 + (std::numbers::pi / (nb - 1)) * i;

    Point &p = a_q[i + q_offset];
    p.x = m_p.x + m_radiusP * std::cos(theta);
    p.y = m_p.y + m_radiusP * std::sin(theta);
    p.z = m_p.z;
  }
  return q_offset + nb - 1;
}

void siconos::fem::cable::Pulley::compute(const Point &a_p, double a_tol, double &g, Point &G,
                                          Point &T, int &c) {
  c = Support::isContact(a_tol, a_p.x - m_p.x, a_p.y - m_p.y, 0, g, G.x, G.y, G.z, T.x, T.y,
                         T.z)
          ? 1
          : 0;
}

const siconos::fem::cable::Point &siconos::fem::cable::Pulley::get_center() const {
  return m_p;
}

double siconos::fem::cable::Pulley::get_L(const Ropeway &a_rope) const {
  return m_radiusP * std::numbers::pi / (1 + m_TR / a_rope.get_meca0().get_EA());
}

void siconos::fem::cable::to_json(ojson &j, const Pulley &p) {
  j["radius"] = p.get_radius();
  j["center"] = p.get_center();
  j["tension"] = p.get_tension();
}

bool siconos::fem::cable::Pulley::isContact(
    const Eigen::Ref<siconos::algebra::SiconosVector3> &a_p, double a_tol) {
  double dx = a_p(0) - m_p.x;
  double dy = a_p(1) - m_p.y;
  double d = sqrt(dx * dx + dy * dy);
  double go = d - get_radius();
  bool isCt = (go <= a_tol);
  if (isCt) {
    (*m_normal)(0) = dx / d;
    (*m_normal)(1) = dy / d;
    (*m_tangent)(0) = -(*m_normal)(1);
    (*m_tangent)(1) = (*m_normal)(0);
  }
  return isCt;
}
