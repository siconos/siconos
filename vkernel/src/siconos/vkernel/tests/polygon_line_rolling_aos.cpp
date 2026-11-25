// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

// AOS kernel reference: a regular polygon (N vertices, circumradius
// Rc, m = 1, g = 9.81) on the line y = 0, exact geometry: one interaction
// per vertex, Lagrangian2dVelocityAngularVelocityR (normal, tangent, rolling rows) between the
// vertex and its projection on the line, NewtonImpactRollingFrictionNSL,
// MoreauJeanOSI, RollingFrictionContact (rolling_friction_2d_driver).  At
// t = 0 the polygon rests on a face, its centre above the origin.
// The slope alpha is modelled by tilting gravity: g = g (sin alpha,
// -cos alpha) in the frame of the line.  The reference of
// rolling_friction_polytope_aos.py.
//
// usage: polygon_line_rolling_aos N Rc mu mu_r alpha vx0 vy0 w0 T h [theta] > out.dat
// columns: t x y theta vx vy w
#include <SiconosKernel.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numbers>
#include <vector>

using siconos::algebra::BlockVector;
using siconos::algebra::SiconosVector2;

// a vertex (vx, vy in the frame of the polygon) against the line y = 0
struct vertex_line_r : siconos::modeling::Lagrangian2dVelocityAngularVelocityR {
  double vx, vy;
  vertex_line_r(double x, double y) : vx(x), vy(y) {}
  void update(const BlockVector& q) {
    const double c = std::cos(q(2)), s = std::sin(q(2));
    SiconosVector2 p1, p2, n;
    p1 << q(0) + c * vx - s * vy, q(1) + s * vx + c * vy;
    p2 << p1(0), 0.;
    n << 0., 1.;
    updateContactPoints(p1, p2, n);
  }
  void computeh(const BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) override {
    update(q);
    Lagrangian2dVelocityAngularVelocityR::computeh(q, y);
  }
  void computeContactMatrix(const BlockVector& q) override {
    update(q);
    Lagrangian2dVelocityAngularVelocityR::computeContactMatrix(q);
  }
};

int main(int argc, char* argv[]) {
  if (argc < 11) {
    std::fprintf(stderr, "usage: %s N Rc mu mu_r alpha vx0 vy0 w0 T h [theta]\n", argv[0]);
    return 2;
  }
  const int N = std::atoi(argv[1]);
  const double Rc = std::atof(argv[2]), mu = std::atof(argv[3]), mu_r = std::atof(argv[4]);
  const double alpha = std::atof(argv[5]);
  const double vx0 = std::atof(argv[6]), vy0 = std::atof(argv[7]), w0 = std::atof(argv[8]);
  const double T = std::atof(argv[9]), h = std::atof(argv[10]);
  const double theta = argc > 11 ? std::atof(argv[11]) : 0.5;
  const double pi = std::numbers::pi;
  const double m = 1., g = 9.81;
  // the inertia of a regular polygon about its centre
  const double I = m * Rc * Rc / 6. * (1. + 2. * std::cos(pi / N) * std::cos(pi / N));

  siconos::algebra::SiconosMatrix mass{3, 3};
  mass.setZero();
  mass(0, 0) = m;
  mass(1, 1) = m;
  mass(2, 2) = I;
  // the vertices, the first and the last on the line at t = 0, mirror
  // images of each other so that both start at a zero gap exactly
  std::vector<double> vx(N), vy(N);
  for (int k = 0; k < (N + 1) / 2; ++k) {
    const double phi = -pi / 2. + pi / N + 2. * pi * k / N;
    vx[k] = Rc * std::cos(phi);
    vy[k] = Rc * std::sin(phi);
    vx[N - 1 - k] = -vx[k];
    vy[N - 1 - k] = vy[k];
  }

  siconos::algebra::SiconosVector q0{3}, vel0{3}, weight{3};
  q0 << 0., -vy[0], 0.;
  vel0 << vx0, vy0, w0;
  weight << m * g * std::sin(alpha), -m * g * std::cos(alpha), 0.;

  auto polygon = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, vel0, mass, siconos::algebra::alias_t);
  polygon->setConstantFext(weight, siconos::algebra::alias_t);

  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., T);
  nsds->insertDynamicalSystem(polygon);
  auto nsl =
      std::make_shared<siconos::modeling::NewtonImpactRollingFrictionNSL>(0., 0., mu, mu_r, 3);
  for (int k = 0; k < N; ++k) {
    auto rel = std::make_shared<vertex_line_r>(vx[k], vy[k]);
    nsds->link(std::make_shared<siconos::modeling::Interaction>(nsl, rel), polygon);
  }

  auto osi = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(0., h);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::RollingFrictionContact>(
      3, SICONOS_ROLLING_FRICTION_2D_NSGS);
  osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  auto s = std::make_shared<siconos::simulation::TimeStepping>(nsds, td, osi, osnspb);

  auto& q = polygon->q_read();
  auto& v = polygon->velocity_read();
  auto print = [&](double t) {
    std::printf("%.10e %.15e %.15e %.15e %.15e %.15e %.15e\n", t, q(0), q(1), q(2), v(0), v(1),
                v(2));
  };
  print(0.);
  while (s->hasNextEvent()) {
    s->computeOneStep();
    print(s->nextTime());
    s->nextStep();
  }
  return 0;
}
