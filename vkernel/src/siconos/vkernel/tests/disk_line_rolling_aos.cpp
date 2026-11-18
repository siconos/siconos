/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
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

// AOS kernel reference: one disk on the line y = 0, exact
// geometry, Lagrangian2d3DR (normal, tangent, rolling rows),
// NewtonImpactRollingFrictionNSL, MoreauJeanOSI, RollingFrictionContact
// (rolling_friction_2d_driver).  The slope alpha is modelled by tilting
// gravity: g = g (sin alpha, -cos alpha) in the frame of the line.
// The reference of rolling_friction_aos.py (R = 1, m = 1,
// I = 1/2, g = 9.81).
//
// usage: disk_line_rolling_aos mu mu_r alpha v0 w0 T h [theta] > out.dat
// columns: t x y theta vx vy w rn rt rr
#include <SiconosKernel.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using siconos::algebra::BlockVector;
using siconos::algebra::SiconosVector2;

struct disk_line_r : siconos::modeling::Lagrangian2d3DR {
  double R;
  explicit disk_line_r(double r) : R(r) {}
  void update(const BlockVector& q) {
    SiconosVector2 p1, p2, n;
    p1 << q(0), q(1) - R;
    p2 << q(0), 0.;
    n << 0., 1.;
    updateContactPoints(p1, p2, n);
  }
  void computeh(const BlockVector& q, Eigen::Ref<siconos::algebra::SiconosVector> y) override {
    update(q);
    Lagrangian2d3DR::computeh(q, y);
  }
  void computeJacobianhOver_q(const BlockVector& q) override {
    update(q);
    Lagrangian2d3DR::computeJacobianhOver_q(q);
  }
};

int main(int argc, char* argv[]) {
  if (argc < 8) {
    std::fprintf(stderr, "usage: %s mu mu_r alpha v0 w0 T h [theta]\n", argv[0]);
    return 2;
  }
  double mu = std::atof(argv[1]), mu_r = std::atof(argv[2]), alpha = std::atof(argv[3]);
  double v0 = std::atof(argv[4]), w0 = std::atof(argv[5]);
  double T = std::atof(argv[6]), h = std::atof(argv[7]);
  double theta = argc > 8 ? std::atof(argv[8]) : 0.5;
  double R = 1., m = 1., g = 9.81, I = 0.5 * m * R * R;

  siconos::algebra::SiconosMatrix mass{3, 3};
  mass.setZero();
  mass(0, 0) = m;
  mass(1, 1) = m;
  mass(2, 2) = I;
  siconos::algebra::SiconosVector q0{3}, vel0{3}, weight{3};
  q0 << 0., R, 0.;
  vel0 << v0, 0., w0;
  weight << m * g * std::sin(alpha), -m * g * std::cos(alpha), 0.;

  auto disk = std::make_shared<siconos::modeling::LagrangianLinearTIDS>(
      q0, vel0, mass, siconos::algebra::alias_t);
  disk->setConstantFext(weight, siconos::algebra::alias_t);

  auto rel = std::make_shared<disk_line_r>(R);
  auto nsl =
      std::make_shared<siconos::modeling::NewtonImpactRollingFrictionNSL>(0., 0., mu, mu_r, 3);
  auto inter = std::make_shared<siconos::modeling::Interaction>(nsl, rel);

  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., T);
  nsds->insertDynamicalSystem(disk);
  nsds->link(inter, disk);

  auto osi = std::make_shared<siconos::integrators::MoreauJeanOSI>(theta);
  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(0., h);
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::RollingFrictionContact>(
      3, SICONOS_ROLLING_FRICTION_2D_NSGS);
  osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-14;
  auto s = std::make_shared<siconos::simulation::TimeStepping>(nsds, td, osi, osnspb);

  auto& q = disk->q_read();
  auto& v = disk->velocity_read();
  auto print = [&](double t) {
    auto lambda = inter->lambda(1);
    std::printf("%.10e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", t, q(0), q(1),
                q(2), v(0), v(1), v(2), (*lambda)(0), (*lambda)(1), (*lambda)(2));
  };
  print(0.);
  while (s->hasNextEvent()) {
    s->computeOneStep();
    print(s->nextTime());
    s->nextStep();
  }
  return 0;
}
