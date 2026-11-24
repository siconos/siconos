/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <print>
#include <vector>

#include "siconos/vkernel.hpp"
#include "siconos/vkernel/collision/space_filter.hpp"
#include "siconos/vkernel/storage/registry.hpp"

// clang-format off
#include "siconos/vkernel/collision/polytopepolytope_r.hpp"
#include "siconos/vkernel/collision/polytopefsegment_r.hpp"
// clang-format on

// The two rigid examples of Acary & Collins-Craft (JTCAM 2025), section 5
// (h = 1e-4, theta = 1/2, e = 1, g = 10, m = 1): a polytope under gravity
// on a fixed ground segment, with Newton and with Frémond.
//
// - the impacting stick (l = 1, I = 1/12, mu = 0.01, contact at its lower
//   end): at the first impact the slip reverses, and Coulomb on u(k + 1)
//   makes friction do positive work, E grows although e = 1;
// - the rocking block (1 x 1, I = 1/6, mu = 0.1, two corner contacts):
//   Newton creates energy at its second impact, Frémond never does.
//
// With theta = 1/2, K = C = 0 and gravity in E, the discrete energy balance
// of the scheme is E(k + 1) - E(k) = u(k + theta)^T p(k + 1), the contact
// work, and its normal part is zero for e = 1. The reference values come
// from the same scheme with an exact local solver.

#define SICONOS_FREMOND_RIGID_CONFIG(NS, LAW)                                               \
  namespace siconos::vkernel::NS {                                                          \
  struct body : model::lagrangian_ds {};                                                    \
  struct nslaw : model::LAW {};                                                             \
  using polytope_shape = collision::shape::polytope;                                        \
  using segment_shape = collision::shape::segment;                                          \
  using diskdisk_r = collision::diskdisk_r;                                                 \
  using diskfsegment_r = collision::diskfsegment_r;                                         \
  struct topo;                                                                              \
  using polytopepolytope_r = collision::polytopepolytope_r<topo>;                           \
  using polytopefsegment_r = collision::polytopefsegment_r<topo>;                           \
  using cxd_narrowphase = collision::cxd_narrowphase<topo>;                                 \
  struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r,                \
                                          polytopepolytope_r, polytopefsegment_r> {};       \
  struct topo : simul::topology<body, interaction> {};                                      \
  struct osi : simul::one_step_integrator<topo>::moreau_jean {};                            \
  struct td : simul::time_discretization<> {};                                              \
  struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};                        \
  struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};                               \
  struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};                       \
  template <typename T>                                                                     \
  struct env : standard_environment<T> {                                                    \
    using params = config::map<config::iparam<"dof", 3>, config::iparam<"ncgroups", 1>>;    \
  };                                                                                        \
  struct data_t : storage::make_runtime_database<                                           \
                      env, simulation, body, polytope_shape, segment_shape, diskdisk_r,     \
                      diskfsegment_r, nslaw, interaction, topo, polytopepolytope_r,         \
                      polytopefsegment_r, cxd_narrowphase, osi, td, osnspb,                 \
                      storage::with_properties<                                             \
                          storage::attached<body, storage::pattern::symbol<"shape">,        \
                                            storage::some::item_ref<polytope_shape>>,       \
                          storage::time_invariant<storage::pattern::attr_t<body, "fext">>,  \
                          storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>, \
                          storage::assembled_diagonal<                                      \
                              storage::pattern::attr_t<osi, "mass_matrix_assembled">>>> {}; \
  /* the types of the config, as one template argument */                                   \
  struct types {                                                                            \
    using data_t = siconos::vkernel::NS::data_t;                                            \
    using body = siconos::vkernel::NS::body;                                                \
    using nslaw = siconos::vkernel::NS::nslaw;                                              \
    using polytope_shape = siconos::vkernel::NS::polytope_shape;                            \
    using segment_shape = siconos::vkernel::NS::segment_shape;                              \
    using polytopefsegment_r = siconos::vkernel::NS::polytopefsegment_r;                    \
    using fc2d = siconos::vkernel::NS::fc2d;                                                \
    using simulation = siconos::vkernel::NS::simulation;                                    \
  };                                                                                        \
  }

SICONOS_FREMOND_RIGID_CONFIG(newton_rigid, newton_impact_friction)
SICONOS_FREMOND_RIGID_CONFIG(fremond_rigid, fremond_impact_friction)

namespace {

constexpr double h = 1e-4, m = 1., g = 10.;

// a rigid body, its contact points given by the vertices of its shape
struct scene {
  std::vector<std::array<double, 2>> vertices;  // counterclockwise
  double inertia;
  std::array<double, 3> q0, v0;
  double mu, tmax;
  int npoints;  // the manifold points of the shape on the ground
};

// the bar of length 1: the lower end (0, -1/2) is its contact point, the
// upper vertices only close the polytope. The manifold is the incident edge,
// not sorted by depth: the lower end may be its point 0 or 1
scene stick() {
  const double gamma = M_PI / 4;
  return {{{0., -0.5}, {0.02, 0.5}, {-0.02, 0.5}},
          1. / 12.,
          {0.5 * std::sin(gamma), 0.5 * std::cos(gamma) + 0.01, gamma},
          {-0.5, 0.1, 0.1},
          0.01,
          0.2,
          2};
}

scene block(double tmax = 1.) {
  return {{{-0.5, -0.5}, {0.5, -0.5}, {0.5, 0.5}, {-0.5, 0.5}},
          1. / 6.,
          {0., 0.6, 0.},
          {0., -0.2, 1.},
          0.1,
          tmax,
          2};
}

// per step n (from 1): dE = E(n) - E(n - 1), and the contact work
// u(n - 1 + theta)^T p(n) of the active contacts, normal and tangential
struct history {
  std::vector<double> de, work_n, work_t;
  double e0, e1;
  std::array<double, 3> q, v;
  int active_steps;
};

template <typename Config>
history run(const scene& s, double theta, double e = 1.) {
  using namespace siconos::vkernel;
  using namespace storage;

  typename Config::data_t data;

  auto shape = storage::add<typename Config::polytope_shape>(data);
  shape.nvertices() = s.vertices.size();
  for (std::size_t i = 0; i < s.vertices.size(); ++i)
    shape.vertices()[i] = {s.vertices[i][0], s.vertices[i][1]};

  auto body = storage::add<typename Config::body>(data);
  storage::prop<"shape", storage::config_declared>(body) = shape;
  body.q() = {s.q0[0], s.q0[1], s.q0[2]};
  body.velocity() = {s.v0[0], s.v0[1], s.v0[2]};
  body.mass_matrix().diagonal() << m, m, s.inertia;
  body.fext() = {0., -m * g, 0.};

  auto nslaw = storage::add<typename Config::nslaw>(data);
  nslaw.e() = e;
  nslaw.mu() = s.mu;

  auto fc2d = storage::add<typename Config::fc2d>(data);
  fc2d.create();

  auto simulation = storage::add<typename Config::simulation>(data);
  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = 0.;
  simulation.time_discretization().h() = h;
  simulation.time_discretization().tmax() = s.tmax;

  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  so.set_dparam(SICONOS_DPARAM_TOL, 1e-14);
  so.set_iparam(SICONOS_IPARAM_MAX_ITER, 1000);
  osnspb.options() = so;

  auto ground = storage::add<typename Config::segment_shape>(data);
  ground.set_p1_p2(std::array<double, 6>{-10., 0., 0., 10., 0., 0.}, 0);
  ground.set_maxpoints(2);
  ground.initialize(0);

  // one interaction per manifold point of the body on the ground
  using relation_t = decltype(storage::add<typename Config::polytopefsegment_r>(data));
  using interaction_t = decltype(simulation.topology().link(body));
  std::vector<std::pair<interaction_t, relation_t>> contacts;
  for (int i = 0; i < s.npoints; ++i) {
    auto rel = storage::add<typename Config::polytopefsegment_r>(data);
    rel.segment() = ground;
    rel.point_index() = i;
    auto inter = simulation.topology().link(body);
    simulation.topology().bind(inter, rel);
    inter.nslaw() = nslaw;
    contacts.emplace_back(inter, rel);
  }

  simulation.initialize();

  auto energy = [&](const auto& q, const auto& v) {
    return 0.5 * m * (v(0) * v(0) + v(1) * v(1)) + 0.5 * s.inertia * v(2) * v(2) +
           m * g * q(1);
  };

  history r{};
  {
    auto step = simulation.current_step();
    r.e0 = energy(storage::attr<"q">(body, step), storage::attr<"velocity">(body, step));
  }
  while (simulation.has_next_event()) {
    auto step = simulation.current_step();
    // copies: the two time levels share their storage
    Eigen::Vector3d q_k = storage::attr<"q">(body, step);
    Eigen::Vector3d v_k = storage::attr<"velocity">(body, step);

    // H(q(k)) and the active set of the scheme, g + h/2 u_N <= 0
    std::vector<Eigen::MatrixXd> jacobians;
    std::vector<bool> active;
    for (auto& [inter, rel] : contacts) {
      Eigen::MatrixXd jac(2, 3);
      rel.compute_jachq(step, body, jac);
      active.push_back(rel.compute_h(step, body) + 0.5 * h * (jac * v_k)(0) <= 0.);
      jacobians.push_back(jac);
    }

    simulation.compute_one_step();

    step = simulation.current_step();
    Eigen::Vector3d q_1 = storage::attr<"q">(body, step);
    Eigen::Vector3d v_1 = storage::attr<"velocity">(body, step);
    Eigen::Vector3d v_theta = theta * v_1 + (1 - theta) * v_k;

    double wn = 0., wt = 0.;
    for (std::size_t i = 0; i < contacts.size(); ++i) {
      if (!active[i]) continue;
      auto lambda = storage::attr<"lambda">(contacts[i].first);
      Eigen::Vector2d u_theta = jacobians[i] * v_theta;
      wn += u_theta(0) * lambda(0);
      wt += u_theta(1) * lambda(1);
    }
    r.active_steps += std::ranges::count(active, true) > 0;
    r.de.push_back(energy(q_1, v_1) - energy(q_k, v_k));
    r.work_n.push_back(wn);
    r.work_t.push_back(wt);
  }
  auto step = simulation.current_step();
  auto q = storage::attr<"q">(body, step);
  auto v = storage::attr<"velocity">(body, step);
  r.e1 = energy(q, v);
  for (int i = 0; i < 3; ++i) {
    r.q[i] = q(i);
    r.v[i] = v(i);
  }
  return r;
}

void print(const char* name, const history& r) {
  auto imax = std::ranges::max_element(r.de) - r.de.begin();
  auto imin = std::ranges::min_element(r.de) - r.de.begin();
  std::print(
      "{}: E(T) - E(0) = {:.9e}, max dE = {:.9e} at step {}, min dE = {:.9e} at step {}, "
      "{} active steps\n",
      name, r.e1 - r.e0, r.de[imax], imax + 1, r.de[imin], imin + 1, r.active_steps);
}

// the largest |dE - contact work| and |normal work| over the run
double balance_defect(const history& r) {
  double d = 0.;
  for (std::size_t i = 0; i < r.de.size(); ++i)
    d = std::max(d, std::abs(r.de[i] - r.work_n[i] - r.work_t[i]));
  return d;
}

double max_normal_work(const history& r) {
  double d = 0.;
  for (auto w : r.work_n) d = std::max(d, std::abs(w));
  return d;
}

class FremondRigidTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(FremondRigidTest);
  CPPUNIT_TEST(testStickNewtonCreatesEnergy);
  CPPUNIT_TEST(testStickFremondDissipates);
  CPPUNIT_TEST(testBlockSecondImpact);
  CPPUNIT_TEST(testEnergyBalance);
  CPPUNIT_TEST(testBlockThetaOneIsNewton);
  CPPUNIT_TEST_SUITE_END();

  // the energy identities hold to round-off; the reference values (the same
  // scheme, the same activation, an exact local solver instead of NSGS) are
  // given to 7 digits, hence the 1e-9 and 1e-8 tolerances
  static constexpr double eps = 1e-12;

 public:
  // first impact, step 604: friction works +7.07e-4, E(0.2) - E(0) > 0
  void testStickNewtonCreatesEnergy() {
    auto r = run<siconos::vkernel::newton_rigid::types>(stick(), 0.5);
    print("stick, Newton", r);
    CPPUNIT_ASSERT_EQUAL(2, r.active_steps);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.070361e-4, r.de[603], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.070361e-4, *std::ranges::max_element(r.de), 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1.083148e-4, r.e1 - r.e0, 1e-9);
  }

  // the tangential work is never positive, E(0.2) - E(0) = -1.303e-3
  void testStickFremondDissipates() {
    auto r = run<siconos::vkernel::fremond_rigid::types>(stick(), 0.5);
    print("stick, Frémond", r);
    CPPUNIT_ASSERT_EQUAL(2, r.active_steps);
    CPPUNIT_ASSERT(*std::ranges::max_element(r.de) <= eps);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-6.768505e-4, r.de[603], 1e-9);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.303179e-3, r.e1 - r.e0, 1e-8);
  }

  // second impact, step 1610, right corner: Newton +0.0552, Frémond -0.0183;
  // Frémond never creates energy
  void testBlockSecondImpact() {
    auto n = run<siconos::vkernel::newton_rigid::types>(block(), 0.5);
    auto f = run<siconos::vkernel::fremond_rigid::types>(block(), 0.5);
    print("block, Newton", n);
    print("block, Frémond", f);
    // first impact, step 892: the slip does not reverse, the laws agree
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.121653e-2, n.de[891], 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.121653e-2, f.de[891], 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.518023e-2, n.de[1609], 1e-8);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.825896e-2, f.de[1609], 1e-8);
    CPPUNIT_ASSERT(*std::ranges::max_element(n.de) > 5e-2);
    CPPUNIT_ASSERT(*std::ranges::max_element(f.de) <= eps);
  }

  // E(k + 1) - E(k) = u(k + theta)^T p(k + 1): lambda is the impulse p(k + 1)
  // with both laws, and the normal work is zero for e = 1, theta = 1/2
  void testEnergyBalance() {
    for (auto r : {run<siconos::vkernel::newton_rigid::types>(block(0.5), 0.5),
                   run<siconos::vkernel::fremond_rigid::types>(block(0.5), 0.5),
                   run<siconos::vkernel::newton_rigid::types>(stick(), 0.5),
                   run<siconos::vkernel::fremond_rigid::types>(stick(), 0.5)}) {
      std::print("balance defect {:.3e}, max |normal work| {:.3e}\n", balance_defect(r),
                 max_normal_work(r));
      CPPUNIT_ASSERT(balance_defect(r) <= eps);
      CPPUNIT_ASSERT(max_normal_work(r) <= eps);
    }
  }

  // two coupled contacts: theta = 1 Frémond is Newton, bit for bit
  void testBlockThetaOneIsNewton() {
    auto n = run<siconos::vkernel::newton_rigid::types>(block(0.5), 1., 0.5);
    auto f = run<siconos::vkernel::fremond_rigid::types>(block(0.5), 1., 0.5);
    print("block theta 1, Newton", n);
    print("block theta 1, Frémond", f);
    CPPUNIT_ASSERT(n.active_steps > 0);
    for (int i = 0; i < 3; ++i) {
      CPPUNIT_ASSERT_EQUAL(n.q[i], f.q[i]);
      CPPUNIT_ASSERT_EQUAL(n.v[i], f.v[i]);
    }
    CPPUNIT_ASSERT(n.de == f.de);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(FremondRigidTest);

}  // namespace
