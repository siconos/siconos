// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cmath>
#include <print>

#include "siconos/vkernel.hpp"

// A disk without gravity falls vertically (vy = -1) on a fixed segment, e = 0,
// mu = 1.  The normal shock law is a shift of the free output derivative q_N:
//
//   Newton, e = 0          q_N += e ydot_N = 0
//   Jean 1999 Eq 25        q_N += max(0, y) / h          (g_pos: no push-out)
//   Anitescu-Tasora        q_N += max(y / h, -v_r)       (signed gap, pushes out)
//
// The scheme is chosen at run time (selectable_scheme.hpp):
//   local_problem (0), local_problem_quasi_inelastic (6),
//   cone_complementarity_inelastic (7).
//
// Closed forms used below (m = 1, no gravity, |v_n^-| = 1):
//   - the disk stops: v_n^+ = 0 (Newton, quasi inelastic) and stays (no gravity);
//   - the normal impulse balances the momentum, sum(lambda_n) = m |v_n^-| = 1,
//     and it is spread over more than one step (the gap is closed first);
//   - Jean closes the gap and never ejects (v_n is never > 0);
//   - Jean with g_pos keeps a small overlap for theta < 1 (no push-out);
//   - the AT stabilisation pushes the overlap out (the gap is >= 0).

#define SICONOS_QI_DISK_CONFIG(NS)                                                   \
  namespace siconos::vkernel::NS {                                                   \
  struct disk : model::lagrangian_ds {};                                             \
  struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};                 \
  struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};                        \
  struct nslaw : model::newton_impact_friction {};                                   \
  using disk_shape = collision::shape::disk;                                         \
  using diskdisk_r = collision::diskdisk_r;                                          \
  using diskfsegment_r = collision::diskfsegment_r;                                  \
  struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r> {};     \
  struct topo : simul::topology<disk, interaction> {};                               \
  struct osi : simul::one_step_integrator<topo>::selectable {};                      \
  struct td : simul::time_discretization<> {};                                       \
  struct simulation : simul::time_stepping<td, osi, osnspb> {};                      \
  template <typename T>                                                              \
  struct env : standard_environment<T> {                                             \
    using params = config::map<config::iparam<"dof", 3>>;                            \
  };                                                                                 \
  struct data_t : storage::make_runtime_database<                                    \
                      env, simulation, disk_shape,                                   \
                      storage::with_properties<                                      \
                          storage::attached<disk, storage::pattern::symbol<"shape">, \
                                            storage::some::item_ref<disk_shape>>,    \
                          storage::time_invariant<storage::attr_t<disk, "fext">>,    \
                          storage::diagonal<storage::attr_t<disk, "mass_matrix">>,   \
                          storage::assembled_diagonal<                               \
                              storage::attr_t<osi, "mass_matrix_assembled">>>> {};   \
  struct types {                                                                     \
    using data_t = siconos::vkernel::NS::data_t;                                     \
    using disk = siconos::vkernel::NS::disk;                                         \
    using nslaw = siconos::vkernel::NS::nslaw;                                       \
    using fc2d = siconos::vkernel::NS::fc2d;                                         \
    using simulation = siconos::vkernel::NS::simulation;                             \
    using diskfsegment_r = siconos::vkernel::NS::diskfsegment_r;                     \
  };                                                                                 \
  }

SICONOS_QI_DISK_CONFIG(qi_config)

namespace {

struct outcome {
  double vx, vy, omega;  // the final velocity of the disk
  double gap;            // the final gap, q_y - radius
  double max_vy;         // the largest normal velocity (ejection if > 0)
  double impulse;        // sum of the normal multiplier over the run
  int impact_steps;      // steps whose normal multiplier is not zero
};

// scheme: one of simul::selectable_scheme's values
template <typename Config>
outcome drop(double theta, std::size_t scheme) {
  using namespace siconos::vkernel;
  using namespace storage;

  typename Config::data_t data;

  double h = 1e-3, radius = 0.1, m = 1.;
  double gap0 = 0.0105;  // 0.5 h |v| above contact: a genuine two-step shock
  double inertia = 0.5 * m * radius * radius;

  auto d1 = storage::add<typename Config::disk>(data);
  d1.q() = {0, radius + gap0, 0};
  d1.velocity() = {0., -1., 0.};
  d1.mass_matrix().diagonal() << m, m, inertia;
  storage::make_handle(data, prop<"shape">(d1)).radius() = radius;
  d1.fext() = {0., 0., 0.};

  auto nslaw = storage::add<typename Config::nslaw>(data);
  nslaw.e() = 0.;
  nslaw.mu() = 1.;

  auto fc2d = storage::add<typename Config::fc2d>(data);
  fc2d.create();

  auto simulation = storage::add<typename Config::simulation>(data);
  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().scheme() = scheme;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = 0.;
  simulation.time_discretization().h() = h;
  simulation.time_discretization().tmax() = 0.05;

  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  so.set_dparam(SICONOS_DPARAM_TOL, 1e-14);
  so.set_iparam(SICONOS_IPARAM_MAX_ITER, 1000);
  osnspb.options() = so;

  auto ground_r = storage::add<typename Config::diskfsegment_r>(data);
  auto segment = ground_r.segment();
  segment.x1() = -10.;
  segment.y1() = 0;
  segment.x2() = 10;
  segment.y2() = 0;
  segment.initialize();
  segment.maxpoints() = 1000;

  auto interaction = simulation.topology().link(d1);
  simulation.topology().bind(interaction, ground_r);
  interaction.nslaw() = nslaw;

  simulation.initialize();

  outcome r{};
  double vy_prev = d1.velocity()(1);
  while (simulation.has_next_event()) {
    simulation.compute_one_step();
    const double vy_now = storage::attr<"velocity">(d1, simulation.current_step())(1);
    // the impact steps: those whose velocity changes (a resting contact at gap 0
    // is a degenerate complementarity whose lambda need not balance anything)
    if (std::abs(vy_now - vy_prev) > 1e-12) {
      ++r.impact_steps;
      r.impulse += interaction.lambda()(0);
    }
    r.max_vy = std::max(r.max_vy, vy_now);
    vy_prev = vy_now;
  }
  auto v = storage::attr<"velocity">(d1, simulation.current_step());
  auto q = storage::attr<"q">(d1, simulation.current_step());
  r.vx = v(0);
  r.vy = v(1);
  r.omega = v(2);
  r.gap = q(1) - radius;
  std::print(
      "theta {} scheme {}: v = ({:.6e}, {:.6e}, {:.6e}), gap {:.6e}, "
      "max vy {:.6e}, impulse {:.6e} over {} steps\n",
      theta, scheme, r.vx, r.vy, r.omega, r.gap, r.max_vy, r.impulse, r.impact_steps);
  return r;
}

class QuasiInelasticDiskTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(QuasiInelasticDiskTest);
  CPPUNIT_TEST(testQuasiInelasticThetaOne);
  CPPUNIT_TEST(testQuasiInelasticThetaHalf);
  CPPUNIT_TEST(testQuasiInelasticTwoStepImpulse);
  CPPUNIT_TEST(testAnitescuTasoraPushesOut);
  CPPUNIT_TEST(testNewtonBaseline);
  CPPUNIT_TEST_SUITE_END();

  static constexpr double tol = 1e-8;

  // theta = 1: the gap closes exactly, v_n = 0, no ejection, at rest
  void testQuasiInelasticThetaOne() {
    namespace s = siconos::vkernel::simul::selectable_scheme;
    auto r = drop<siconos::vkernel::qi_config::types>(1., s::local_problem_quasi_inelastic);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., r.vy, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., r.gap, 1e-6);
    CPPUNIT_ASSERT(r.max_vy <= tol);  // never ejects
    CPPUNIT_ASSERT(r.impact_steps >= 2);
  }

  // theta = 1/2: the gap closes, no ejection, at rest, and g_pos keeps a small
  // permanent overlap (no push-out) - within the documented range
  void testQuasiInelasticThetaHalf() {
    namespace s = siconos::vkernel::simul::selectable_scheme;
    auto r = drop<siconos::vkernel::qi_config::types>(0.5, s::local_problem_quasi_inelastic);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., r.vy, 1e-3);
    CPPUNIT_ASSERT(r.max_vy <= 1e-4);               // never ejects
    CPPUNIT_ASSERT(r.gap > -1e-2 && r.gap < 1e-6);  // overlap kept, not pushed out
    CPPUNIT_ASSERT(r.impact_steps >= 2);
  }

  // the two-step property: sum(lambda_n) = m |v_n^-| = 1, over more than one step
  void testQuasiInelasticTwoStepImpulse() {
    namespace s = siconos::vkernel::simul::selectable_scheme;
    for (double theta : {0.5, 1.}) {
      auto r =
          drop<siconos::vkernel::qi_config::types>(theta, s::local_problem_quasi_inelastic);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1., r.impulse, 1e-6);
    }
  }

  // the AT stabilisation uses the signed gap: an overlap is pushed out, the
  // final gap is non negative (where Jean keeps it)
  void testAnitescuTasoraPushesOut() {
    namespace s = siconos::vkernel::simul::selectable_scheme;
    auto r = drop<siconos::vkernel::qi_config::types>(0.5, s::cone_complementarity_inelastic);
    CPPUNIT_ASSERT(r.gap >= -1e-6);   // pushed out, no residual overlap
    CPPUNIT_ASSERT(r.max_vy > 1e-3);  // the recovery push-out
  }

  // the Newton baseline, e = 0: stops, no restitution, a small overlap
  void testNewtonBaseline() {
    namespace s = siconos::vkernel::simul::selectable_scheme;
    auto r = drop<siconos::vkernel::qi_config::types>(0.5, s::local_problem);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., r.vy, 1e-6);
    CPPUNIT_ASSERT(r.gap <= 1e-6);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(QuasiInelasticDiskTest);

}  // namespace
