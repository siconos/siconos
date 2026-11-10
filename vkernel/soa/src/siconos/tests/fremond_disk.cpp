#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cmath>
#include <print>

#include "siconos/siconos.hpp"

// A disk without gravity hits a fixed segment obliquely: vx = 1, vy = -1,
// no spin, I = 2/5 m R^2, so that W_NN = 1/m, W_TT = 7/(2m), W_NT = 0. With
// enough friction the contact sticks, and the velocities after the impact
// are known in closed form (Acary & Collins-Craft, JTCAM 2025):
//
//   Newton:  u_T+ = 0                              vx+ = 5/7
//   Fremond: theta u_T+ + (1 - theta) u_T- = 0     vx+ = 1 - 2/(7 theta)
//
// i.e. vx+ = 3/7 with theta = 1/2, where the kinetic energy is conserved
// for e = 1, and the Newton value with theta = 1. vy+ = e in both.

#define SICONOS_FREMOND_DISK_CONFIG(NS, LAW)                                         \
  namespace siconos::NS {                                                            \
  struct disk : model::lagrangian_ds {};                                             \
  struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};                 \
  struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};                        \
  struct nslaw : model::LAW {};                                                      \
  using disk_shape = collision::shape::disk;                                         \
  using diskdisk_r = collision::diskdisk_r;                                          \
  using diskfsegment_r = collision::diskfsegment_r;                                  \
  struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r> {};     \
  struct topo : simul::topology<disk, interaction> {};                               \
  struct osi : simul::one_step_integrator<topo>::moreau_jean {};                     \
  struct td : simul::time_discretization<> {};                                       \
  struct simulation : simul::time_stepping<td, osi, osnspb> {};                      \
  template <typename T>                                                              \
  struct env : standard_environment<T> {                                             \
    using params = config::map<config::iparam<"dof", 3>>;                            \
  };                                                                                 \
  struct data_t : storage::make_compile_time_database<                               \
                      env, simulation, disk_shape,                                   \
                      storage::with_properties<                                      \
                          storage::attached<disk, storage::pattern::symbol<"shape">, \
                                            storage::some::item_ref<disk_shape>>,    \
                          storage::time_invariant<storage::attr_t<disk, "fext">>,    \
                          storage::diagonal<storage::attr_t<disk, "mass_matrix">>,   \
                          storage::assembled_diagonal<                               \
                              storage::attr_t<osi, "mass_matrix_assembled">>>> {};   \
  /* the types of the config, as one template argument */                            \
  struct types {                                                                     \
    using data_t = siconos::NS::data_t;                                              \
    using disk = siconos::NS::disk;                                                  \
    using nslaw = siconos::NS::nslaw;                                                \
    using fc2d = siconos::NS::fc2d;                                                  \
    using simulation = siconos::NS::simulation;                                      \
    using diskfsegment_r = siconos::NS::diskfsegment_r;                              \
  };                                                                                 \
  }

SICONOS_FREMOND_DISK_CONFIG(newton_config, newton_impact_friction)
SICONOS_FREMOND_DISK_CONFIG(fremond_config, fremond_impact_friction)

namespace {

struct outcome {
  double vx, vy, omega, kinetic_energy_before, kinetic_energy_after;
  int impacts;
};

template <typename Config>
outcome oblique_impact(double theta, double e, double mu) {
  using namespace siconos;
  using namespace storage;

  typename Config::data_t data;

  double h = 1e-3, radius = 0.1, m = 1.;
  double inertia = 2. / 5. * m * radius * radius;

  auto d1 = storage::add<typename Config::disk>(data);
  d1.q() = {0, radius + 0.0105, 0};
  d1.velocity() = {1., -1., 0};
  d1.mass_matrix().diagonal() << m, m, inertia;
  storage::make_handle(data, prop<"shape">(d1)).radius() = radius;
  d1.fext() = {0., 0., 0.};

  auto nslaw = storage::add<typename Config::nslaw>(data);
  nslaw.e() = e;
  nslaw.mu() = mu;

  auto fc2d = storage::add<typename Config::fc2d>(data);
  fc2d.create();

  auto simulation = storage::add<typename Config::simulation>(data);
  simulation.one_step_integrator().theta() = theta;
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

  auto kinetic_energy = [&](const auto& v) {
    return 0.5 * m * (v(0) * v(0) + v(1) * v(1)) + 0.5 * inertia * v(2) * v(2);
  };

  outcome r{};
  r.kinetic_energy_before = kinetic_energy(d1.velocity());
  while (simulation.has_next_event()) {
    // a copy: the velocities at two time levels share their storage
    double vy_k = storage::attr<"velocity">(d1, simulation.current_step())(1);
    simulation.compute_one_step();
    if (storage::attr<"velocity">(d1, simulation.current_step())(1) != vy_k) ++r.impacts;
  }
  auto v = storage::attr<"velocity">(d1, simulation.current_step());
  r.vx = v(0);
  r.vy = v(1);
  r.omega = v(2);
  r.kinetic_energy_after = kinetic_energy(v);
  std::print("theta {} e {} mu {}: v = ({:.15e}, {:.15e}, {:.15e}), KE {:.15e} -> {:.15e}\n",
             theta, e, mu, r.vx, r.vy, r.omega, r.kinetic_energy_before,
             r.kinetic_energy_after);
  return r;
}

class FremondDiskTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(FremondDiskTest);
  CPPUNIT_TEST(testNewtonSticks);
  CPPUNIT_TEST(testFremondSticksOnAverage);
  CPPUNIT_TEST(testFremondThetaOneIsNewton);
  CPPUNIT_TEST(testFremondSlidingDissipates);
  CPPUNIT_TEST_SUITE_END();

  static constexpr double tol = 1e-10;

 public:
  void testNewtonSticks() {
    auto r = oblique_impact<siconos::newton_config::types>(0.5, 1., 1.);
    CPPUNIT_ASSERT_EQUAL(1, r.impacts);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5. / 7., r.vx, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1., r.vy, tol);
  }

  // u_T reverses: the kinetic energy is conserved (e = 1, theta = 1/2)
  void testFremondSticksOnAverage() {
    auto r = oblique_impact<siconos::fremond_config::types>(0.5, 1., 1.);
    CPPUNIT_ASSERT_EQUAL(1, r.impacts);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3. / 7., r.vx, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1., r.vy, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(r.kinetic_energy_before, r.kinetic_energy_after, tol);
  }

  void testFremondThetaOneIsNewton() {
    auto n = oblique_impact<siconos::newton_config::types>(1., 0.5, 1.);
    auto f = oblique_impact<siconos::fremond_config::types>(1., 0.5, 1.);
    CPPUNIT_ASSERT_EQUAL(n.vx, f.vx);
    CPPUNIT_ASSERT_EQUAL(n.vy, f.vy);
    CPPUNIT_ASSERT_EQUAL(n.omega, f.omega);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5. / 7., f.vx, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5, f.vy, tol);
  }

  // small mu: the contact slides, |p_T| = mu p_N = 2 mu, vx+ = 1 - 2 mu,
  // and friction dissipates
  void testFremondSlidingDissipates() {
    auto r = oblique_impact<siconos::fremond_config::types>(0.5, 1., 0.1);
    CPPUNIT_ASSERT_EQUAL(1, r.impacts);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8, r.vx, tol);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(1., r.vy, tol);
    CPPUNIT_ASSERT(r.kinetic_energy_after < r.kinetic_energy_before - 1e-3);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(FremondDiskTest);

}  // namespace
