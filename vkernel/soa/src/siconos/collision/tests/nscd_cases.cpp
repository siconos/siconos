/* Analytic test cases of the NSCD scheme (Moreau-Jean, velocity level,
   Coulomb friction), the reference numbers for the comparisons with the
   Anitescu-Tasora and explicit DEM schemes (design note "Plan: where NSCD
   wins, reference schemes for comparison").

   Rigid bodies on a ground segment y = 0, the space filter, the cxd
   narrowphase stage, Moreau-Jean theta = 0.5, the 2D friction problem
   solved by NSGS at a tolerance of 1e-12, so that the errors below are the
   scheme's. Not Lemke: on the step where the sliding square stops,
   fc2d_lemke returns a wrong solution (a penetrating contact, a friction
   of the wrong sign) with info 0. An incline of angle alpha is the ground with the
   gravity tilted by alpha: the same equations, in the frame of the incline.

   1. sliding square: launched at v0 on the ground with friction mu; the
      exact stopping distance v0^2 / (2 mu g), no lift-off (y and vy stay at
      their contact values: the convexified cone of Anitescu-Tasora would
      lift the block by ~ h mu |u_T| per step).
   2. square on an incline below the friction angle (tan alpha < mu): it
      sticks, no creep.
   3. square on an incline above the friction angle: the exact acceleration
      g (sin alpha - mu cos alpha).
   4. a column of N disks at rest in contact: no interpenetration at any
      contact, whatever N (a penalty scheme has a bottom overlap ~ N m g /
      k_n).

   Each case prints its measured numbers. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>
#include <cmath>
#include <numbers>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
#include "siconos/collision/polytopedisk_r.hpp"
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

struct body : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using disk_shape = collision::shape::disk;
using polytope_shape = collision::shape::polytope;
using segment_shape = collision::shape::segment;
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using polytopedisk_r = collision::polytopedisk_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, polytopepolytope_r,
                                        polytopefsegment_r, polytopedisk_r> {};
struct topo : simul::topology<body, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<body, collision::empty_shape>;
using pointl = collision::point<storage::pattern::empty_item, collision::shape::segment>;
struct neighborhood : collision::neighborhood<pointd, pointl> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, interaction_manager, body, disk_shape, polytope_shape,
          segment_shape, diskdisk_r, diskfsegment_r, nslaw, interaction, topo,
          polytopepolytope_r, polytopefsegment_r, cxd_narrowphase, osi, td, osnspb,
          neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopefsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopedisk_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>,
              storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                         some::item_ref<polytope_shape>>>,
              storage::time_invariant<storage::pattern::attr_t<body, "fext">>,
              storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::pattern::attr_t<osi, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {

constexpr double m = 1.;
constexpr double g = 9.81;
constexpr double r = 0.5;  // disk radius, square half side
constexpr double deg = std::numbers::pi / 180.;

struct scene {
  double mu = 0.;
  double h = 0.005;
  double tmax = 1.;
  double alpha = 0.;  // incline angle: the gravity tilted towards +x
  // a contact resting at gap 0 has y and ydot of the order of the rounding
  // errors: with a threshold of 0 its activation is decided by their sign,
  // and a dropped row lets the body fall by g h^2 on the other contact
  double threshold = 1e-10;
};

struct world {
  siconos::config::data_t data;
  decltype(siconos::storage::add<siconos::config::simulation>(data)) simulation;
  decltype(siconos::storage::add<siconos::config::space_filter>(data)) spacef;
  decltype(siconos::storage::add<siconos::config::neighborhood>(data)) ngbh;
  decltype(siconos::storage::add<siconos::config::cxd_narrowphase>(data)) np;
  decltype(siconos::storage::add<siconos::config::disk_shape>(data)) dsh;
  decltype(siconos::storage::add<siconos::config::polytope_shape>(data)) psh;
  scene sc;

  explicit world(scene sc_)
      : simulation(siconos::storage::add<siconos::config::simulation>(data)),
        spacef(siconos::storage::add<siconos::config::space_filter>(data)),
        ngbh(siconos::storage::add<siconos::config::neighborhood>(data)),
        np(siconos::storage::add<siconos::config::cxd_narrowphase>(data)),
        dsh(siconos::storage::add<siconos::config::disk_shape>(data)),
        psh(siconos::storage::add<siconos::config::polytope_shape>(data)),
        sc(sc_) {
    using namespace siconos;
    dsh.radius() = r;
    const std::array<std::array<double, 2>, 4> square = {{{-r, -r}, {r, -r}, {r, r}, {-r, r}}};
    for (int i = 0; i < 4; ++i) psh.vertices()[i] = {square[i][0], square[i][1]};
    psh.nvertices() = 4;

    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = 0.;
    nslaw.mu() = sc.mu;

    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

    simulation.one_step_integrator().theta() = 0.5;
    simulation.one_step_integrator().constraint_activation_threshold() = sc.threshold;
    simulation.time_discretization().t0() = 0.;
    simulation.time_discretization().h() = sc.h;
    simulation.time_discretization().tmax() = sc.tmax;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;
    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_NSGS);
    so.set_iparam(SICONOS_IPARAM_MAX_ITER, 100000);
    so.set_dparam(SICONOS_DPARAM_TOL, 1e-12);
    osnspb.options() = so;

    // ground y = 0, sampled every 0.25 so that a body centre in contact
    // (y = 0.5) always has a sample point inside the radius
    auto segment = storage::add<config::segment_shape>(data);
    segment.x1() = -20.;
    segment.y1() = 0.;
    segment.x2() = 20.;
    segment.y2() = 0.;
    segment.initialize();
    segment.maxpoints() = 161;

    ngbh.create(1.5);
    auto dsds = storage::add<config::diskdisk_r>(data);  // generic ds-ds binder
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = dsds;
    spacef.insert_nonsmooth_law(nslaw, 0, 0);
    np.topology() = spacef.topology();
    np.dsds_r() = dsds;
  }

  auto add_body(auto shape, double x, double y, double vx, double inertia) {
    using namespace siconos;
    auto b = storage::add<config::body>(data);
    b.q() = {x, y, 0};
    b.velocity() = {vx, 0, 0};
    b.mass_matrix().diagonal() << m, m, inertia;
    b.fext() = {m * g * std::sin(sc.alpha), -m * g * std::cos(sc.alpha), 0.};
    collision::set_shape(b, shape);
    return b;
  }
  auto add_disk(double x, double y) { return add_body(dsh, x, y, 0., 0.5 * m * r * r); }
  auto add_square(double x, double y, double vx = 0.) {
    return add_body(psh, x, y, vx, m * (4 * r * r + 4 * r * r) / 12.);
  }

  void finalize() {
    spacef.make_points();
    ngbh.add_point_sets(0);
    simulation.initialize();
  }

  void step() {
    const auto s = simulation.current_step();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    np.ensure(s);
    simulation.compute_one_step();
  }
};

/* 1. the sliding square: the friction impulse is mu m g h at every sliding
   step, so the velocity decreases exactly by mu g h per step and the
   trapezoidal position update (theta = 0.5) is exact on it; only the last,
   partial step (the block sticks inside it) differs from the analytic motion,
   by at most v h with v < mu g h. */
void sliding_square(double h) {
  using namespace siconos;
  const double v0 = 2.;
  const double mu = 0.3;
  const double t_stop = v0 / (mu * g);
  const double d_exact = v0 * v0 / (2 * mu * g);
  world w({.mu = mu, .h = h, .tmax = t_stop + 0.3});
  auto p = w.add_square(0., r, v0);
  w.finalize();

  double max_dy = 0., max_vy = 0., max_vx_err = 0., max_theta = 0.;
  double t_rest = -1.;
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    const auto& q = storage::attr<"q">(p, s);
    const auto& v = storage::attr<"velocity">(p, s);
    const double t = s * h;
    max_dy = std::max(max_dy, std::abs(q(1) - r));
    max_vy = std::max(max_vy, std::abs(v(1)));
    max_theta = std::max(max_theta, std::abs(q(2)));
    max_vx_err = std::max(max_vx_err, std::abs(v(0) - std::max(0., v0 - mu * g * t)));
    if (t_rest < 0 && std::abs(v(0)) <= 1e-12) t_rest = t;
  }
  const auto s = w.simulation.current_step();
  const double d = storage::attr<"q">(p, s)(0);
  std::print(stderr,
             "sliding square h = {:.0e}: distance {:.9f} (exact {:.9f}, error {:.2e}), "
             "stops at t = {:.4f} (exact {:.4f}), max |y - r| {:.2e}, max |vy| {:.2e}, "
             "max |vx - vx_exact| {:.2e}, max |theta| {:.2e}\n",
             h, d, d_exact, d - d_exact, t_rest, t_stop, max_dy, max_vy, max_vx_err,
             max_theta);
  const double bound = mu * g * h * h + 1e-12;
  CPPUNIT_ASSERT_MESSAGE("the stopping distance is v0^2 / (2 mu g)",
                         std::abs(d - d_exact) <= bound);
  CPPUNIT_ASSERT_MESSAGE("the velocity decreases by mu g h per step",
                         max_vx_err <= mu * g * h + 1e-12);
  CPPUNIT_ASSERT_MESSAGE("the block stops at v0 / (mu g)",
                         t_rest >= t_stop - 1e-12 && t_rest <= t_stop + h + 1e-12);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: y stays at r", max_dy <= 1e-12);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: vy stays 0", max_vy <= 1e-12);
  CPPUNIT_ASSERT_MESSAGE("no tipping (mu < 1): theta stays 0", max_theta <= 1e-12);
}

/* 2. and 3. a square at rest on an incline of angle alpha: below the
   friction angle it sticks, above it slides with the acceleration
   a = g (sin alpha - mu cos alpha), exactly at every step (the same
   argument as for the sliding square) */
void incline(double alpha, double mu, double h) {
  using namespace siconos;
  const double tmax = 1.;
  const bool sticks = std::tan(alpha) < mu;
  const double a = sticks ? 0. : g * (std::sin(alpha) - mu * std::cos(alpha));
  world w({.mu = mu, .h = h, .tmax = tmax, .alpha = alpha});
  auto p = w.add_square(0., r);
  w.finalize();

  double max_dy = 0., max_vy = 0., max_vx_err = 0., max_x_err = 0., max_theta = 0.;
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    const auto& q = storage::attr<"q">(p, s);
    const auto& v = storage::attr<"velocity">(p, s);
    const double t = s * h;
    max_dy = std::max(max_dy, std::abs(q(1) - r));
    max_vy = std::max(max_vy, std::abs(v(1)));
    max_theta = std::max(max_theta, std::abs(q(2)));
    max_vx_err = std::max(max_vx_err, std::abs(v(0) - a * t));
    max_x_err = std::max(max_x_err, std::abs(q(0) - 0.5 * a * t * t));
  }
  const auto s = w.simulation.current_step();
  std::print(stderr,
             "incline {:.0f} deg, mu {} ({}), h = {:.0e}: x(T) = {:.12f} (exact {:.12f}), "
             "max |x - x_exact| {:.2e}, max |vx - vx_exact| {:.2e}, max |y - r| {:.2e}, "
             "max |vy| {:.2e}, max |theta| {:.2e}\n",
             alpha / deg, mu, sticks ? "sticks" : "slides", h, storage::attr<"q">(p, s)(0),
             0.5 * a * tmax * tmax, max_x_err, max_vx_err, max_dy, max_vy, max_theta);
  // on the trapezoidal rule x is exact for a linear v: rounding only
  CPPUNIT_ASSERT_MESSAGE("the velocity is a t", max_vx_err <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("the position is a t^2 / 2", max_x_err <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: y stays at r", max_dy <= 1e-12);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: vy stays 0", max_vy <= 1e-12);
  CPPUNIT_ASSERT_MESSAGE("no tipping (mu < 1): theta stays 0", max_theta <= 1e-12);
}

/* 4. a column of n disks, in contact and at rest: the free velocity -g h of
   every disk is cancelled at the velocity level, so no contact ever closes
   further, at any depth of the column */
void column(std::size_t n, double h) {
  using namespace siconos;
  world w({.mu = 0.3, .h = h, .tmax = 1.});
  std::vector<decltype(w.add_disk(0., 0.))> disks;
  for (std::size_t i = 0; i < n; ++i) disks.push_back(w.add_disk(0., r + 2 * r * i));
  w.finalize();

  double max_pen = 0., max_v = 0.;
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    double below = 0.;  // the top of the ground, then of the disk below
    for (auto& d : disks) {
      const double y = storage::attr<"q">(d, s)(1);
      max_pen = std::max(max_pen, (below + r) - y);
      max_v = std::max(max_v, storage::attr<"velocity">(d, s).norm());
      below = y + r;
    }
  }
  std::print(stderr,
             "column of {} disks, h = {:.0e}: max interpenetration {:.2e}, max |v| {:.2e}\n",
             n, h, max_pen, max_v);
  // the rounding errors and the solver tolerance (1e-12) only
  CPPUNIT_ASSERT_MESSAGE("no interpenetration in the column", max_pen <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("the column stays at rest", max_v <= 1e-10);
}

class NscdCasesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(NscdCasesTest);
  CPPUNIT_TEST(testSlidingSquare);
  CPPUNIT_TEST(testInclineSticks);
  CPPUNIT_TEST(testInclineSlides);
  CPPUNIT_TEST(testColumn);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testSlidingSquare() {
    for (double h : {1e-2, 5e-3, 1e-3}) sliding_square(h);
  }
  void testInclineSticks() {
    for (double h : {1e-2, 5e-3, 1e-3}) incline(20 * deg, 0.5, h);
  }
  void testInclineSlides() {
    for (double h : {1e-2, 5e-3, 1e-3}) incline(30 * deg, 0.3, h);
  }
  void testColumn() {
    for (std::size_t n : {1, 5, 20, 50}) column(n, 5e-3);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(NscdCasesTest);

}  // namespace
