/* Analytic test cases of the NSCD scheme (Moreau-Jean, velocity level,
   Coulomb friction), the numbers for the comparisons with the
   Anitescu-Tasora and explicit DEM schemes (design note, "Plan: where NSCD
   wins").

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

   The integrator is the selectable one: the same cases checks the
   Anitescu-Tasora scheme (simul/anitescu_tasora.hpp, its sweeps at the same
   tolerance) where it must agree with Moreau-Jean, mu = 0 (a column of
   disks dropped with gaps, step by step; the frictionless incline), and
   where it must not: the sliding square lifts off, u_N = mu |u_T|. And the
   explicit DEM scheme (simul/dem.hpp, penalty springs, theta = 1: the
   symplectic Euler scheme), on its analytic numbers: the overlaps of a
   column at rest (N - i) m g / k_n, the stability limit h < 2 sqrt(m /
   k_n), the stopping distance, the elastic offset m g sin(alpha) / (2 k_t)
   of a square on an incline.

   Each case prints its measured numbers. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <print>
#include <vector>

#include "siconos/collision/tests/nscd_world.hpp"

namespace {

using namespace nscd;

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
void incline(double alpha, double mu, double h, std::size_t scheme = schemes::local_problem) {
  using namespace siconos;
  const double tmax = 1.;
  const bool sticks = std::tan(alpha) < mu;
  const double a = sticks ? 0. : g * (std::sin(alpha) - mu * std::cos(alpha));
  world w({.mu = mu, .h = h, .tmax = tmax, .alpha = alpha, .scheme = scheme});
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
             "{}incline {:.0f} deg, mu {} ({}), h = {:.0e}: x(T) = {:.12f} (exact {:.12f}), "
             "max |x - x_exact| {:.2e}, max |vx - vx_exact| {:.2e}, max |y - r| {:.2e}, "
             "max |vy| {:.2e}, max |theta| {:.2e}\n",
             scheme == schemes::local_problem ? "" : "AT: ", alpha / deg, mu,
             sticks ? "sticks" : "slides", h, storage::attr<"q">(p, s)(0),
             0.5 * a * tmax * tmax, max_x_err, max_vx_err, max_dy, max_vy, max_theta);
  // on the trapezoidal rule x is exact for a linear v: the rounding errors
  // and the solver tolerance (1e-12) only
  CPPUNIT_ASSERT_MESSAGE("the velocity is a t", max_vx_err <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("the position is a t^2 / 2", max_x_err <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: y stays at r", max_dy <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("no lift-off: vy stays 0", max_vy <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("no tipping (mu < 1): theta stays 0", max_theta <= 1e-10);
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

// the Anitescu-Tasora solvers, in the measured lines
const char* name(std::size_t scheme) {
  return scheme == schemes::cone_complementarity_apgd ? "AT (APGD)" : "AT (PGS)";
}

/* Anitescu-Tasora, mu = 0: the same law as Moreau-Jean. A column of disks
   dropped with gaps of 0.1 (impacts at different times, then at rest), the
   two schemes step by step: their difference is the solvers' tolerance */
void frictionless_column(std::size_t n, std::size_t scheme) {
  using namespace siconos;
  auto run = [n](std::size_t scheme) {
    world w({.mu = 0., .h = 5e-3, .tmax = 1., .scheme = scheme});
    std::vector<decltype(w.add_disk(0., 0.))> disks;
    for (std::size_t i = 0; i < n; ++i)
      disks.push_back(w.add_disk(0., r + 0.1 + (2 * r + 0.1) * i));
    w.finalize();
    std::vector<double> trajectory;
    while (w.simulation.has_next_event()) {
      w.step();
      const auto s = w.simulation.current_step();
      for (auto& d : disks) {
        trajectory.push_back(storage::attr<"q">(d, s)(1));
        trajectory.push_back(storage::attr<"velocity">(d, s)(1));
      }
    }
    return trajectory;
  };
  const auto mj = run(schemes::local_problem);
  const auto at = run(scheme);
  double diff = 0., lowest = 1.;
  for (std::size_t k = 0; k < mj.size(); ++k) diff = std::max(diff, std::abs(mj[k] - at[k]));
  for (std::size_t k = 0; k < mj.size(); k += 2 * n) lowest = std::min(lowest, mj[k]);
  std::print(stderr,
             "{}: frictionless column of {} disks dropped with gaps: max |AT - MJ| on y, vy "
             "{:.2e} over {} steps (lowest disk down to y = {:.6f})\n",
             name(scheme), n, diff, mj.size() / (2 * n), lowest);
  CPPUNIT_ASSERT_MESSAGE("the steps ran", mj.size() == at.size() && !mj.empty());
  CPPUNIT_ASSERT_MESSAGE("the disks landed", lowest < r + 1e-3);
  CPPUNIT_ASSERT_MESSAGE("mu = 0: Anitescu-Tasora is Moreau-Jean", diff <= 1e-9);
}

/* Anitescu-Tasora, mu > 0: the square launched on the ground. In the first
   step it slides, and the cone complementarity gives u_N = mu |u_T| instead
   of 0: vy = mu vx with the normal impulse P = (g h + mu v0) / (1 + mu^2),
   the block lifts off (Moreau-Jean: vy = 0, see sliding_square) */
void separation_when_sliding(std::size_t scheme) {
  using namespace siconos;
  const double v0 = 2.;
  const double mu = 0.3;
  const double h = 5e-3;
  world w({.mu = mu, .h = h, .tmax = 0.5, .scheme = scheme});
  auto p = w.add_square(0., r, v0);
  w.finalize();
  w.step();
  const auto& v1 = storage::attr<"velocity">(p, w.simulation.current_step());
  const double vx1 = v1(0), vy1 = v1(1), om1 = v1(2);
  double max_lift = 0.;
  while (w.simulation.has_next_event()) {
    w.step();
    max_lift = std::max(max_lift, storage::attr<"q">(p, w.simulation.current_step())(1) - r);
  }
  std::print(
      stderr,
      "{}: sliding square, first step vx {:.9f} vy {:.9f} (mu vx {:.9f}), omega {:.2e}, "
      "max lift-off {:.3e} (h mu v0 = {:.3e})\n",
      name(scheme), vx1, vy1, mu * vx1, om1, max_lift, h * mu * v0);
  // the normal impulse P solves -g h + P = mu (v0 - mu P)
  const double impulse = (g * h + mu * v0) / (1 + mu * mu);
  CPPUNIT_ASSERT_MESSAGE("the cone complementarity: vy = mu vx",
                         std::abs(vy1 - mu * vx1) <= 1e-9);
  CPPUNIT_ASSERT_MESSAGE("vx = v0 - mu P", std::abs(vx1 - (v0 - mu * impulse)) <= 1e-9);
  CPPUNIT_ASSERT_MESSAGE("no rotation", std::abs(om1) <= 1e-9);
  CPPUNIT_ASSERT_MESSAGE("the block lifts off", max_lift > 0.1 * h * mu * v0);
}

/* Anitescu-Tasora, the warm start: a column of n disks at rest. The first
   step starts from lambda = 0 (cold), the next ones from the lambda of the
   previous step, close to the solution. Not at it: the sweeps stop when
   lambda changes less than the tolerance, and Gauss-Seidel converges slowly
   on a chain, so a small change still leaves an error that the next step
   corrects. 20 disks: PGS 10332 sweeps cold, at most 337 warm; APGD 404
   iterations cold, at most 42 warm */
void warm_started_column(std::size_t n, std::size_t scheme) {
  using namespace siconos;
  world w({.mu = 0.3, .h = 5e-3, .tmax = 0.5, .scheme = scheme});
  std::vector<decltype(w.add_disk(0., 0.))> disks;
  for (std::size_t i = 0; i < n; ++i) disks.push_back(w.add_disk(0., r + 2 * r * i));
  w.finalize();
  auto osnspb = w.simulation.one_step_nonsmooth_problem();
  std::size_t first = 0, later = 0;
  double max_pen = 0., max_v = 0.;
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    const std::size_t sweeps = storage::attr<"number_of_iterations">(osnspb);
    (s == 1 ? first : later) = std::max(s == 1 ? first : later, sweeps);
    double below = 0.;
    for (auto& d : disks) {
      const double y = storage::attr<"q">(d, s)(1);
      max_pen = std::max(max_pen, (below + r) - y);
      max_v = std::max(max_v, storage::attr<"velocity">(d, s).norm());
      below = y + r;
    }
  }
  std::print(stderr,
             "{}: column of {} disks at rest: {} iterations at the first step (cold), at most "
             "{} after (warm), max interpenetration {:.2e}, max |v| {:.2e}\n",
             name(scheme), n, first, later, max_pen, max_v);
  CPPUNIT_ASSERT_MESSAGE("the warm start needs far fewer iterations", 4 * later <= first);
  CPPUNIT_ASSERT_MESSAGE("no interpenetration in the column", max_pen <= 1e-10);
  CPPUNIT_ASSERT_MESSAGE("the column stays at rest", max_v <= 1e-10);
}

/* The explicit DEM scheme: k_n = 1e5 (m g / k_n ~ 1e-4), c_n for a damping
   ratio of 0.5 on a contact with the ground, k_t = 2/7 k_n, theta = 1 (the
   symplectic Euler scheme), h = 1e-4 (a tenth of the period of a contact) */
constexpr double kn = 1e5;
constexpr double kt = 2. / 7. * kn;
const double cn = 2 * 0.5 * std::sqrt(m * kn);

scene dem_scene(double mu, double h, double tmax, double alpha = 0., double damping = cn) {
  return {.mu = mu,
          .h = h,
          .tmax = tmax,
          .alpha = alpha,
          .scheme = schemes::explicit_penalty,
          .theta = 1.,
          .kn = kn,
          .cn = damping,
          .kt = kt};
}

/* DEM: a column of n disks at rest settles with the overlap of the i-th
   contact from the bottom (n - i) m g / k_n: the weight above it over the
   stiffness (Moreau-Jean: none, see column) */
void dem_column(std::size_t n) {
  using namespace siconos;
  world w(dem_scene(0.3, 1e-4, 1.));
  std::vector<decltype(w.add_disk(0., 0.))> disks;
  for (std::size_t i = 0; i < n; ++i) disks.push_back(w.add_disk(0., r + 2 * r * i));
  w.finalize();
  while (w.simulation.has_next_event()) w.step();
  const auto s = w.simulation.current_step();
  double below = 0., worst = 0., bottom = 0.;
  for (std::size_t i = 0; i < n; ++i) {
    const double y = storage::attr<"q">(disks[i], s)(1);
    const double overlap = (below + r) - y;
    const double expected = (n - i) * m * g / kn;
    worst = std::max(worst, std::abs(overlap - expected) / expected);
    if (i == 0) bottom = overlap;
    below = y + r;
  }
  std::print(stderr,
             "DEM: column of {} disks at rest: bottom overlap {:.6e} (n m g / k_n = {:.6e}), "
             "worst relative error on the overlaps {:.2e}\n",
             n, bottom, n * m * g / kn, worst);
  CPPUNIT_ASSERT_MESSAGE("the overlaps are the weights over k_n", worst <= 2e-2);
}

/* DEM: an undamped disk at rest on the ground. The symplectic Euler scheme
   is stable on the contact spring for h < 2 sqrt(m / k_n): below, the disk
   oscillates in the overlap [0, 2 m g / k_n] and stays on the ground;
   above, each contact gives energy and it bounces away */
void dem_stability() {
  using namespace siconos;
  const double h_crit = 2 * std::sqrt(m / kn);
  auto run = [](double h) {
    world w(dem_scene(0., h, 0.5, 0., 0.));
    auto d = w.add_disk(0., r);
    w.finalize();
    double lift = 0., overlap = 0.;
    while (w.simulation.has_next_event()) {
      w.step();
      const double y = storage::attr<"q">(d, w.simulation.current_step())(1);
      lift = std::max(lift, std::isfinite(y) ? y - r : INFINITY);
      overlap = std::max(overlap, r - y);
    }
    return std::pair{lift, overlap};
  };
  const auto [stable_lift, stable_overlap] = run(0.5 * h_crit);
  const auto [unstable_lift, unstable_overlap] = run(1.5 * h_crit);
  std::print(
      stderr,
      "DEM: disk at rest, undamped, h_crit = {:.4e}: at 0.5 h_crit max lift {:.2e}, max "
      "overlap {:.4e} (2 m g / k_n = {:.4e}); at 1.5 h_crit max lift {:.3e}, max overlap "
      "{:.3e}\n",
      h_crit, stable_lift, stable_overlap, 2 * m * g / kn, unstable_lift, unstable_overlap);
  CPPUNIT_ASSERT_MESSAGE(
      "stable below h_crit: in the overlap",
      stable_overlap <= 2.5 * m * g / kn && stable_lift <= 0.1 * m * g / kn);
  CPPUNIT_ASSERT_MESSAGE("unstable above h_crit: it bounces away",
                         unstable_lift >= 100 * m * g / kn);
}

/* DEM: the sliding square at its equilibrium overlap. While it slides the
   springs hold the Coulomb force: the stopping distance v0^2 / (2 mu g),
   up to the elastic offset of the springs and the normal vibrations */
void dem_sliding_square() {
  using namespace siconos;
  const double v0 = 2.;
  const double mu = 0.3;
  const double d_exact = v0 * v0 / (2 * mu * g);
  world w(dem_scene(mu, 1e-4, v0 / (mu * g) + 0.3));
  // two corners, each carrying m g / 2
  auto p = w.add_square(0., r - m * g / (2 * kn), v0);
  w.finalize();
  double max_dy = 0.;
  while (w.simulation.has_next_event()) {
    w.step();
    max_dy = std::max(max_dy, std::abs(storage::attr<"q">(p, w.simulation.current_step())(1) -
                                       (r - m * g / (2 * kn))));
  }
  const double d = storage::attr<"q">(p, w.simulation.current_step())(0);
  std::print(stderr,
             "DEM: sliding square: distance {:.6f} (exact {:.6f}, relative error {:.2e}), max "
             "|y - y_eq| {:.2e} (m g / (2 k_n) = {:.2e})\n",
             d, d_exact, (d - d_exact) / d_exact, max_dy, m * g / (2 * kn));
  CPPUNIT_ASSERT_MESSAGE("the stopping distance", std::abs(d - d_exact) <= 1e-2 * d_exact);
}

/* DEM: a square on an incline below the friction angle (tan(alpha) =
   0.36 < mu = 0.5). The tangential springs (two corners) hold m g
   sin(alpha) with an elastic offset x_eq = m g sin(alpha) / (2 k_t), but
   from unloaded springs the first oscillation overshoots: it needs twice
   that force, 2 tan(alpha) = 0.73 > mu, and the square slips. The
   oscillation then dies out (through the damped normal springs) and the
   square rests at x_eq plus the slip: 1.67e-4 for x_eq = 5.9e-5
   (Moreau-Jean: x = 0, see incline) */
void dem_incline_sticks() {
  using namespace siconos;
  const double alpha = 20 * deg;
  world w(dem_scene(0.5, 1e-4, 1., alpha));
  auto p = w.add_square(0., r - m * g * std::cos(alpha) / (2 * kn));
  w.finalize();
  const double x_eq = m * g * std::sin(alpha) / (2 * kt);
  // the range of x before and after t = 0.5
  double x_max = 0., late_min = INFINITY, late_max = -INFINITY;
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    const double x = storage::attr<"q">(p, s)(0);
    x_max = std::max(x_max, x);
    if (s * 1e-4 > 0.5) {
      late_min = std::min(late_min, x);
      late_max = std::max(late_max, x);
    }
  }
  const double x_end = storage::attr<"q">(p, w.simulation.current_step())(0);
  std::print(stderr,
             "DEM: incline 20 deg, mu 0.5 (sticks): x up to {:.3e} in the transient, at rest "
             "x(T) {:.4e} (elastic offset m g sin(alpha) / (2 k_t) = {:.3e}, slip {:.3e}), x "
             "range after t = 0.5 {:.2e}\n",
             x_max, x_end, x_eq, x_end - x_eq, late_max - late_min);
  CPPUNIT_ASSERT_MESSAGE("after the transient, no creep", late_max - late_min <= 1e-2 * x_eq);
  CPPUNIT_ASSERT_MESSAGE("it moved: at least the elastic offset", x_end >= 0.99 * x_eq);
}

class NscdCasesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(NscdCasesTest);
#ifndef NSCD_DEM_ONLY
  CPPUNIT_TEST(testSlidingSquare);
  CPPUNIT_TEST(testInclineSticks);
  CPPUNIT_TEST(testInclineSlides);
  CPPUNIT_TEST(testColumn);
  CPPUNIT_TEST(testAnitescuTasoraFrictionless);
  CPPUNIT_TEST(testAnitescuTasoraSeparation);
  CPPUNIT_TEST(testAnitescuTasoraWarmStart);
#endif
  CPPUNIT_TEST(testDemColumn);
  CPPUNIT_TEST(testDemStability);
  CPPUNIT_TEST(testDemFriction);
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
  void testAnitescuTasoraFrictionless() {
    for (auto scheme : {schemes::cone_complementarity, schemes::cone_complementarity_apgd}) {
      frictionless_column(5, scheme);
      for (double h : {1e-2, 5e-3}) incline(30 * deg, 0., h, scheme);
    }
  }
  void testAnitescuTasoraSeparation() {
    separation_when_sliding(schemes::cone_complementarity);
    separation_when_sliding(schemes::cone_complementarity_apgd);
  }
  void testAnitescuTasoraWarmStart() {
    warm_started_column(20, schemes::cone_complementarity);
    warm_started_column(20, schemes::cone_complementarity_apgd);
  }
  void testDemColumn() {
    for (std::size_t n : {1, 5}) dem_column(n);
  }
  void testDemStability() { dem_stability(); }
  void testDemFriction() {
    dem_sliding_square();
    dem_incline_sticks();
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(NscdCasesTest);

}  // namespace
