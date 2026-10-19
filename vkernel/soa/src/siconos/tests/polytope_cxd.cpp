/* Unit tests for the cxd-backed polytope narrowphase (Tasks 2..6).
   Task 2: shape::polytope item — vertices, nvertices, rmax, convexity check. */

#include <cassert>
#include <cmath>
#include <print>

#include "siconos/siconos.hpp"

// clang-format off
/* before the config block: the world wires these types directly. The
   load-bearing cxd include order lives inside cxd_narrowphase.hpp; the
   relation header transitively brings the engine */
#include "siconos/collision/polytopexpolytope_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
struct polytope_ds : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using diskdisk_r = collision::diskdisk_r;
/* topo is fwd-declared as the ACYCLIC marker for the poly_r/narrowphase
   refs; defined right after interaction (Task 5 sim config, modeled on
   bouncing_disks) */
struct topo;
using polytopexpolytope_r = collision::polytopexpolytope_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct interaction : simul::interaction<nslaw, diskdisk_r, polytopexpolytope_r> {};
struct topo : simul::topology<polytope_ds, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
struct data_t
    : storage::make<env, polytope_ds, polytope_shape, diskdisk_r, nslaw, interaction, topo,
                    polytopexpolytope_r, cxd_narrowphase, osi, td, osnspb,
                    storage::with_properties<
                        storage::wrapped<polytope_ds, some::unbounded_collection>,
                        storage::wrapped<interaction, some::unbounded_collection>,
                        storage::attached<polytope_ds, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<polytope_shape>>>> {};
}  // namespace siconos::config

#include <stdexcept>

/* ---- Task 3: narrowphase batch engine (exact pin values) ---- */
int test_polytope_item(siconos::config::data_t& data);

namespace {

auto check = [](bool ok, auto&& what) {
  if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
};
// relation vector attrs come back as Eigen 2-vectors; the near() helpers
// below take std::array, so funnel them through this
auto to_arr = [](const auto& v) {
  return std::array<double, 2>{static_cast<double>(v[0]), static_cast<double>(v[1])};
};

std::vector<std::array<double, 2>> square_at(double cx, double cy) {
  return {{cx - 1, cy - 1}, {cx + 1, cy - 1}, {cx + 1, cy + 1}, {cx - 1, cy + 1}};
}
std::vector<std::array<double, 2>> diamond_at_origin() {
  const double r = 1.4142135623730951;
  return {{r, 0}, {0, r}, {-r, 0}, {0, -r}};
}
auto near = [](std::array<double, 2> a, std::array<double, 2> e, double tol = 1e-9) {
  return std::abs(a[0] - e[0]) <= tol && std::abs(a[1] - e[1]) <= tol;
};
auto nears = [](double a, double e, double tol = 1e-9) { return std::abs(a - e) <= tol; };

siconos::collision::narrowphase_input make_input(
    std::vector<std::vector<std::array<double, 2>>> shapes,
    std::vector<std::array<std::uint32_t, 2>> pairs_in) {
  siconos::collision::narrowphase_input in;
  for (auto& s : shapes) {
    in.entries.push_back({static_cast<std::uint32_t>(in.vertices.size()),
                          static_cast<std::uint32_t>(s.size())});
    in.vertices.insert(in.vertices.end(), s.begin(), s.end());
  }
  for (auto [a, b] : pairs_in) {
    in.ids_a.push_back(a);
    in.ids_b.push_back(b);
    in.init_dirs.push_back({1.0, 0.0});
  }
  return in;
}
}  // namespace

int test_polytope_item(siconos::config::data_t& data) {
  using siconos::storage::add;
  // side-2 axis-aligned square, CCW; assert() is compiled out under Release,
  // an explicit throwing check keeps RED/GREEN honest
  auto check = [](bool ok, auto&& what) {
    if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
  };
  auto sh = add<siconos::config::polytope_shape>(data);
  std::array<std::array<double, 2>, 4> sq = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
  for (int i = 0; i < 4; ++i) sh.vertices()[i] = {sq[i][0], sq[i][1]};
  sh.nvertices() = 4;
  check(std::abs(sh.rmax() - std::sqrt(2.0)) < 1e-12, "rmax == sqrt(2)");

  // non-convex polygon rejected by __init__
  auto sh2 = add<siconos::config::polytope_shape>(data);
  std::array<std::array<double, 2>, 5> arrow = {{{0, 0}, {2, 0}, {1, 1}, {2, 2}, {0, 2}}};
  for (int i = 0; i < 5; ++i) sh2.vertices()[i] = {arrow[i][0], arrow[i][1]};
  sh2.nvertices() = 5;
  {
    bool threw = false;
    try {
      sh2.__init__();
    } catch (const std::runtime_error&) {
      threw = true;
    }
    check(threw, "non-convex polygon rejected by __init__");
  }
  std::print("polytope item ok, rmax={}\n", sh.rmax());
  return 0;
}

void test_batch_engine() {
  // {0,0} square x {3,0} square: separated
  {
    auto in = make_input({square_at(0, 0), square_at(3, 0)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(out.size() == 1, "one pair out");
    check(out[0].status == 0 /* separated */, "s(0,0)-s(3,0) separated");
    check(nears(out[0].h, 1.0), "h == 1");
    check(near(out[0].cp1, {1, 0}), "cp1 == (1,0)");
    check(near(out[0].cp2, {2, 0}), "cp2 == (2,0)");
    check(near(out[0].normal, {1, 0}), "normal == (1,0)");
  }
  // {0,0} square x {3,3} square: separated on the diagonal
  {
    auto in = make_input({square_at(0, 0), square_at(3, 3)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(nears(out[0].h, std::sqrt(2.0)), "h == sqrt(2)");
    check(near(out[0].cp1, {1, 1}), "cp1 == (1,1)");
    check(near(out[0].cp2, {2, 2}), "cp2 == (2,2)");
    const double c = 0.7071067811865476;
    check(near(out[0].normal, {c, c}), "normal == diag");
  }
  // diamond x square(3,0): separated, vertex-edge closest pair
  {
    auto in = make_input({diamond_at_origin(), square_at(3, 0)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(out[0].status == 0, "diamond-square separated");
    check(nears(out[0].h, 2 - 1.4142135623730951, 1e-9), "h == 2-sqrt(2)");
    check(near(out[0].cp1, {1.4142135623730951, 0}), "cp1 == diamond east");
    check(near(out[0].cp2, {2, 0}), "cp2 == (2,0)");
    check(near(out[0].normal, {1, 0}), "normal == (1,0)");
  }
  // {0,0} square x {1.9,0} square: enclosed, overlap depth 0.1
  {
    auto in = make_input({square_at(0, 0), square_at(1.9, 0)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(std::abs(out[0].h + 0.1) <= 1e-6, "h == -0.1");
    check(near(out[0].cp1, {1, 0}, 1e-6), "cp1 == (1,0)");
    check(near(out[0].cp2, {0.9, 0}, 1e-6), "cp2 == (0.9,0)");
    check(near(out[0].normal, {1, 0}, 1e-6), "normal == (1,0)");
  }
  // two identical squares: enclosed, and the depth is the minimum
  // translation that separates them (the half-width of the difference
  // square, 2) -- EPA answers it, so the h is -2
  {
    auto in = make_input({square_at(0, 0), square_at(0, 0)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(out[0].status == 1, "identical squares enclosed");
    check(std::abs(out[0].h + 2.0) <= 1e-6, "h == -2");
    for (double x : {out[0].h, out[0].normal[0], out[0].normal[1], out[0].cp1[0],
                     out[0].cp1[1], out[0].cp2[0], out[0].cp2[1]})
      check(std::isfinite(x), "all outputs finite");
  }
  // one batch call over {separated, enclosed}: statuses preserved
  {
    auto in =
        make_input({square_at(0, 0), square_at(3, 0), square_at(1.9, 0)}, {{0, 1}, {0, 2}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(out.size() == 2, "two pairs out");
    check(out[0].status == 0 && out[1].status == 1, "statuses preserved");
    check(nears(out[1].h, -0.1, 1e-6), "second pair h == -0.1");
  }
  std::print("narrowphase_batch engine: all exact cases pass\n");
}

/* ---- Task 4: polytopexpolytope_r relation, exact G pins ----
   Standalone world: no narrowphase stage item exists, the relation's cached
   h/normal/cps are set manually and ensure() must no-op (ledger ruling). */
void test_polytopexpolytope_relation(siconos::config::data_t& data) {
  using siconos::storage::add;

  // body 1 center (0,0), body 2 center (0,1.9): flush vertical overlap of
  // 0.1 -> n = (0,1), cp1 = (0,1), cp2 = (0,0.9)
  auto ds1 = add<siconos::config::polytope_ds>(data);
  auto ds2 = add<siconos::config::polytope_ds>(data);
  ds1.q() = {0, 0, 0};
  ds2.q() = {0, 1.9, 0};

  auto rel = add<siconos::config::polytopexpolytope_r>(data);
  // manual cache (the narrowphase ref stays unbound)
  rel.status() = 1;
  rel.h() = -0.1;
  rel.normal() = {0, 1};
  rel.cp1() = {0, 1};
  rel.cp2() = {0, 0.9};
  rel.init_dir() = {0, -1};

  check(std::abs(rel.compute_h(0, ds1, ds2) + 0.1) <= 1e-12, "compute_h == -0.1");

  /* hand-derived (plan Task 4, with one corrected pin): n = (0,1),
     t = z^n = (-1,0); arm1 = cp1 - c1 = (0,1) -> z^arm1 = (-1,0);
     arm2 = cp2 - c2 = (0,0.9) - (0,1.9) = (0,-1.0) -> z^arm2 = (1.0,0).
     (the plan text pins arm2 = (0,-1.1): arithmetic slip, 0.9 - 1.9 = -1.0;
     the structural rule arm_i = cp_i - c(q_i) is unchanged — ledger) */
  Eigen::Matrix<double, 2, 3> g1 = Eigen::Matrix<double, 2, 3>::Zero();
  Eigen::Matrix<double, 2, 3> g2 = Eigen::Matrix<double, 2, 3>::Zero();
  rel.compute_jachq(0, ds1, ds2, g1, g2);
  const double pin_g1[2][3] = {{0, -1, 0}, {1, 0, -1}};
  const double pin_g2[2][3] = {{0, 1, 0}, {-1, 0, -1.0}};
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 3; ++j) {
      check(std::abs(g1(i, j) - pin_g1[i][j]) <= 1e-12, "G1 pin");
      check(std::abs(g2(i, j) - pin_g2[i][j]) <= 1e-12, "G2 pin");
    }
  }
  std::print("polytopexpolytope_r: G exact pins pass\n");
}

/* ---- Task 5: cxd_narrowphase stage — rebind + batch + write-back + cache ---- */
void test_narrowphase_stage(siconos::config::data_t& data) {
  using siconos::storage::add;

  // body1: axis-aligned square, body2: diamond r=1 (CCW locals)
  auto sh1 = add<siconos::config::polytope_shape>(data);
  std::array<std::array<double, 2>, 4> sq = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
  for (int i = 0; i < 4; ++i) sh1.vertices()[i] = {sq[i][0], sq[i][1]};
  sh1.nvertices() = 4;
  auto sh2 = add<siconos::config::polytope_shape>(data);
  std::array<std::array<double, 2>, 4> dm = {{{1, 0}, {0, 1}, {-1, 0}, {0, -1}}};
  for (int i = 0; i < 4; ++i) sh2.vertices()[i] = {dm[i][0], dm[i][1]};
  sh2.nvertices() = 4;

  // body1 center (0,0), body2 center (0,1.9): overlap 0.1 along y
  auto ds1 = add<siconos::config::polytope_ds>(data);
  ds1.q() = {0, 0, 0};
  siconos::storage::prop<"shape">(ds1) = sh1;
  auto ds2 = add<siconos::config::polytope_ds>(data);
  ds2.q() = {0, 1.9, 0};
  siconos::storage::prop<"shape">(ds2) = sh2;

  // generic space_filter-style binding: topo.link + generic diskdisk_r
  auto topo_h = add<siconos::config::topo>(data);
  auto dsds = add<siconos::config::diskdisk_r>(data);
  auto inter = topo_h.link(ds1, ds2);
  inter.relation() = dsds;

  // the narrowphase stage
  auto np = add<siconos::config::cxd_narrowphase>(data);
  np.topology() = topo_h;
  np.dsds_r() = dsds;

  np.ensure(0);

  // the interaction's relation is now a polytopexpolytope_r with the exact
  // cached pair (square x diamond, flush y-overlap 0.1)
  using poly_r_index =
      siconos::storage::index<siconos::config::polytopexpolytope_r, std::size_t>;
  auto check_rebound = [&]() {
    check(std::holds_alternative<poly_r_index>(inter.relation()),
          "relation rebound to polytopexpolytope_r");
    return siconos::storage::make_handle(data, std::get<poly_r_index>(inter.relation()));
  };
  {
    auto rel = check_rebound();
    check(rel.status() == 1 /* enclosed */, "rebound relation enclosed");
    check(std::abs(rel.h() + 0.1) <= 1e-6, "rebound h == -0.1");
    check(near(to_arr(rel.normal()), {0, 1}, 1e-6), "rebound normal == (0,1)");
    check(near(to_arr(rel.cp1()), {0, 1}, 1e-6), "rebound cp1 == (0,1)");
    check(near(to_arr(rel.cp2()), {0, 0.9}, 1e-6), "rebound cp2 == (0,0.9)");
    // the warm start written back for the next batch is the contact normal
    check(near(to_arr(rel.init_dir()), {0, 1}, 1e-6), "init_dir == the contact normal");
  }

  // mutate q2 at the SAME step tag: the checksum must re-run the batch
  ds2.q() = {0, 1.95, 0};
  np.ensure(0);
  {
    auto rel = check_rebound();
    check(std::abs(rel.h() + 0.05) <= 1e-6, "checksum re-run: h == -0.05");
  }

  // rotation convention pin: local (1,0) under theta=pi/2 -> world (0,1) CCW
  ds2.q() = {0, 0, M_PI / 2};
  const auto wv = np.world_vertices(ds2, 0);
  check(wv.size() == 4, "four world vertices");
  check(near(wv[0], {0, 1}), "R(pi/2)*(1,0) == (0,1) (CCW convention)");
  // and the rotated diamond still reads the same overlap when centered on
  // body1's box gap... verify invariance with the square back at (0,1.9)
  ds2.q() = {0, 1.9, M_PI / 2};
  np.ensure(0);
  {
    auto rel = check_rebound();
    check(std::abs(rel.h() + 0.1) <= 1e-6, "rotated pair h == -0.1");
    check(near(to_arr(rel.normal()), {0, 1}, 1e-6), "rotated pair normal == (0,1)");
  }
  ds2.q() = {0, 1.9, 0};

  std::print("cxd_narrowphase stage: rebound, cached exact pair, checksum re-runs\n");
}

int main() {
  siconos::config::data_t data;
  test_polytope_item(data);
  test_batch_engine();
  test_polytopexpolytope_relation(data);
  test_narrowphase_stage(data);
}
