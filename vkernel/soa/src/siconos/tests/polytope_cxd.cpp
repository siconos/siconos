/* Unit tests for the cxd-backed polytope narrowphase (Tasks 2..6).
   Task 2: shape::polytope item — vertices, nvertices, rmax, convexity check. */

#include <cassert>
#include <cmath>
#include <print>

#include "siconos/siconos.hpp"

namespace siconos::config {
namespace some = siconos::storage::some;
struct polytope_ds : model::lagrangian_ds {};
using polytope_shape = collision::shape::polytope;
template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
struct data_t
    : storage::make<env, polytope_ds, polytope_shape,
                    storage::with_properties<
                        storage::wrapped<polytope_ds, some::unbounded_collection>,
                        storage::attached<polytope_ds, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<polytope_shape>>>> {};
}  // namespace siconos::config

#include <stdexcept>

/* ---- Task 3: narrowphase batch engine (exact pin values) ---- */
int test_polytope_item(siconos::config::data_t& data);

// clang-format off
#include "siconos/collision/cxd_narrowphase.hpp"
// clang-format on

namespace {

auto check = [](bool ok, auto&& what) {
  if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
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
  // two identical squares: enclosed, near-zero depth, finite outputs
  {
    auto in = make_input({square_at(0, 0), square_at(0, 0)}, {{0, 1}});
    auto out = siconos::collision::narrowphase_batch(std::move(in));
    check(out[0].status == 1, "identical squares enclosed");
    check(std::abs(out[0].h) <= 1e-6, "|h| < 1e-6");
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

int main() {
  siconos::config::data_t data;
  test_polytope_item(data);
  test_batch_engine();
}
