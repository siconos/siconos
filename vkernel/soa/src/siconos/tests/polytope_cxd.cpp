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

int main() {
  using namespace siconos;
  config::data_t data;

  // assert() is compiled out under Release (-DNDEBUG); an explicit throwing
  // check keeps RED/GREEN honest.
  auto check = [](bool ok, auto&& what) {
    if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
  };

  // side-2 axis-aligned square, CCW
  auto sh = storage::add<config::polytope_shape>(data);
  std::array<std::array<double, 2>, 4> sq = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
  for (int i = 0; i < 4; ++i) {
    sh.vertices()[i] = {sq[i][0], sq[i][1]};
  }
  sh.nvertices() = 4;
  check(std::abs(sh.rmax() - std::sqrt(2.0)) < 1e-12, "rmax == sqrt(2)");  // exact
  std::print("polytope item ok, rmax={}\n", sh.rmax());

  // convexity check: a non-convex polygon must make __init__ throw
  auto sh2 = storage::add<config::polytope_shape>(data);
  std::array<std::array<double, 2>, 5> arrow = {
      {{0, 0}, {2, 0}, {1, 1}, {2, 2}, {0, 2}}};  // arrowhead: sign flip
  for (int i = 0; i < 5; ++i) {
    sh2.vertices()[i] = {arrow[i][0], arrow[i][1]};
  }
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
  std::print("polytope item ok, rmax={}, convexity check ok\n", sh.rmax());
}
