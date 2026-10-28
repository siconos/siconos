/* Smoke test for the vendored cxd batch narrowphase API (externals/cxd).
   Two axis-aligned side-2 squares, one separated pair and one enclosed pair.
   Pool point type: std::array<double,2> — cxd's batch pool requires
   trivially-copyable vertices, which Eigen::Vector2d is not; the Eigen
   backend carries the algebra elsewhere. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>
#include <cassert>
#include <print>
#include <vector>

// clang-format off
/* load-bearing order — mirrors cxd's own test harness (test.cpp); do not
   regroup or sort: backends must be registered before cxd_batch.hpp pulls
   in cxd.hpp, which resolves cxd::num::inner_prod through them. */
#include "cxd_std.hpp"
#include "cxd_std_array_vec.hpp"
#include "cxd_ublas.hpp"
#include "cxd_eigen.hpp"
#include "cxd_batch.hpp"
// clang-format on

static std::vector<std::array<double, 2>> square(double cx, double cy) {
  // side 2, axis-aligned: (±1,±1) around (cx,cy), CCW
  return {{cx - 1, cy - 1}, {cx + 1, cy - 1}, {cx + 1, cy + 1}, {cx - 1, cy + 1}};
}

namespace {

class CxdBatchSmokeTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(CxdBatchSmokeTest);
  CPPUNIT_TEST(testCxdBatchSmoke);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testCxdBatchSmoke() {
    using P = std::array<double, 2>;
    std::vector<std::vector<P>> shapes = {square(0, 0), square(3, 0), square(1.9, 0)};
    auto pool = cxd::pack_shapes(shapes);  // pool ids 0,1,2
    std::vector<cxd::pair_record<double, 2>> pairs;
    pairs.push_back({0, 1, {1.0, 0.0}});  // separated: gap = 1 (touching at x=1..2)
    pairs.push_back({0, 2, {1.0, 0.0}});  // enclosed: overlap depth 0.1
    auto results = cxd::gjk_batch<2>(pool, pairs, {});
    CPPUNIT_ASSERT(results[0].st == cxd::status::separated);
    CPPUNIT_ASSERT(results[1].st == cxd::status::enclosed);
    for (auto& r : results) std::print("st={}\n", static_cast<int>(r.st));
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(CxdBatchSmokeTest);

}  // namespace
