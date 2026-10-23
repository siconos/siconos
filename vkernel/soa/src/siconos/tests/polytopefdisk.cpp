/* Unit test for polytopefdisk_r core geometry (polytope x fixed disk).

   Tests the analytic signed distance from a point to a convex polygon
   (point_polytope_distance), which is the core of compute_h. Storage-free.

   Square polytope centred at origin, side 2, vertices CCW:
     (1,-1), (1,1), (-1,1), (-1,-1)
   Disk centre (0,3), radius 1:
     dist((0,3), square) = 2  (to top edge y=1)
     h = dist - r = 1  (separated), normal (0,1), cp1=(0,1), cp2=(0,2).
   Disk centre (0,0), radius 1 (centre inside square):
     dist = 0, enclosed, h = -1, normal outward of nearest edge.
*/

#include <cmath>
#include <print>
#include <stdexcept>
#include <vector>

#include "siconos/collision/polytopefdisk_r.hpp"

namespace {
auto check = [](bool ok, auto&& what) {
  if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
};
auto near = [](double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; };
}  // namespace

int main() {
  using namespace siconos::collision;
  using Scalar = double;
  using Pt = std::array<Scalar, 2>;

  /* Square CCW: (1,-1), (1,1), (-1,1), (-1,-1) */
  std::vector<Pt> v = {{1, -1}, {1, 1}, {-1, 1}, {-1, -1}};
  std::size_t n = v.size();

  /* Case 1: disk centre (0,3) outside, radius 1 -> h = 1, separated */
  {
    Pt c = {0, 3};
    Scalar dist;
    bool outside;
    Pt nrm, closest;
    point_polytope_distance(c, v, n, dist, outside, nrm, closest);
    std::print("case1: dist={} outside={} nrm=({},{}) closest=({},{})\n", dist, outside,
               nrm[0], nrm[1], closest[0], closest[1]);
    check(near(dist, 2), "dist == 2");
    check(outside, "outside");
    check(near(nrm[0], 0) && near(nrm[1], 1), "normal (0,1)");
    check(near(closest[0], 0) && near(closest[1], 1), "closest (0,1)");
    Scalar h = dist - 1;
    check(near(h, 1), "h == 1");
  }

  /* Case 2: disk centre (0,0) inside square, radius 1 -> enclosed, h = -1 */
  {
    Pt c = {0, 0};
    Scalar dist;
    bool outside;
    Pt nrm, closest;
    point_polytope_distance(c, v, n, dist, outside, nrm, closest);
    std::print("case2: dist={} outside={} nrm=({},{}) closest=({},{})\n", dist, outside,
               nrm[0], nrm[1], closest[0], closest[1]);
    check(near(dist, 0), "dist == 0 (inside)");
    check(!outside, "inside");
    Scalar h = dist - 1;
    check(near(h, -1), "h == -1 (enclosed)");
  }

  std::print("polytopefdisk core geometry: all checks pass\n");
  return 0;
}
