#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Signed distance from a point c to a convex polygon (vertices CCW).
   Returns dist >= 0, outside flag, outward normal (1->2, from polygon
   towards c when outside) and closest point on the polygon boundary. */
template <typename Scalar>
void point_polytope_distance(const std::array<Scalar, 2>& c,
                             const std::vector<std::array<Scalar, 2>>& v, std::size_t n,
                             Scalar& dist, bool& outside, std::array<Scalar, 2>& nrm,
                             std::array<Scalar, 2>& closest) {
  // inside test: for a CCW polygon, c is inside iff it is to the left of
  // every directed edge (cross(edge, c - vi) >= 0).
  bool inside = true;
  for (std::size_t i = 0; i < n && inside; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    const Scalar cr = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    if (cr < Scalar{0}) inside = false;
  }
  if (inside) {
    dist = Scalar{0};
    outside = false;
    closest = c;
    // outward normal of the nearest edge (right of edge dir for CCW)
    Scalar best = Scalar{0};
    bool have = false;
    for (std::size_t i = 0; i < n; ++i) {
      const auto& a = v[i];
      const auto& b = v[(i + 1) % n];
      const Scalar ex = b[0] - a[0];
      const Scalar ey = b[1] - a[1];
      // outward (right of edge dir) = (ey, -ex)
      const Scalar ox = ey;
      const Scalar oy = -ex;
      const Scalar onorm = std::hypot(ox, oy);
      if (onorm <= Scalar{0}) continue;
      const Scalar oox = ox / onorm;
      const Scalar ooy = oy / onorm;
      const Scalar dd = (c[0] - a[0]) * oox + (c[1] - a[1]) * ooy;
      const Scalar ad = std::fabs(dd);
      if (!have || ad < best) {
        best = ad;
        have = true;
        nrm = {oox, ooy};
        closest = {a[0] + dd * oox, a[1] + dd * ooy};
      }
    }
    if (!have) {
      nrm = {Scalar{1}, Scalar{0}};
      closest = c;
    }
    return;
  }
  // outside: closest point on polygon boundary (vertex or edge projection)
  dist = Scalar{0};
  outside = true;
  bool have = false;
  std::array<Scalar, 2> bestp{c[0], c[1]};
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    const Scalar ex = b[0] - a[0];
    const Scalar ey = b[1] - a[1];
    const Scalar ee = ex * ex + ey * ey;
    const Scalar ax = c[0] - a[0];
    const Scalar ay = c[1] - a[1];
    Scalar tt = (ax * ex + ay * ey) / ee;
    tt = std::min(Scalar{1}, std::max(Scalar{0}, tt));
    const Scalar px = a[0] + tt * ex;
    const Scalar py = a[1] + tt * ey;
    const Scalar dx = c[0] - px;
    const Scalar dy = c[1] - py;
    const Scalar dd = std::hypot(dx, dy);
    if (!have || dd < dist) {
      dist = dd;
      have = true;
      bestp = {px, py};
    }
  }
  closest = bestp;
  const Scalar nx = c[0] - bestp[0];
  const Scalar ny = c[1] - bestp[1];
  const Scalar nn = std::hypot(nx, ny);
  if (nn > Scalar{0})
    nrm = {nx / nn, ny / nn};
  else
    nrm = {Scalar{1}, Scalar{0}};
}

/* Polytope x fixed-disk contact relation (2D convex polygon vs a fixed
   disk/circle). Analytic: closed-form signed distance from the disk centre
   to the convex polygon, minus the disk radius. Consistent with the analytic
   diskfsegment_r/diskfdisk_r and with polytopefsegment_r (item 4); this avoids
   the still-open "circle shape in the cxd pool" ruling.

   Relation1: single dynamic body (the polytope). The fixed disk carries its
   own shape::disk (radius) and position q2. The cached status/h/normal/cp1
   /cp2 are computed per pair in evaluate().

   Convention (matches diskdisk_r/diskfsegment_r):
     h = +distance when separated / -depth when penetrating,
     normal points body1 (polytope) -> body2 (disk),
     cp1 on the polytope, cp2 on the disk surface.
   The signed distance from the disk centre c to the polytope is positive
   outside, zero on/inside the boundary; h = dist - r.
   When c is inside the polytope the pair reports as enclosed with depth = r
   (distance-based convention, same as diskfsegment_r). */

template <typename Topo>
struct polytopefdisk_r : item, model::relation1, model::any_lagrangian_relation {
  struct attributes {
    // the fixed disk this relation is bound to
    some::item_ref<shape::disk> disk;
    some::indice status;  // 0 separated, 1 enclosed (cxd::status value)
    some::scalar h;       // +distance / -depth
    some::vector<some::scalar, some::indice_value<2>> normal;  // 1 -> 2
    some::vector<some::scalar, some::indice_value<2>> cp1;     // on the polytope
    some::vector<some::scalar, some::indice_value<2>> cp2;     // on the disk
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) disk() {
      return storage::make_ref_handle(self()->data(), storage::attr<"disk">(*self()));
    }
    decltype(auto) shape() { return self()->disk(); }
    decltype(auto) status() { return storage::attr<"status">(*self()); }
    decltype(auto) h() { return storage::attr<"h">(*self()); }
    decltype(auto) normal() { return storage::attr<"normal">(*self()); }
    decltype(auto) cp1() { return storage::attr<"cp1">(*self()); }
    decltype(auto) cp2() { return storage::attr<"cp2">(*self()); }

    /* run the analytic disk x polytope gap and cache the result */
    template <match::handle<model::lagrangian_ds> DS1, match::handle<model::lagrangian_ds> DS2>
    void evaluate(auto step, DS1& ds1, DS2& ds2) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();

      auto& q1 = storage::attr<"q">(ds1, step);
      auto& q2 = storage::attr<"q">(ds2, step);
      auto sh = collision::shape_of<shape::polytope>(ds1);
      const scalar c = std::cos(q1(2)), s = std::sin(q1(2));
      const auto n = sh.nvertices();

      std::vector<std::array<scalar, 2>> v;
      v.reserve(n);
      for (std::size_t i = 0; i < n; ++i) {
        const auto& vi = sh.vertices()[i];
        v.push_back(std::array<scalar, 2>{q1(0) + c * vi[0] - s * vi[1],
                                          q1(1) + s * vi[0] + c * vi[1]});
      }

      // disk centre and radius
      const scalar cx = q2(0);
      const scalar cy = q2(1);
      auto dsk = self()->disk();
      const scalar r = dsk.radius();

      std::array<scalar, 2> cc{cx, cy};
      scalar dist;
      bool outside;
      std::array<scalar, 2> nrm, closest;
      point_polytope_distance(cc, v, n, dist, outside, nrm, closest);

      storage::attr<"h">(*self()) = dist - r;
      storage::attr<"status">(*self()) = outside ? 0 : 1;
      storage::attr<"normal">(*self()) = {nrm[0], nrm[1]};
      storage::attr<"cp1">(*self()) = {closest[0], closest[1]};
      // disk surface point: separated -> towards polytope (c - r*n);
      // enclosed -> in penetration direction (c + r*n)
      const scalar h = dist - r;
      std::array<scalar, 2> cp2;
      if (h > scalar{0})
        cp2 = {cx - r * nrm[0], cy - r * nrm[1]};
      else
        cp2 = {cx + r * nrm[0], cy + r * nrm[1]};
      storage::attr<"cp2">(*self()) = {cp2[0], cp2[1]};
    }

    template <match::handle<model::lagrangian_ds> DS1, match::handle<model::lagrangian_ds> DS2>
    [[gnu::always_inline]] decltype(auto) compute_h(auto step, DS1& ds1, DS2& ds2) {
      evaluate(step, ds1, ds2);
      return storage::attr<"h">(*self());
    }

    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::lagrangian_ds> DS2, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS1& ds1, DS2& ds2, M& g1) {
      evaluate(step, ds1, ds2);

      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = storage::attr<"q">(ds1, step);
      auto& n = storage::attr<"normal">(*self());
      auto& c1 = storage::attr<"cp1">(*self());

      const scalar nx = n[0];
      const scalar ny = n[1];
      const scalar tx = -ny;
      const scalar ty = nx;

      const scalar a1x = c1[0] - q1(0);
      const scalar a1y = c1[1] - q1(1);
      const scalar za1x = -a1y;
      const scalar za1y = a1x;

      g1(0, 0) = -nx;
      g1(0, 1) = -ny;
      g1(0, 2) = -(nx * za1x + ny * za1y);
      g1(1, 0) = -tx;
      g1(1, 1) = -ty;
      g1(1, 2) = -(tx * za1x + ty * za1y);
    }
  };
};

}  // namespace siconos::collision
