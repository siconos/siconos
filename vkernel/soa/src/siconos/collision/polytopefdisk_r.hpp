#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/collision/translated.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Signed distance from a point c to the boundary of a convex polygon of n
   vertices, of either orientation: dist > 0 outside, dist < 0 inside (minus
   the distance to the nearest edge). normal is the direction in which c
   leaves the polygon fastest -- from the closest boundary point towards c
   outside, the outward normal of the nearest edge inside -- and closest is
   that boundary point. */
template <typename Scalar>
struct point_polytope_result {
  Scalar dist;
  std::array<Scalar, 2> normal;
  std::array<Scalar, 2> closest;
};

template <typename Scalar, typename Vertices>
point_polytope_result<Scalar> point_polytope_distance(const std::array<Scalar, 2>& c,
                                                      const Vertices& v, std::size_t n) {
  // orientation: the outward normal of edge (a, b) is (ey, -ex) for a CCW
  // polygon and its opposite for a CW one
  Scalar area2 = Scalar{0};
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    area2 += a[0] * b[1] - b[0] * a[1];
  }
  const Scalar orient = area2 >= Scalar{0} ? Scalar{1} : Scalar{-1};

  bool inside = true;
  Scalar best_edge = Scalar{0};
  bool have_edge = false;
  point_polytope_result<Scalar> res{Scalar{0}, {Scalar{1}, Scalar{0}}, c};
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    const Scalar ex = b[0] - a[0];
    const Scalar ey = b[1] - a[1];
    const Scalar len = std::hypot(ex, ey);
    if (len <= Scalar{0}) continue;
    const Scalar ox = orient * ey / len;
    const Scalar oy = -orient * ex / len;
    // signed distance of c to the edge line, positive on the outer side
    const Scalar sd = (c[0] - a[0]) * ox + (c[1] - a[1]) * oy;
    if (sd > Scalar{0}) inside = false;
    if (!have_edge || -sd < best_edge) {
      best_edge = -sd;
      have_edge = true;
      res.normal = {ox, oy};
      res.closest = {c[0] - sd * ox, c[1] - sd * oy};
    }
  }
  if (inside) {
    res.dist = -best_edge;
    return res;
  }

  // outside: closest point of the boundary (vertex or edge projection)
  bool have = false;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    const Scalar ex = b[0] - a[0];
    const Scalar ey = b[1] - a[1];
    const Scalar ee = ex * ex + ey * ey;
    Scalar t = ee > Scalar{0} ? ((c[0] - a[0]) * ex + (c[1] - a[1]) * ey) / ee : Scalar{0};
    t = std::min(Scalar{1}, std::max(Scalar{0}, t));
    const Scalar px = a[0] + t * ex;
    const Scalar py = a[1] + t * ey;
    const Scalar d = std::hypot(c[0] - px, c[1] - py);
    if (!have || d < res.dist) {
      have = true;
      res.dist = d;
      res.closest = {px, py};
    }
  }
  if (res.dist > Scalar{0}) {
    res.normal = {(c[0] - res.closest[0]) / res.dist, (c[1] - res.closest[1]) / res.dist};
  }
  return res;
}

/* Polytope x fixed-disk contact relation (2D convex polygon vs a fixed
   disk). Analytic: the signed distance from the disk centre to the polygon
   minus the disk radius, as diskfdisk_r is for two disks.

   Relation1: the polytope is the only dynamical system; the fixed disk is a
   translated disk shape (centre = translation, radius of the disk).

   Convention (matches diskdisk_r / diskfsegment_r / diskfdisk_r):
     h = +distance when separated, -depth when penetrating, also when the
       disk centre is inside the polygon (depth = r + distance from the
       centre to the nearest edge);
     normal points body1 (polytope) -> body2 (disk);
     cp1 on the polytope boundary, cp2 = c - r n on the disk surface.

   One relation item serves every polytope near this disk, so it holds no
   per-pair state: contact() evaluates the pair and returns the result. */
template <typename Topo>
struct polytopefdisk_r : item, model::relation1, model::any_lagrangian_relation {
  struct attributes {
    some::item_ref<collision::translated<shape::disk>> translated_disk_shape;
  };

  template <typename Scalar>
  struct contact_t {
    int status;  // 0 separated, 1 penetrating (cxd::status values)
    Scalar h;
    std::array<Scalar, 2> normal;
    std::array<Scalar, 2> cp1;
    std::array<Scalar, 2> cp2;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) translated_disk_shape() {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"translated_disk_shape">(*self()));
    }
    decltype(auto) shape() { return self()->translated_disk_shape(); }

    template <match::handle<model::lagrangian_ds> DS>
    auto contact(auto step, DS& ds) {
      using scalar = typename decltype(self()->env())::scalar;

      auto& q = storage::attr<"q">(ds, step);
      auto sh = collision::shape_of<shape::polytope>(ds);
      const scalar c = std::cos(q(2)), s = std::sin(q(2));
      const std::size_t n = sh.nvertices();

      std::array<std::array<scalar, 2>, 16> v;
      for (std::size_t i = 0; i < n; ++i) {
        const auto& vi = sh.vertices()[i];
        v[i] = {q(0) + c * vi[0] - s * vi[1], q(1) + s * vi[0] + c * vi[1]};
      }

      auto fdisk = self()->translated_disk_shape();
      const auto& t = fdisk.translation();
      const std::array<scalar, 2> centre{t[0], t[1]};
      const scalar r = fdisk.translated().radius();

      const auto pp = point_polytope_distance(centre, v, n);
      contact_t<scalar> ct;
      ct.h = pp.dist - r;
      ct.status = ct.h > scalar{0} ? 0 : 1;
      ct.normal = pp.normal;
      ct.cp1 = pp.closest;
      ct.cp2 = {centre[0] - r * pp.normal[0], centre[1] - r * pp.normal[1]};
      return ct;
    }

    template <match::handle<model::lagrangian_ds> DS>
    [[gnu::always_inline]] auto compute_h(auto step, DS& ds) {
      return contact(step, ds).h;
    }

    template <typename S, match::handle<model::lagrangian_ds> DS, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS& ds, M& g1) {
      using scalar = typename decltype(self()->env())::scalar;
      const auto ct = contact(step, ds);

      auto& q = storage::attr<"q">(ds, step);
      const scalar nx = ct.normal[0];
      const scalar ny = ct.normal[1];
      const scalar tx = -ny;
      const scalar ty = nx;

      const scalar a1x = ct.cp1[0] - q(0);
      const scalar a1y = ct.cp1[1] - q(1);
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
