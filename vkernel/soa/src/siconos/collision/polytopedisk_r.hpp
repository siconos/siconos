#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/polytopefdisk_r.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Disk x polytope contact relation: two dynamical systems, ds1 the disk, ds2
   the convex polygon (the space filter links a mixed pair in that order).
   Analytic: the signed distance from the disk centre to the polygon minus
   the disk radius, as polytopefdisk_r for a fixed disk. The disk side is
   curved, so the contact is one point.

   Convention (matches diskdisk_r):
     h = +distance when separated, -depth when penetrating, also when the
       disk centre is inside the polygon (depth = r + distance from the
       centre to the nearest edge);
     normal points body1 (disk) -> body2 (polytope);
     cp1 = c + r n on the disk, cp2 on the polygon boundary.
   G rows are the normal and the tangent t = (-n_y, n_x); the theta columns
   use the lever arms cp_i - c(q_i).

   One relation item serves every disk / polytope pair, so it holds no
   per-pair state: contact() evaluates the pair and returns the result. */
template <typename Topo>
struct polytopedisk_r : item, model::relation2, model::any_lagrangian_relation {
  struct attributes {};  // stateless

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

    template <match::handle<model::lagrangian_ds> DS1, match::handle<model::lagrangian_ds> DS2>
    auto contact(auto step, DS1& ds1, DS2& ds2) {
      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = storage::attr<"q">(ds1, step);
      const std::array<scalar, 2> centre{q1(0), q1(1)};
      const scalar r = collision::disk_radius(ds1);

      auto& q2 = storage::attr<"q">(ds2, step);
      auto sh = collision::shape_of<shape::polytope>(ds2);
      const scalar c = std::cos(q2(2)), s = std::sin(q2(2));
      const std::size_t n = sh.nvertices();
      std::array<std::array<scalar, 2>, 16> v;
      for (std::size_t i = 0; i < n; ++i) {
        const auto& vi = sh.vertices()[i];
        v[i] = {q2(0) + c * vi[0] - s * vi[1], q2(1) + s * vi[0] + c * vi[1]};
      }

      // pp.normal leaves the polygon towards the centre: polytope -> disk
      const auto pp = point_polytope_distance(centre, v, n);
      contact_t<scalar> ct;
      ct.h = pp.dist - r;
      ct.status = ct.h > scalar{0} ? 0 : 1;
      ct.normal = {-pp.normal[0], -pp.normal[1]};
      ct.cp1 = {centre[0] + r * ct.normal[0], centre[1] + r * ct.normal[1]};
      ct.cp2 = pp.closest;
      return ct;
    }

    template <match::handle<model::lagrangian_ds> DS1, match::handle<model::lagrangian_ds> DS2>
    [[gnu::always_inline]] auto compute_h(auto step, DS1& ds1, DS2& ds2) {
      return contact(step, ds1, ds2).h;
    }

    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::lagrangian_ds> DS2, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS1& ds1, DS2& ds2, M& g1, M& g2) {
      using scalar = typename decltype(self()->env())::scalar;
      const auto ct = contact(step, ds1, ds2);

      auto& q1 = storage::attr<"q">(ds1, step);
      auto& q2 = storage::attr<"q">(ds2, step);
      const scalar nx = ct.normal[0];
      const scalar ny = ct.normal[1];
      const scalar tx = -ny;
      const scalar ty = nx;

      // z ^ arm = (-arm_y, arm_x)
      const scalar za1x = -(ct.cp1[1] - q1(1));
      const scalar za1y = ct.cp1[0] - q1(0);
      const scalar za2x = -(ct.cp2[1] - q2(1));
      const scalar za2y = ct.cp2[0] - q2(0);

      g1(0, 0) = -nx;
      g1(0, 1) = -ny;
      g1(0, 2) = -(nx * za1x + ny * za1y);
      g1(1, 0) = -tx;
      g1(1, 1) = -ty;
      g1(1, 2) = -(tx * za1x + ty * za1y);

      g2(0, 0) = nx;
      g2(0, 1) = ny;
      g2(0, 2) = nx * za2x + ny * za2y;
      g2(1, 0) = tx;
      g2(1, 1) = ty;
      g2(1, 2) = tx * za2x + ty * za2y;
    }
  };
};

}  // namespace siconos::collision
