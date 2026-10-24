#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/cxd_narrowphase.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Polytope x fixed-segment contact relation (2D convex polygon vs a fixed
   segment). The segment is treated as a 2-vertex polytope with STATIC world
   vertices (no per-step transform). Analytic-ish: geometry is evaluated per
   pair via a one-pair cxd batch, so this relation does not depend on the
   cxd_narrowphase stage. Consistent with diskmesh_r / polytopemesh_r.
   One relation item per (segment, manifold point) serves every polytope near
   the segment, so it holds no per-pair state: contact() evaluates the pair
   and returns the point_index-th point of its contact manifold. */

template <typename Topo>
struct polytopefsegment_r : item, model::relation1, model::any_lagrangian_relation {
  struct attributes {
    // the fixed segment this relation is bound to
    some::item_ref<shape::segment> segment;
    // which point of the contact manifold (0 or 1) this relation stands for:
    // a body near the segment has one interaction per manifold point
    some::indice point_index;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) segment() {
      return storage::make_ref_handle(self()->data(), storage::attr<"segment">(*self()));
    }
    decltype(auto) shape() { return self()->segment(); }
    decltype(auto) point_index() { return storage::attr<"point_index">(*self()); }

    /* run cxd on this single polytope x segment pair */
    template <match::handle<model::lagrangian_ds> DS>
    auto contact(auto step, DS& ds) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();

      auto& q = storage::attr<"q">(ds, step);
      auto sh = collision::shape_of<shape::polytope>(ds);
      const scalar c = std::cos(q(2)), s = std::sin(q(2));
      const auto n = sh.nvertices();

      narrowphase_input<scalar> in;
      // entry 0: polytope (world vertices)
      const std::uint32_t off0 = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const auto& v = sh.vertices()[i];
        in.vertices.push_back(
            std::array<scalar, 2>{q(0) + c * v[0] - s * v[1], q(1) + s * v[0] + c * v[1]});
      }
      in.entries.push_back({off0, static_cast<std::uint32_t>(n)});
      const std::uint32_t slot0 = 0;
      // entry 1: segment endpoints (static world coords)
      auto seg = self()->segment();
      const std::uint32_t off1 = static_cast<std::uint32_t>(in.vertices.size());
      in.vertices.push_back(std::array<scalar, 2>{seg.x1(), seg.y1()});
      in.vertices.push_back(std::array<scalar, 2>{seg.x2(), seg.y2()});
      in.entries.push_back({off1, 2});
      const std::uint32_t slot1 = 1;

      in.ids_a.push_back(slot0);
      in.ids_b.push_back(slot1);
      in.init_dirs.push_back({1.0, 0.0});

      return manifold_point(narrowphase_batch(std::move(in))[0], point_index());
    }

    template <match::handle<model::lagrangian_ds> DS>
    [[gnu::always_inline]] auto compute_h(auto step, DS& ds) {
      return contact(step, ds).h;
    }

    template <typename S, match::handle<model::lagrangian_ds> DS, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS& ds, M& g1) {
      const auto ct = contact(step, ds);

      using scalar = typename decltype(self()->env())::scalar;

      auto& q = storage::attr<"q">(ds, step);
      const auto& n = ct.normal;
      const auto& c1 = ct.cp1;

      const scalar nx = n[0];
      const scalar ny = n[1];
      const scalar tx = -ny;
      const scalar ty = nx;

      const scalar a1x = c1[0] - q(0);
      const scalar a1y = c1[1] - q(1);
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
