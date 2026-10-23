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
   cxd_narrowphase stage. Consistent with diskmesh_r / polytopemesh_r. */

template <typename Topo>
struct polytopefsegment_r : item,
                            mesh_relation,
                            model::relation1,
                            model::any_lagrangian_relation {
  struct attributes {
    // the fixed segment this relation is bound to
    some::item_ref<shape::segment> segment;
    some::indice status;
    some::scalar h;
    some::vector<some::scalar, some::indice_value<2>> normal;
    some::vector<some::scalar, some::indice_value<2>> cp1;
    some::vector<some::scalar, some::indice_value<2>> cp2;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) segment() {
      return storage::make_ref_handle(self()->data(), storage::attr<"segment">(*self()));
    }
    decltype(auto) shape() { return self()->segment(); }
    decltype(auto) status() { return storage::attr<"status">(*self()); }
    decltype(auto) h() { return storage::attr<"h">(*self()); }
    decltype(auto) normal() { return storage::attr<"normal">(*self()); }
    decltype(auto) cp1() { return storage::attr<"cp1">(*self()); }
    decltype(auto) cp2() { return storage::attr<"cp2">(*self()); }

    /* run cxd on this single polytope x segment pair and cache the result */
    template <match::handle<model::lagrangian_ds> DS>
    void evaluate(auto step, DS& ds) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();

      auto& q = storage::attr<"q">(ds, step);
      auto sh =
          storage::make_handle(data, storage::prop<"shape", storage::config_declared>(ds));
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

      auto out = narrowphase_batch(std::move(in));
      const auto& res = out[0];
      storage::attr<"status">(*self()) = res.status;
      storage::attr<"h">(*self()) = res.h;
      storage::attr<"normal">(*self()) = {res.normal[0], res.normal[1]};
      storage::attr<"cp1">(*self()) = {res.cp1[0], res.cp1[1]};
      storage::attr<"cp2">(*self()) = {res.cp2[0], res.cp2[1]};
    }

    template <match::handle<model::lagrangian_ds> DS>
    [[gnu::always_inline]] decltype(auto) compute_h(auto step, DS& ds) {
      evaluate(step, ds);
      return storage::attr<"h">(*self());
    }

    template <typename S, match::handle<model::lagrangian_ds> DS, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS& ds, M& g1) {
      evaluate(step, ds);

      using scalar = typename decltype(self()->env())::scalar;

      auto& q = storage::attr<"q">(ds, step);
      auto& n = storage::attr<"normal">(*self());
      auto& c1 = storage::attr<"cp1">(*self());
      auto& c2 = storage::attr<"cp2">(*self());

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
