#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/cxd_narrowphase.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Polytope x fixed-segment contact relation (2D convex polygon vs a fixed
   segment). The segment is treated as a 2-vertex polytope with STATIC world
   vertices (no per-step transform) and is batched by the cxd_narrowphase
   stage alongside the polytope x polytope pairs.

   Structural reference: diskfsegment_r (item + model::relation1 +
   any_lagrangian_relation, single body, fixed segment stored as an
   item_ref<shape::segment>). The relation is batch-BACKED: the stage fills
   the cached status/h/normal/cp1/cp2, and compute_h/compute_jachq read the
   cache. An unbound narrowphase ref (standalone hand-cached relation) no-ops. */

template <typename Topo>
struct polytopefsegment_r : item, model::relation1, model::any_lagrangian_relation {
  using dof = some::indice_parameter<"dof">;
  using narrowphase_t = cxd_narrowphase<Topo>;

  // no attached properties: the segment ref lives on the relation side

  struct attributes {
    // loaded once at rebind time by the stage; unbound = standalone
    some::item_ref<narrowphase_t> narrowphase;
    // the fixed segment this relation is bound to
    some::item_ref<shape::segment> segment;
    some::indice status;  // 0 separated, 1 enclosed (cxd::status value)
    some::scalar h;       // +distance when separated / -depth when enclosed
    some::vector<some::scalar, some::indice_value<2>> normal;  // 1 -> 2 (polytope -> segment)
    some::vector<some::scalar, some::indice_value<2>> cp1;     // on the polytope
    some::vector<some::scalar, some::indice_value<2>> cp2;     // on the segment
    some::vector<some::scalar, some::indice_value<2>> init_dir;  // batch warm start
    some::indice cache_pair_checksum;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) status() { return storage::attr<"status">(*self()); }
    decltype(auto) h() { return storage::attr<"h">(*self()); }
    decltype(auto) narrowphase() { return storage::attr<"narrowphase">(*self()); }
    decltype(auto) segment() {
      return storage::make_ref_handle(self()->data(), storage::attr<"segment">(*self()));
    }
    decltype(auto) normal() { return storage::attr<"normal">(*self()); }
    decltype(auto) cp1() { return storage::attr<"cp1">(*self()); }
    decltype(auto) cp2() { return storage::attr<"cp2">(*self()); }
    decltype(auto) init_dir() { return storage::attr<"init_dir">(*self()); }
    decltype(auto) cache_pair_checksum() {
      return storage::attr<"cache_pair_checksum">(*self());
    }

    void ensure(auto step) {
      auto& data = self()->data();
      auto& ref = storage::attr<"narrowphase">(*self());
      using np_t = typename std::decay_t<decltype(ref)>::type;
      const auto n_stage = std::size(storage::handles<np_t>(data));
      if (ref.value() >= n_stage) return;  // unbound: standalone relation
      storage::make_ref_handle(data, ref).ensure(step);
    }

    template <match::handle<model::lagrangian_ds> DS>
    [[gnu::always_inline]] decltype(auto) compute_h(auto step, DS& ds) {
      ensure(step);
      return storage::attr<"h">(*self());
    }

    template <typename S, match::handle<model::lagrangian_ds> DS, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS& ds, M& g1) {
      ensure(step);

      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& q = storage::attr<"q">(ds, step);
      auto& n = storage::attr<"normal">(*self());
      auto& c1 = storage::attr<"cp1">(*self());

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
