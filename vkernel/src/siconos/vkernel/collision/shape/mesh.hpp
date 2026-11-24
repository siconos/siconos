// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/collision/collision.hpp"
#include "siconos/vkernel/collision/collision_head.hpp"
#include "siconos/vkernel/collision/shape/chained_segment.hpp"

namespace siconos::vkernel::collision::shape {

struct mesh : item {
  using items = gather<chained_segment>;
  // its chained_segment would have to go with it, and the segments refs of
  // the other meshes follow the one moved by the removal: it goes with the
  // store
  static constexpr bool removable = false;

  struct attributes {
    some::item_ref<chained_segment> segments;
    some::unbounded_vector<some::indice> global_indices;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__() {
      // Create the associated chained_segment
      auto& data = self()->data();
      auto seg = storage::add<chained_segment>(data);
      attr<"segments">(*self()) = seg.index();
    }

    decltype(auto) segments() {
      return storage::make_ref_handle(self()->data(), storage::attr<"segments">(*self()));
    }

    decltype(auto) global_indices() { return storage::attr<"global_indices">(*self()); }

    decltype(auto) set_nodes(auto& nodes) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      auto segments = self()->segments();

      // segments.nodes().clear();

      for (indice i = 0; i < algebra::rows(nodes); ++i) {
        segments.insert(nodes(i, 0), nodes(i, 1), nodes(i, 2));
      }
    }

    decltype(auto) indices() { return storage::attr<"indices">(*self()); }

    auto methods() {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using unbounded_matrix = typename env_t::template unbounded_matrix<scalar>;

      return collect(method("set_nodes", &interface<Handle>::set_nodes<unbounded_matrix>));
    }
  };
};

}  // namespace siconos::vkernel::collision::shape
