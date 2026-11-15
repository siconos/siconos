/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/chained_segment.hpp"

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
