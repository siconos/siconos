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
// translation only (for symmetric items)

#include "collision.hpp"
#include "collision_head.hpp"

namespace siconos::collision {

template <match::item Item>
struct translated : item {
  using item_t = Item;

  struct attributes {
    some::item_ref<item_t> translated;
    some::vector<some::scalar, some::indice_value<3>> translation;
    some::integer collision_group;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    // access to the translated shape
    decltype(auto) translated() {
      return storage::make_ref_handle(self()->data(), storage::attr<"translated">(*self()));
    };

    // translated points coordinates
    decltype(auto) point_coord(auto point_index) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using vector_t = typename env_t::template vector<scalar, 3>;

      vector_t coord = translation() + translated().point_coord(point_index);

      return coord;
    }

    decltype(auto) points_coords() {
      using env_t = decltype(self()->env());
      using indice_t = typename env_t::indice;

      return view::iota((indice_t)0, translated().maxpoints()) |
             view::transform([&](auto i) { return point_coord(i); });
    }

    decltype(auto) translation() { return storage::attr<"translation">(*self()); }
    decltype(auto) collision_group() { return storage::attr<"collision_group">(*self()); }
  };
};
}  // namespace siconos::collision
