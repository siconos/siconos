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

#include <type_traits>

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/polytope.hpp"
#include "siconos/collision/shape_of.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* A rigid body whose shape kind is per body: declare it with the per-kind
   sets of shape_of.hpp,
     storage::dynamic_storage<storage::attached<body, symbol<"disk_shape">,
                                                some::item_ref<disk>>>
     storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                some::item_ref<polytope>>>
   (wrapped in storage::without_binding, so that the bindings do not expose
   their get_or_create accessors). The methods give the bindings a way to
   set the shape that keeps one shape kind per body. */
struct rigid_body : model::lagrangian_ds {
  template <typename Handle>
  struct interface : model::lagrangian_ds::template interface<Handle> {
    using model::lagrangian_ds::template interface<Handle>::self;

    template <typename SH>
    void set_shape(SH sh) {
      collision::set_shape(*self(), sh);
    }

    bool is_disk() { return collision::has_shape<shape::disk>(*self()); }
    bool is_polytope() { return collision::has_shape<shape::polytope>(*self()); }

    auto methods() {
      // types only: methods() is called on a bare interface object, without
      // a world (see storage::methods)
      using data_t = std::decay_t<decltype(self()->data())>;
      using indice = typename decltype(self()->env())::indice;
      using disk_h = storage::handle<storage::handle_base, shape::disk, indice, data_t>;
      using polytope_h =
          storage::handle<storage::handle_base, shape::polytope, indice, data_t>;
      return collect(
          method("set_disk_shape", &interface<Handle>::template set_shape<disk_h>),
          method("set_polytope_shape", &interface<Handle>::template set_shape<polytope_h>),
          method("is_disk", &interface<Handle>::is_disk),
          method("is_polytope", &interface<Handle>::is_polytope));
    }
  };
};

}  // namespace siconos::collision
