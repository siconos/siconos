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

#include "siconos/vkernel/simul/simul_head.hpp"

namespace siconos::vkernel::simul {

template <typename SpaceFilter>
struct interaction_manager : item {
  using space_filter = SpaceFilter;
  using nslaw = typename space_filter::nslaw;

  using ncgroups = some::indice_parameter<"ncgroups">;

  using attributes =
      gather<attribute<"space_filter", some::item_ref<space_filter>>,
             attribute<"nslaws", some::matrix<some::item_ref<nslaw>, ncgroups, ncgroups>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) insert_nonsmooth_law(auto nslaw, auto gid1, auto gid2) {
      attr<"nslaws"> (*self())(gid1, gid2) = nslaw;
      return *self();
    }

    decltype(auto) get_nonsmooth_law(auto gid1, auto gid2) {
      auto& data = self()->data();
      return storage::make_handle(data, attr<"nslaws">(*self())(gid1, gid2));
    }

    void update_interactions(auto step) {
      auto& data = self()->data();
      storage::make_handle(data, attr<"space_filter">(*self())).update_index_set0(step);
    }

    auto methods() {
      using data_t = std::decay_t<decltype(self()->data())>;
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      // a handle, not an index: python passes handles, and pybind would
      // reinterpret a handle (data reference first) as an index
      using ref_nslaw_t = std::decay_t<decltype(storage::make_handle(
          std::declval<data_t&>(), storage::index<nslaw, indice>{}))>;

      return collect(
          method("insert_nonsmooth_law",
                 &interface<Handle>::insert_nonsmooth_law<ref_nslaw_t, indice, indice>),
          method("update_interactions", &interface<Handle>::update_interactions<indice>));
    }
  };
};
}  // namespace siconos::vkernel::simul
