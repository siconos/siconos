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

#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/mp/mp.hpp"

namespace siconos::vkernel::storage {

template <typename Handle>
struct default_interface {
  decltype(auto) self() {
    return static_cast<Handle*>(this);  // handle inherits from default_interface
  };

  auto env() {
    auto& data = self()->data();
    using info_t = get_info_t<decltype(data)>;
    return typename info_t::template env<typename Handle::type>{};
  }

  auto params() { return typename decltype(env())::params{}; }

  template <string_literal S>
  constexpr auto env_param() {
    return mp::get<pattern::param<S>>(self()->params()).value;
  }

  auto make_handle(auto& index) { return storage::make_handle(self()->data(), index); }
};

}  // namespace siconos::vkernel::storage
