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

#include <assert.h>

#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/properties.hpp"

namespace siconos::vkernel::storage {

using namespace siconos::vkernel::storage::pattern;

template <match::item Item>
decltype(auto) add(auto&& data) {
  // the columns of Item were listed at registration, so they are walked at
  // run time, with no compile-time scan of the world
  using indice = typename get_info_t<std::decay_t<decltype(data)>>::template env<Item>::indice;
  const auto index = static_cast<indice>(data.store()._registry.template push_item<Item>());
  auto h = make_full_handle<Item>(data, index);
  // handle initialization is not mandatory
  if constexpr (has_init(h)) {
    h.__init__();
  }
  return h;
}

// Multiple add at once.
template <match::item Item>
decltype(auto) add(auto&& data, std::size_t count) {
  return view::iota(std::size_t{0}, count) |
         view::transform([&data](auto&&) { return add<Item>(data); });
}

}  // namespace siconos::vkernel::storage
