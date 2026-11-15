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

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/refs.hpp"
#include "siconos/storage/remap_refs.hpp"

namespace siconos::storage {

static auto remove = [](auto& data, auto&& h) {
  using item_t = typename std::decay_t<decltype(h)>::type;
  static_assert(!match::not_removable<item_t>,
                "storage::remove: this item lives as long as the store (removable = false)");

  // walk the columns listed for item_t at registration, dynamic properties
  // included (their entries follow the swap and pop)
  auto h1 = h;  // clang 19
  if constexpr (has_del(h1)) {
    h1.__del__();
  }
  data.store()._registry.template remove_item<item_t>(h.index().value());
};

/**
 * Remove the items of Item of the given indices: as many swaps and pops,
 * with the refs to them rewritten once instead of once per removal. The indices are those
 * before the call; duplicates are ignored.
 */
template <typename Item>
void remove_many(auto& data, std::vector<std::size_t> indices) {
  static_assert(
      !match::not_removable<Item>,
      "storage::remove_many: this item lives as long as the store (removable = false)");
  std::sort(indices.begin(), indices.end(), std::greater<>{});
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  using indice = typename get_info_t<std::decay_t<decltype(data)>>::template env<Item>::indice;
  using handle_t = decltype(make_handle(data, index<Item, indice>{}));
  if constexpr (decltype(has_del(std::declval<handle_t&>()))::value) {
    for (auto i : indices) {
      auto h = make_handle(data, index<Item, indice>(static_cast<indice>(i)));
      h.__del__();
    }
  }
  data.store()._registry.template remove_items<Item>(std::move(indices));
}

}  // namespace siconos::storage
