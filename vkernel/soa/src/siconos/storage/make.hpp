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

#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage {

namespace match = siconos::storage::pattern::match;

/**
 * @brief Applies a wrapper to a storage type based on the item type and
 * available wrappers
 *
 * @tparam Item The item type to find a wrapper for
 * @tparam Wrappers Tuple of wrapper definitions
 * @tparam Storage The underlying storage type
 * @param storage The storage instance to potentially wrap
 * @return The wrapped storage if a matching wrapper is found, otherwise the
 * original storage
 */
namespace wrapper_detail {
// Index of the first wrapped<> property for Item (sizeof...(Ws) if none), as a
// constant-folded scan with __is_same builtins, instead of a hana filter per
// column over the whole wrapper list.
template <typename Item, typename Wrappers>
struct wrapper_of;

template <typename Item, typename... Ws>
struct wrapper_of<Item, boost::hana::tuple<Ws...>> {
  static constexpr std::size_t index = [] {
    constexpr bool same[] = {false, __is_same(Item, typename Ws::type)...};
    for (std::size_t i = 1; i <= sizeof...(Ws); ++i) {
      if (same[i]) return i - 1;
    }
    return sizeof...(Ws);
  }();
  static constexpr bool found = index < sizeof...(Ws);
};
}  // namespace wrapper_detail

template <match::item Item, typename Wrappers, typename Storage>
static constexpr auto apply_wrapper(Storage storage) {
  using of = wrapper_detail::wrapper_of<Item, std::decay_t<Wrappers>>;
  // If matching wrapper found, apply it; otherwise return storage as-is
  if constexpr (of::found) {
    return []<typename... Ws>(boost::hana::tuple<Ws...>) {
      return typename std::tuple_element_t<of::index,
                                           std::tuple<Ws...>>::template wrapper<Storage>{};
    }(std::decay_t<Wrappers>{});
  } else {
    // without wrapper
    return Storage{};
  }
};

namespace build_detail {

/**
 * @brief Storage of one attribute column of Item, before its memory depth.
 *
 * Its three stages: attribute refinement, item wrapper or default storage,
 * and the time_levels<> depth (registry::column_of).
 */
template <typename Info, match::attribute A, typename Item>
struct column_storage {
  using all_properties_t = typename Info::all_properties_t;
  using item_t = Item;
  using env_t = typename Info::template env<item_t>;

  // attribute level: user refinements over the recursively rebuilt attribute
  using recursed_t = typename recursive_rebuild<A>::type;
  using refines_t =
      decltype(mp::filter(mp::filter(all_properties_t{}, mp::derive_from<property::refine>),
                          mp::is_inside_type_parent<recursed_t>));
  using refined_t = decltype([] {
    if constexpr (mp::size(refines_t{}) > mp::size_c<0_c>) {
      return typename nth_t<0, refines_t>::template refine<recursed_t>{};
    } else {
      return recursed_t{};
    }
  }());
  using value_t = typename traits::config<env_t>::template convert<refined_t>::type;

  // item level: item wrapper, wrapped<> property, or the env default storage
  using wrappers_t =
      decltype(mp::filter(all_properties_t{}, mp::derive_from<property::wrapped>));
  using storage_t = decltype([] {
    if constexpr (match::wrap<item_t>) {
      return typename traits::config<env_t>::template convert<
          typename item_t::template wrapper<value_t>>::type{};
    } else {
      using wrapped_t = decltype(apply_wrapper<item_t, wrappers_t>(value_t{}));
      if constexpr (!std::is_same_v<wrapped_t, value_t>) {
        return typename traits::config<env_t>::template convert<wrapped_t>::type{};
      } else {
        return typename env_t::template default_storage<value_t>{};
      }
    }
  }());

  // memory level: time_levels<> depth
  using time_levels_list =
      decltype(mp::filter(all_properties_t{}, mp::derive_from<property::time_levels>));
  static constexpr std::size_t depth = memory_size<A, time_levels_list>();
};

}  // namespace build_detail

}  // namespace siconos::storage
