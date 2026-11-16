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

#include <array>

#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"

namespace siconos::vkernel::storage {

using namespace pattern;

template <typename T, std::size_t N>
using memory_t = std::array<T, N>;

template <typename T>
constexpr bool is_std_array = false;
template <typename T, std::size_t N>
constexpr bool is_std_array<std::array<T, N>> = true;

static constexpr auto memory = []<typename T>(typename std::decay_t<T>::size_type step,
                                              T&& mem) constexpr -> decltype(auto) {
  if constexpr (is_std_array<std::decay_t<T>>) {
    // compile-time depth: the modulo folds away
    return static_cast<T&&>(mem)[step % std::size(static_cast<T&&>(mem))];
  } else {
    // runtime depth (registry columns): time_levels<> depths are 1 or 2, so avoid an
    // integer division at every access
    const auto n = std::size(static_cast<T&&>(mem));
    return static_cast<T&&>(mem)[n == 1 ? 0 : (n == 2 ? (step & 1) : step % n)];
  }
};

namespace memory_detail {
// the largest N of the time_levels<Attr, N> on Attr, 1 if none (several
// schemes may declare the levels they need, see simul/selectable_scheme.hpp):
// a constant-folded scan with __is_same builtins (a hana filter per attribute
// would instantiate one filter per column, over the whole list).  The leading
// dummy entry covers an empty list.
template <typename Attr, typename TimeLevels>
struct time_levels_size;

template <typename Attr, typename... Ks>
struct time_levels_size<Attr, boost::hana::tuple<Ks...>> {
  static constexpr std::size_t value = [] {
    constexpr bool same[] = {false, __is_same(Attr, typename Ks::type)...};
    constexpr std::size_t sizes[] = {1, Ks::size...};
    std::size_t levels = 1;  // not specified: one level
    for (std::size_t i = 1; i <= sizeof...(Ks); ++i) {
      if (same[i] && sizes[i] > levels) levels = sizes[i];
    }
    return levels;
  }();
};
}  // namespace memory_detail

template <match::attribute Attr, typename TimeLevels>
static constexpr std::size_t memory_size() {
  return memory_detail::time_levels_size<Attr, std::decay_t<TimeLevels>>::value;
};

}  // namespace siconos::vkernel::storage
