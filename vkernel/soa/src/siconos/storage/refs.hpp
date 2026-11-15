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
#include <atomic>
#include <cstddef>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "siconos/storage/sparse_set.hpp"

/**
 * References between items held in the store, and their fix-up when items
 * move: storage::remove (swap and pop) and storage::permute move the rows of
 * an item, the refs to these rows held in other columns follow through
 * remap_refs.
 *
 * A ref is a storage::index<T, R> (some::item_ref<T> after conversion). It
 * may sit in a plain value, a variant (polymorphic attribute), a tuple or a
 * pair, a range (the per-step storages, std::vector, small_vector,
 * std::array), the mapped values of a map, the values of a sparse set
 * (dynamic property) or an Eigen matrix. The refs held outside the store (the
 * handles, the maps of the user code, the indices kept as plain integers) do
 * not follow.
 */
namespace siconos::vkernel::storage {

template <typename... Ts>
struct type_list {};

namespace registry_detail {
inline std::size_t next_slot() {
  // atomic: two threads may meet two new keys for the first time at once
  static std::atomic<std::size_t> n = 0;
  return n.fetch_add(1, std::memory_order_relaxed);
}
}  // namespace registry_detail

/// Dense process-wide id of a column key, assigned on first use (entt's
/// type_index pattern). Shared by every translation unit, so a registry filled
/// in one TU is readable from another.
template <typename Key>
std::size_t slot_id() {
  static const std::size_t id = registry_detail::next_slot();
  return id;
}

/// the kind of an item: refs compare their target kind with the kind of the
/// moved item
template <typename Item>
struct item_kind {};

/// the process-wide id of the kind of item T (the registry's kind slot)
template <typename T>
struct kind_id {
  static std::size_t get() { return slot_id<item_kind<T>>(); }
};

namespace refs_detail {

template <typename... Ls>
struct concat;
template <>
struct concat<> {
  using type = type_list<>;
};
template <typename... As>
struct concat<type_list<As...>> {
  using type = type_list<As...>;
};
template <typename... As, typename... Bs, typename... Ls>
struct concat<type_list<As...>, type_list<Bs...>, Ls...> {
  using type = typename concat<type_list<As..., Bs...>, Ls...>::type;
};

template <typename T>
concept index_like = requires {
  typename T::index_t;
  typename T::type;
  typename T::value_t;
} && requires(T t) { t.value(); };

template <typename T>
concept eigen_like = requires(T t) {
  typename T::Scalar;
  t.data();
  t.size();
};

template <typename T>
concept map_like = requires {
  typename T::key_type;
  typename T::mapped_type;
};

template <typename T>
concept range_like = requires(T t) {
  typename T::value_type;
  t.begin();
  t.end();
} && !std::is_same_v<typename T::value_type, T>;

template <typename T>
struct is_variant : std::false_type {};
template <typename... Ts>
struct is_variant<std::variant<Ts...>> : std::true_type {};

template <typename T>
struct is_tuple : std::false_type {};
template <typename... Ts>
struct is_tuple<std::tuple<Ts...>> : std::true_type {};
template <typename A, typename B>
struct is_tuple<std::pair<A, B>> : std::true_type {};

template <typename V>
struct targets;

template <typename V>
using targets_t = typename targets<std::remove_cv_t<V>>::type;

template <typename... Ts>
struct targets_of_all {
  using type = typename concat<targets_t<Ts>...>::type;
};

template <typename V>
struct targets {
  static constexpr auto get() {
    if constexpr (index_like<V>) {
      return type_list<typename V::type>{};
    } else if constexpr (is_sparse_set_v<V>) {
      return targets_t<typename V::value_type>{};
    } else if constexpr (eigen_like<V>) {
      return targets_t<typename V::Scalar>{};
    } else if constexpr (map_like<V>) {
      return targets_t<typename V::mapped_type>{};
    } else if constexpr (range_like<V>) {
      return targets_t<typename V::value_type>{};
    } else {
      return type_list<>{};
    }
  }
  using type = decltype(get());
};

template <typename... Ts>
struct targets<std::variant<Ts...>> {
  using type = typename targets_of_all<Ts...>::type;
};
template <typename... Ts>
struct targets<std::tuple<Ts...>> {
  using type = typename targets_of_all<Ts...>::type;
};
template <typename A, typename B>
struct targets<std::pair<A, B>> {
  using type = typename targets_of_all<A, B>::type;
};

template <typename V>
struct has_no_key_refs : std::true_type {};
template <map_like V>
struct has_no_key_refs<V> : std::is_same<targets_t<typename V::key_type>, type_list<>> {};

}  // namespace refs_detail

/// the item types V refers to (with repetitions), as a type_list
template <typename V>
using ref_targets_t = refs_detail::targets_t<V>;

/// whether V holds refs to some item
template <typename V>
inline constexpr bool holds_refs_v = !std::is_same_v<ref_targets_t<V>, type_list<>>;

namespace refs_detail {
template <typename T, typename L>
struct contains;
template <typename T, typename... Ts>
struct contains<T, type_list<Ts...>> : std::bool_constant<(std::is_same_v<T, Ts> || ...)> {};
}  // namespace refs_detail

/// whether V holds refs to items of type Item
template <typename V, typename Item>
inline constexpr bool refers_to_v = refs_detail::contains<Item, ref_targets_t<V>>::value;

/**
 * How the rows of an item kind moved. Swap and pop (storage::remove): the
 * row removed goes away, the row last takes its index. Table
 * (storage::permute, a batch removal): row k goes to old_to_new[k], which is
 * the invalid index for a removed row.
 */
struct ref_remap {
  static constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

  std::size_t removed = invalid;
  std::size_t last = invalid;
  const std::vector<std::size_t>* old_to_new = nullptr;

  static ref_remap swap_and_pop(std::size_t removed, std::size_t last) {
    return {removed, last, nullptr};
  }
  static ref_remap table(const std::vector<std::size_t>& old_to_new) {
    return {invalid, invalid, &old_to_new};
  }

  std::size_t operator()(std::size_t k) const {
    if (old_to_new) {
      return k < old_to_new->size() ? (*old_to_new)[k] : k;
    }
    if (k == removed) return invalid;
    if (k == last) return removed;
    return k;
  }
};

namespace refs_detail {

/// rewrite, in v, the refs whose target kind (item_kind<T> through KindId)
/// is kind
template <template <typename> typename KindId, typename V>
void remap(V& v, std::size_t kind, const ref_remap& r) {
  using v_t = std::remove_cv_t<V>;
  if constexpr (!holds_refs_v<v_t>) {
    return;
  } else if constexpr (index_like<v_t>) {
    if (KindId<typename v_t::type>::get() == kind && v.value() != ref_remap::invalid) {
      v.value() = static_cast<typename v_t::value_t>(r(v.value()));
    }
  } else if constexpr (is_variant<v_t>::value) {
    std::visit([&](auto& alt) { remap<KindId>(alt, kind, r); }, v);
  } else if constexpr (is_tuple<v_t>::value) {
    std::apply([&](auto&... e) { (remap<KindId>(e, kind, r), ...); }, v);
  } else if constexpr (is_sparse_set_v<v_t>) {
    for (auto& e : v.dense_) remap<KindId>(e, kind, r);
  } else if constexpr (eigen_like<v_t>) {
    auto* p = v.data();
    for (decltype(v.size()) k = 0; k < v.size(); ++k) remap<KindId>(p[k], kind, r);
  } else if constexpr (map_like<v_t>) {
    static_assert(has_no_key_refs<v_t>::value,
                  "storage::remap_refs: refs in the keys of a map are not supported");
    for (auto& kv : v) remap<KindId>(kv.second, kind, r);
  } else if constexpr (range_like<v_t>) {
    for (auto& e : v) remap<KindId>(e, kind, r);
  }
}

/// the distinct kinds of the targets of a type_list, through KindId
template <template <typename> typename KindId, typename... Ts>
std::vector<std::size_t> kinds_of(type_list<Ts...>) {
  std::vector<std::size_t> kinds;
  (
      [&kinds] {
        const std::size_t k = KindId<Ts>::get();
        for (auto x : kinds) {
          if (x == k) return;
        }
        kinds.push_back(k);
      }(),
      ...);
  return kinds;
}

}  // namespace refs_detail
}  // namespace siconos::vkernel::storage
