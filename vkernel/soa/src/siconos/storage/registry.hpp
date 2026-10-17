#pragma once

#include <any>
#include <atomic>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/traits/traits.hpp"

namespace siconos::storage {

namespace match = siconos::storage::pattern::match;

/**
 * @brief Type-erased column slot.
 *
 * The registry owns columns without knowing their types; the typed view is
 * recovered at the access site, which computes the column type from the
 * attribute key alone (see column_of).
 */
struct column_base {
  virtual ~column_base() = default;
  virtual std::unique_ptr<column_base> clone() const = 0;
};

template <typename Memory>
struct column final : column_base {
  Memory value{};

  std::unique_ptr<column_base> clone() const override {
    return std::make_unique<column>(*this);
  }
};

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

/**
 * @brief World-independent column storage.
 *
 * Not templated on the items: it is a vector of type-erased columns indexed by
 * slot_id<Key>, created lazily on first access.
 */
struct registry {
  std::vector<std::unique_ptr<column_base>> slots;

  registry() = default;
  registry(registry&&) = default;
  registry& operator=(registry&&) = default;

  // deep copy, as the storage::make database it replaces is copyable
  registry(const registry& other) : slots(other.slots.size()) {
    for (std::size_t i = 0; i < slots.size(); ++i) {
      if (other.slots[i]) slots[i] = other.slots[i]->clone();
    }
  }

  registry& operator=(const registry& other) {
    if (this != &other) *this = registry(other);
    return *this;
  }

  template <typename Key, typename Memory>
  Memory& column() {
    const std::size_t id = slot_id<Key>();
    if (id >= slots.size()) slots.resize(id + 1);
    auto& slot = slots[id];
    if (!slot) slot = std::make_unique<storage::column<Memory>>();
    return static_cast<storage::column<Memory>&>(*slot).value;
  }
};

/**
 * @brief Compile-time description of a world, without its storage.
 *
 * Only what the access sites ask for through get_info_t: the item closure,
 * the gathered properties and the per-item environment. Unlike
 * item_storage_info it does not gather all attributes, and unlike
 * item_storage it builds no unit storage for any of them.
 */
template <template <typename> typename Env, match::item... Items>
struct world_info {
  using items = gather<Items...>;

  using all_items_t =
      decltype(mp::tuple_unique(mp::concat_all(all_items(std::declval<Items>())...)));

  using all_properties_t =
      decltype(mp::flatten(mp::concat_all(all_properties(std::declval<Items>())...)));

  template <typename Item>
  using env = Env<Item>;
};

/**
 * @brief Concrete column type of one attribute key.
 *
 * The same three stages as storage::make::internal_build (attribute
 * refinement, item wrapper or default storage, memory depth), applied to a
 * single key instead of to the whole database. It is a class template, so it
 * is computed once per key however many access sites ask for it.
 */
template <typename Info, match::attribute A>
struct column_of {
  using all_properties_t = typename Info::all_properties_t;
  using item_t = decltype(item_attribute<A>(typename Info::all_items_t{}));
  using env_t = typename Info::template env<item_t>;

  // attribute level: user refinements over the recursively rebuilt attribute
  using recursed_t = typename recursive_rebuild<A>::type;
  using refines_t = decltype(mp::filter(
      mp::filter(all_properties_t{}, mp::derive_from<property::refine>),
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

  // memory level: keep<> depth
  using keeps_t = decltype(mp::filter(all_properties_t{}, mp::derive_from<property::keep>));
  using type = memory_t<storage_t, memory_size<A, keeps_t>()>;
};

/**
 * @brief A world: a registry plus the compile-time Info to type its columns.
 *
 * It is its own store, so every existing access site, which goes through
 * mp::get<A>(data.store()), works unchanged. Handles carry world<Info>
 * instead of the database type of storage::make.
 */
template <typename Info>
struct world {
  registry _registry;
  std::unordered_map<std::string, std::any> _dynamic_properties;

  world& store() { return *this; }

  template <typename K>
  decltype(auto) column() {
    if constexpr (std::is_same_v<K, info>) {
      return Info{};
    } else {
      return _registry.template column<K, typename column_of<Info, K>::type>();
    }
  }
};

/// storage::make counterpart that builds a world instead of a database.
template <template <typename> typename Env, match::item... Items>
struct make_world : world<world_info<Env, Items...>> {};

}  // namespace siconos::storage
