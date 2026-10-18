#pragma once

#include <algorithm>
#include <any>
#include <atomic>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/sparse_set.hpp"
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
  /// returned by push_default for columns that do not number the items
  /// (sparse, dynamic properties)
  static constexpr std::size_t no_index = static_cast<std::size_t>(-1);

  virtual ~column_base() = default;
  virtual std::unique_ptr<column_base> clone() const = 0;

  /// storage::add: append a default value for a new item, in every kept step;
  /// returns the new item's index, or no_index
  virtual std::size_t push_default() = 0;

  /// storage::remove: item i goes away and the last item, last, takes its
  /// index (swap and pop), in every kept step
  virtual void move_back(std::size_t i, std::size_t last) = 0;

  /// number of items stored, or no_index for columns that do not number them
  virtual std::size_t item_count() const = 0;
};

namespace registry_detail {
template <typename T>
constexpr bool is_sparse_set = false;
template <typename K, typename V>
constexpr bool is_sparse_set<storage::sparse_set<K, V>> = true;
}  // namespace registry_detail

/**
 * @brief A typed column: the memory of one attribute (one storage per kept
 * step), or the sparse set of a dynamic property.
 *
 * push_default and move_back are the per-type code of storage::add and
 * storage::remove, so that both can walk the columns of an item at run time.
 */
template <typename Memory>
struct column final : column_base {
  Memory value{};

  column() = default;
  // runtime memory depth (keep<>): one storage per kept time step
  explicit column(std::size_t depth) {
    if constexpr (!registry_detail::is_sparse_set<Memory>) value.resize(depth);
  }

  std::unique_ptr<column_base> clone() const override {
    return std::make_unique<column>(*this);
  }

  std::size_t push_default() override {
    if constexpr (registry_detail::is_sparse_set<Memory>) {
      return no_index;  // absent means empty
    } else {
      // as storage::add does for the static database
      std::size_t index = 0;
      for (auto& storage : value) {
        using storage_t = std::decay_t<decltype(storage)>;
        if constexpr (match::push_back<storage_t>) {
          storage.push_back(typename storage_t::value_type{});
          index = std::size(storage) - 1;
        } else {
          // std::array or similar: a single item
          storage[0] = typename storage_t::value_type{};
          index = 0;
        }
      }
      return index;
    }
  }

  std::size_t item_count() const override {
    if constexpr (registry_detail::is_sparse_set<Memory>) {
      return no_index;
    } else {
      using storage_t = typename Memory::value_type;
      if constexpr (match::push_back<storage_t>) {
        return value.empty() ? 0 : std::size(value[0]);
      } else {
        return 1;
      }
    }
  }

  void move_back(std::size_t i, std::size_t last) override {
    if constexpr (registry_detail::is_sparse_set<Memory>) {
      // the removed item's entry goes, the last item's entry is re-keyed to i
      using key_t = typename Memory::key_type;
      value.erase(static_cast<key_t>(i));
      if (last != i && value.contains(static_cast<key_t>(last))) {
        auto moved = std::move(value.at(static_cast<key_t>(last)));
        value.erase(static_cast<key_t>(last));
        value.emplace(static_cast<key_t>(i), std::move(moved));
      }
    } else {
      // as storage::remove does for the static database (move_back)
      for (auto& storage : value) {
        using storage_t = std::decay_t<decltype(storage)>;
        if constexpr (match::push_back<storage_t>) {
          storage[i] = std::move(storage.back());
          storage.pop_back();
        }
      }
    }
  }
};

namespace registry_detail {
/// The slow path of registry::column: every column is registered when its
/// world is constructed (world::register_columns), so reaching it is a bug.
/// Not a template, so that access sites instantiate no creation code.
[[noreturn, gnu::noinline, gnu::cold]] inline void column_not_registered() {
  throw std::logic_error("storage::registry: column accessed before registration");
}

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

namespace registry_detail {
/// slot_id<Key>() + 1, read on the hot path without the guard of slot_id's
/// static: that guard is an atomic acquire load, which keeps the compiler
/// from hoisting or merging column lookups in inlined kernels.  Initialised
/// from slot_id when the program or library is loaded, so both agree; 0 means
/// "not initialised yet" (a lookup during static initialisation), which the
/// lookup turns into its checked path.
template <typename Key>
inline const std::size_t slot_plus_one = slot_id<Key>() + 1;
}  // namespace registry_detail

/**
 * @brief World-independent column storage.
 *
 * Not templated on the items: it is a vector of type-erased columns indexed by
 * slot_id<Key>, created lazily on first access.
 */
struct registry {
  std::vector<std::unique_ptr<column_base>> slots;

  /// for each item kind (indexed by slot_id<kind<Item>>), the slots of its
  /// columns, in registration order; slot ids, not pointers, so that a copy
  /// of the registry needs no remapping
  std::vector<std::vector<std::size_t>> kinds;
  std::vector<bool> kind_registered;

  template <typename Item>
  struct kind {};

  registry() = default;
  registry(registry&&) = default;
  registry& operator=(registry&&) = default;

  // deep copy, as the storage::make database it replaces is copyable
  registry(const registry& other)
      : slots(other.slots.size()), kinds(other.kinds), kind_registered(other.kind_registered) {
    for (std::size_t i = 0; i < slots.size(); ++i) {
      if (other.slots[i]) slots[i] = other.slots[i]->clone();
    }
  }

  registry& operator=(const registry& other) {
    if (this != &other) *this = registry(other);
    return *this;
  }

  // Hot path, small enough to be inlined at every access: the column exists.
  // Profiling bench-disks (900 disks) showed the former single function,
  // creation included, called out of line at each access: 3.7% of the run.
  template <typename Key, typename Memory>
  Memory& column() {
    // the slot id without the guard of slot_id's static (see slot_plus_one);
    // wraps to SIZE_MAX before slot_plus_one is initialised
    std::size_t id = registry_detail::slot_plus_one<Key> - 1;
#ifdef VKERNEL_UNCHECKED_LOOKUP
    // opt-in (CMake option VKERNEL_UNCHECKED_LOOKUP): no check at all; a
    // column that was never registered is undefined behaviour
    return static_cast<storage::column<Memory>*>(slots[id].get())->value;
#else
    if (id >= slots.size()) [[unlikely]] {
      id = slot_id<Key>();  // static initialisation order: ask slot_id itself
    }
    if (id < slots.size()) [[likely]] {
      if (column_base* slot = slots[id].get()) [[likely]] {
        return static_cast<storage::column<Memory>*>(slot)->value;
      }
    }
    registry_detail::column_not_registered();
#endif
  }

  /// Create the column of Key with its memory depth (keep<>), and list it
  /// among the columns of Item; done once per column and item, when the world
  /// registers its columns.
  template <typename Item, typename Key, typename Memory>
  Memory& create_column(std::size_t depth) {
    const std::size_t id = slot_id<Key>();
    if (id >= slots.size()) slots.resize(id + 1);
    auto& slot = slots[id];
    if (!slot) slot = std::make_unique<storage::column<Memory>>(depth);
    auto& columns = kind_columns<Item>();
    if (std::find(columns.begin(), columns.end(), id) == columns.end()) columns.push_back(id);
    return static_cast<storage::column<Memory>&>(*slot).value;
  }

  /// the column list of Item, marked as registered (possibly empty)
  template <typename Item>
  std::vector<std::size_t>& kind_columns() {
    const std::size_t k = slot_id<kind<Item>>();
    if (k >= kinds.size()) {
      kinds.resize(k + 1);
      kind_registered.resize(k + 1, false);
    }
    kind_registered[k] = true;
    return kinds[k];
  }

  /// storage::add: a default value in every column of Item; the new index
  template <typename Item>
  std::size_t push_item() {
    std::size_t index = 0;
    for (std::size_t id : registered_kind<Item>()) {
      if (std::size_t i = slots[id]->push_default(); i != column_base::no_index) index = i;
    }
    return index;
  }

  /// storage::remove: item i of Item goes away, the last one takes index i
  template <typename Item>
  void remove_item(std::size_t i) {
    const auto& columns = registered_kind<Item>();
    std::size_t count = column_base::no_index;
    for (std::size_t id : columns) {
      if (std::size_t c = slots[id]->item_count(); c != column_base::no_index) {
        count = c;
        break;
      }
    }
    if (count == column_base::no_index || count == 0) return;  // nothing stored
    for (std::size_t id : columns) slots[id]->move_back(i, count - 1);
  }

 private:
  template <typename Item>
  const std::vector<std::size_t>& registered_kind() const {
    const std::size_t k = slot_id<kind<Item>>();
    if (k >= kind_registered.size() || !kind_registered[k]) [[unlikely]] {
      throw std::logic_error("storage::registry: item kind not registered");
    }
    return kinds[k];
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

  /// the static attached storages of the world (dynamic_storage excluded),
  /// filtered once for world::register_columns
  using static_attached_t =
      decltype(mp::filter(all_properties_t{}, mp::is_a_model<[]<typename T>() {
                            return std::derived_from<T, some::attached_storage> &&
                                   !requires { typename T::dynamic_storage_t; };
                          }>));

  /// the dynamic properties of the world (dynamic_storage<attached<...>>)
  using dynamic_attached_t =
      decltype(mp::filter(all_properties_t{}, mp::is_a_model<[]<typename T>() {
                            return requires { typename T::dynamic_storage_t; };
                          }>));
};

namespace registry_detail {
template <typename T>
concept configuration_item = match::item<T> && requires { typename T::properties; } &&
                             std::derived_from<T, with_properties_base>;
}  // namespace registry_detail

/**
 * @brief What handles of a registry world know about it: the per-item
 * environment and the configuration's own properties, without the item
 * closure (step (a) P4).
 *
 * all_properties_t holds ONLY the configuration's properties (the
 * with_properties<...> given to make_world): the user refinements, wrappers
 * and configuration-declared attached storages.  Properties declared by items
 * (topology's attached storages, keep<> depths) are not there on purpose: a
 * lookup of those must name its declarer, and a world-wide one fails to
 * compile ("attached storage not found") instead of silently differing.
 */
template <template <typename> typename Env, typename ConfigProperties>
struct world_config {
  template <typename Item>
  using env = Env<Item>;

  using all_properties_t = ConfigProperties;
};

/// the with_properties<...> lists among the items given to make_world
template <match::item... Items>
using configuration_properties_t = decltype(mp::concat_all([] {
  if constexpr (registry_detail::configuration_item<Items>) {
    return typename Items::properties{};
  } else {
    return mp::tuple<>{};
  }
}()...));

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
  // the item of the key: searched in the item closure when the Info has one
  // (world_info, at registration), taken from the key otherwise (world_config,
  // at access sites; paired<> and attached<> both name their item)
  using item_t = decltype([] {
    if constexpr (requires { typename Info::all_items_t; }) {
      return item_attribute<A>(typename Info::all_items_t{});
    } else {
      return typename A::item{};
    }
  }());
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

  // memory level: keep<> depth.  The depth is a runtime size, not part of
  // the type (unlike memory_t's std::array), so that the type of a column
  // does not depend on which item declares a keep<> on it.  memory(step, m)
  // indexes m[step % size] either way.  Only meaningful with a world_info
  // (registration): keep<> is declared by items.
  using keeps_t = decltype(mp::filter(all_properties_t{}, mp::derive_from<property::keep>));
  static constexpr std::size_t depth = memory_size<A, keeps_t>();
  using type = std::vector<storage_t>;
};

/**
 * @brief Column type of a dynamic property P (dynamic_storage<attached<...>>):
 * a sparse set from the item index to the value, as prop_dynamic_memory uses
 * for the storage::make database (there, in a map keyed by name).
 */
template <typename Info, typename P>
struct dynamic_column_of {
  using env_t = typename Info::template env<typename P::item>;
  using value_t = typename traits::config<env_t>::template convert<typename P::type>::type;
  using type = storage::sparse_set<typename env_t::indice, value_t>;
};

/// P is a dynamic property attached to Item
template <typename Item>
static constexpr auto is_attached_dynamic = mp::is_a_model<[]<typename T>() constexpr {
  return match::attached_storage<T, Item> && requires { typename T::dynamic_storage_t; };
}>;

/**
 * @brief A world: a registry plus the compile-time Info to type its columns.
 *
 * It is its own store, so every existing access site, which goes through
 * mp::get<A>(data.store()), works unchanged. Handles carry world<Info>
 * instead of the database type of storage::make.  Info is a world_config
 * (make_world): the environment and the configuration's properties, not the
 * item closure, which only the constructor of make_world uses.
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
    } else if constexpr (requires { typename K::dynamic_storage_t; }) {
      return _registry.template column<K, typename dynamic_column_of<Info, K>::type>();
    } else {
      return _registry.template column<K, typename column_of<Info, K>::type>();
    }
  }

  /// Create every column of the world up front, with its depth: the
  /// attributes of every item and the static attached storages declared for
  /// it, as storage::make enumerates them.  Afterwards every access takes the
  /// fast path of registry::column.
  ///
  /// WorldInfo is the full world_info (item closure, all properties): only
  /// here, where the world is built, are the keep<> depths and the item lists
  /// needed.  The column types are those the access sites compute from Info
  /// alone, checked against the full computation.
  template <typename WorldInfo>
  void register_columns() {
    using static_attached_t = typename WorldInfo::static_attached_t;
    using dynamic_all_t = typename WorldInfo::dynamic_attached_t;
    mp::for_each(typename WorldInfo::all_items_t{}, [this]<typename Item>(Item) {
      // the static attached storages of Item
      using attached_t = decltype(mp::filter(static_attached_t{}, is_attached_storage<Item>));
      mp::for_each(concat(attributes(Item{}), attached_t{}), [this]<typename A>(A) {
        using access_t = typename column_of<Info, A>::type;
        static_assert(std::is_same_v<access_t, typename column_of<WorldInfo, A>::type>,
                      "the column type of an access site differs from the world's");
        (void)_registry.template create_column<Item, A, access_t>(
            column_of<WorldInfo, A>::depth);
      });
      // its dynamic properties: sparse columns, kept consistent by remove
      using dynamic_t = decltype(mp::filter(dynamic_all_t{}, is_attached_dynamic<Item>));
      mp::for_each(dynamic_t{}, [this]<typename P>(P) {
        using access_t = typename dynamic_column_of<Info, P>::type;
        static_assert(std::is_same_v<access_t, typename dynamic_column_of<WorldInfo, P>::type>,
                      "the dynamic column type of an access site differs from the world's");
        (void)_registry.template create_column<Item, P, access_t>(0);
      });
      (void)_registry.template kind_columns<Item>();  // registered, even without columns
    });
  }
};

/// storage::make counterpart that builds a world instead of a database; its
/// columns are registered on construction.  Handles of it carry
/// world<world_config<...>>: the configuration, not the item closure.
template <template <typename> typename Env, match::item... Items>
struct make_world : world<world_config<Env, configuration_properties_t<Items...>>> {
  /// the full description (item closure, every property): for registration
  /// and for assembly code such as the python bindings, never in handles
  using world_info_t = world_info<Env, Items...>;

  make_world() { this->template register_columns<world_info_t>(); }
};

}  // namespace siconos::storage
