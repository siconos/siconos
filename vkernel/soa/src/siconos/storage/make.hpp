#pragma once

#include <any>

#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
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
      return typename __type_pack_element<of::index, Ws...>::template wrapper<Storage>{};
    }(std::decay_t<Wrappers>{});
  } else {
    // without wrapper
    return Storage{};
  }
};

/**
 * @brief Defines storage unit for a set of attributes in a given environment
 *
 * @tparam Env The environment type
 * @tparam Attributes Variadic list of attribute types
 */
template <typename Env, match::attribute... Attributes>
struct unit_storage {
  using env = Env;

  using type = mp::pre_map<mp::key_value<
      Attributes, typename traits::config<env>::template convert<Attributes>::type>...>;
};

/**
 * @brief Type alias to prepend info key-value pair to a map
 *
 * @tparam Info The info type to prepend
 * @tparam M The map type
 */
template <typename Info, typename M>
using with_info_t = decltype(mp::prepend(M{}, mp::key_value<info, Info>{}));

/**
 * @brief Main storage structure for managing items and their metadata
 *
 * @tparam Env The environment type
 * @tparam Items Variadic list of item types
 */
template <template <typename> typename EnvTemplate, match::item... Items>
struct item_storage_info {
  using items = gather<Items...>;

  using all_items_t =
      decltype(mp::tuple_unique(mp::concat_all(all_items(std::declval<Items>())...)));

  using all_attributes_t =
      decltype(mp::flatten(mp::concat_all(all_attributes(std::declval<Items>())...)));

  using all_properties_t =
      decltype(mp::flatten(mp::concat_all(all_properties(std::declval<Items>())...)));

  using all_attached_storages_t =
      decltype(mp::filter(all_properties_t{}, mp::derive_from<some::attached_storage>));
};

template <template <typename> typename EnvTemplate, match::item... Items>
struct item_storage {
  struct iinfo : item_storage_info<EnvTemplate, Items...> {
    template <typename Item>
    using env = EnvTemplate<Item>;
  };

  template <match::item Item>
  struct unit_type_for;

  template <match::item Item>
  struct unit_type_for {
    using attributes_t = decltype(attributes(Item{}));
    using attached_storages_t = decltype(mp::filter(
        typename iinfo::all_properties_t{}, mp::is_a_model<[]<typename T>() {
          return match::attached_storage<T, Item> &&
                 !requires { typename T::dynamic_storage_t; };
        }>));
    using all_attrs_t = decltype(concat(attributes_t{}, attached_storages_t{}));
    using type = decltype(mp::unpack(all_attrs_t{}, []<typename... As>(As...) {
      return typename unit_storage<EnvTemplate<Item>, As...>::type{};
    }));
  };

  // final storage with info
  using type =
      with_info_t<iinfo,
                  decltype(mp::unpack(
                      typename iinfo::all_items_t{}, []<typename... AllItems>(AllItems...) {
                        return mp::concat_all(typename unit_type_for<AllItems>::type{}...);
                      }))>;
};

/**
 * @brief Transform function that processes storage items based on their
 * attributes
 *
 * @tparam Info The info type
 */
template <typename Info>
static constexpr auto item_storage_transform =
    []<typename D>(D&& d, auto&& f) constexpr -> decltype(auto) {
  using info_t = Info;

  // Apply transformation to each key-value pair in the storage
  return mp::transform(d, [&f]<typename P>(P&& key_value) {
    auto&& key = mp::first(key_value);
    auto&& value = mp::second(key_value);
    using key_t = std::decay_t<decltype(key)>;
    using value_t = std::decay_t<decltype(value)>;

    // Process attributes differently from other key-value pairs
    if constexpr (match::attribute<typename key_t::type>) {
      using attr_t = typename key_t::type;
      return f(item_attribute<attr_t>(typename info_t::all_items_t{}), attr_t{},
               std::forward<value_t>(value));
    } else {
      return std::forward<std::decay_t<P>>(key_value);
    }
  });
};

/**
 * @brief Transform function for attribute storage, applying a function to
 * each attribute
 *
 * @param d The storage data to transform
 * @param f The function to apply to each attribute
 * @return The transformed storage
 */
static constexpr auto attribute_storage_transform =
    []<typename D, typename F>(D&& d, F&& f) constexpr -> decltype(auto) {
  return mp::pre_map_value_transform(
      std::forward<D>(d), [&f]<typename K, typename S>(K, S&& s) {
        // Only apply function to attribute types
        if constexpr (match::attribute<typename K::type>) {
          return std::forward<F>(f)(typename K::type{}, std::forward<S>(s));
        } else {
          return std::move(s);
        }
      });
};

/**
 * @brief Factory struct for building a complete compile-time typed database
 *
 * The whole database type (every column of every item, refined, wrapped and
 * sized for its keep<> depth) is computed here, at compile time, and every
 * handle carries it: slowest compile, fastest runtime access.
 * For the counterpart where columns are compile-time typed per key too, but
 * only created (registered) at run time - much faster to compile, at the
 * cost of a small runtime column lookup - see storage::make_runtime_database.
 *
 * @tparam Env The environment type
 * @tparam Items Variadic list of item types to include
 */
template <template <typename> typename Env, match::item... Items>
struct make_compile_time_database {
  /**
   * @brief Internal build method that constructs the complete storage
   * @return The built storage structure
   */
  static constexpr decltype(auto) internal_build() {
    using item_storage_t = item_storage<Env, Items...>;
    using info_t = typename item_storage_t::iinfo;
    auto base_storage = typename item_storage_t::type{};
    using all_wrappers_t =
        decltype(pre_map_all_properties_as<property::wrapped>(base_storage));
    using all_keeps_t = decltype(pre_map_all_properties_as<property::keep>(base_storage));

    // Build the storage through a multi-stage transformation pipeline
    return mp::to_database(attribute_storage_transform(
        item_storage_transform<info_t>(
            attribute_storage_transform(
                base_storage,
                // Attribute level: refine attributes based on specifications
                []<match::attribute Attribute, typename Storage>(
                    Attribute, Storage& s) -> decltype(auto) {
                  // Convert attributes using recursive refinement
                  using item_t =
                      decltype(item_attribute<Attribute>(typename info_t::all_items_t{}));
                  using env_t = typename info_t::template env<item_t>;

                  return typename traits::config<env_t>::template convert<
                      decltype(refine_recursively_attribute(base_storage,
                                                            Attribute{}))>::type{};
                }),
            // Item level: handle collection wrappers for items
            []<match::item Item, match::attribute Attr>(Item item, Attr attr, auto s) {
              using storage_t = std::decay_t<decltype(s)>;
              using env_t = typename info_t::template env<Item>;

              // Handle wrapped items by applying their wrapper
              if constexpr (match::wrap<Item>) {
                return mp::key_value<Attr,
                                     typename traits::config<env_t>::template convert<
                                         typename Item::template wrapper<storage_t>>::type>{};
              } else {
                using storage_wrapped = decltype(apply_wrapper<Item, all_wrappers_t>(s));

                // Apply wrapper if different from original storage type
                if constexpr (!std::is_same_v<storage_wrapped, storage_t>) {
                  using storage_wrapped_and_converted =
                      typename traits::config<env_t>::template convert<storage_wrapped>::type;

                  return mp::key_value<Attr, storage_wrapped_and_converted>{};
                } else {
                  // Use environment's default storage if no wrapper
                  return mp::key_value<
                      Attr, typename Env<Item>::template default_storage<storage_t>>{};
                }
              }
            }),
        // Memory level: determine memory size based on keeps properties
        []<match::attribute Attribute>(Attribute attr, auto s) {
          using storage_t = std::decay_t<decltype(s)>;

          return memory_t<storage_t, memory_size<Attribute, all_keeps_t>()>{};
        }));
  };

  using database_t = decltype(internal_build());

  // The actual storage member variable
  struct store_t : database_t {
    std::unordered_map<std::string, std::any> _dynamic_properties;
  };

  store_t _store;

  // Accessor for the storage
  store_t& store() { return _store; }
};

/// previous names of make_compile_time_database, kept so existing call sites
/// compile
template <template <typename> typename Env, match::item... Items>
using make = make_compile_time_database<Env, Items...>;
template <template <typename> typename Env, match::item... Items>
using make_database = make_compile_time_database<Env, Items...>;

template <typename Make1, typename Make2>
struct make_merged {
  using db1_t = typename Make1::database_t;
  using db2_t = typename Make2::database_t;

  using merged_db_t =
      decltype(mp::merge_databases(std::declval<db1_t>(), std::declval<db2_t>()));

  struct store_t : merged_db_t {
    std::unordered_map<std::string, std::any> _dynamic_properties;
  };

  store_t _store;

  store_t& store() { return _store; }
};

}  // namespace siconos::storage
