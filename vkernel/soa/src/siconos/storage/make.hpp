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

  // unique: a property declared for several item kinds that resolve to the
  // same type (simul::topology's dfinteraction and ddinteraction are both
  // empty_item without dynamic dof systems) would otherwise be a duplicate
  // column
  using all_properties_t = decltype(mp::tuple_unique(
      mp::flatten(mp::concat_all(all_properties(std::declval<Items>())...))));

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
    // no column attached to empty_item: it stands for an absent item kind
    // (a topology without dynamic dof systems), no item of it is ever added
    using attached_storages_t = decltype([] {
      if constexpr (match::empty_item<Item>) {
        return mp::tuple<>{};
      } else {
        return mp::filter(typename iinfo::all_properties_t{}, mp::is_a_model<[]<typename T>() {
                            return match::attached_storage<T, Item> &&
                                   !requires { typename T::dynamic_storage_t; };
                          }>);
      }
    }());
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

namespace build_detail {

/**
 * @brief Storage of one attribute column of Item, before its memory depth.
 *
 * The three stages shared by both databases: attribute refinement, item
 * wrapper or default storage, and the keep<> depth. column_type stores the
 * depth in the type (memory_t), registry::column_of at run time.
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

  // memory level: keep<> depth
  using keeps_t = decltype(mp::filter(all_properties_t{}, mp::derive_from<property::keep>));
  static constexpr std::size_t depth = memory_size<A, keeps_t>();
};

/**
 * @brief Column type of one attribute key of the compile-time database.
 */
template <typename Info, match::attribute A>
struct column_type
    : column_storage<Info, A, decltype(item_attribute<A>(typename Info::all_items_t{}))> {
  using type = memory_t<typename column_type::storage_t, column_type::depth>;
};

// non-attribute pair (info, attached storages): kept as-is
template <typename Info, typename Pair>
struct pair_type {
  using type = Pair;
};

// attribute pair: its column type (a class partial specialization takes its
// requires clause after the template parameter list only)
template <typename Info, typename A, typename B>
  requires(match::attribute<typename A::type>)
struct pair_type<Info, boost::hana::pair<A, B>> {
  using type =
      mp::key_value<typename A::type, typename column_type<Info, typename A::type>::type>;
};

/**
 * @brief Single-pass computation of the compile-time database columns.
 *
 * Every column goes through column_type in one parameter-pack pass over the
 * base storage column keys. The previous pipeline applied its stages as
 * three whole-tuple hana::transform passes, each rebuilding the complete
 * tuple of key-value pairs.
 */
template <typename Info, typename Pairs>
struct pack_columns;

template <typename Info, typename... Ps>
struct pack_columns<Info, boost::hana::tuple<Ps...>> {
  using type = boost::hana::tuple<typename pair_type<Info, Ps>::type...>;

  static constexpr type build() { return type(typename pair_type<Info, Ps>::type{}...); }
};

}  // namespace build_detail

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
    using base_storage_t = typename item_storage_t::type;
    return mp::to_database(build_detail::pack_columns<info_t, base_storage_t>::build());
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
