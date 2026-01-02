#pragma once

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
template <match::item Item, typename Wrappers, typename Storage>
static constexpr auto apply_wrapper(Storage storage)
{
  // Filter wrappers to find those matching the item type
  auto tpl =
      mp::filter(Wrappers{}, mp::is_a_model<[]<typename T>() consteval {
                   return std::is_same_v<Item, typename T::type>;
                 }>);
  // If matching wrapper found, apply it; otherwise return storage as-is
  if constexpr (mp::size(tpl) > mp::size_c<0>) {
    return typename std::decay_t<decltype(tpl[0_c])>::template wrapper<
        Storage>{};
  }
  else {
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
      Attributes,
      typename traits::config<env>::template convert<Attributes>::type>...>;
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
template <typename Env, match::item... Items>
struct item_storage {
  // Internal info structure containing computed type aliases
  struct iinfo {
    using env = Env;
    using items = gather<Items...>;

    // Compute all unique items through a fold-left operation
    using all_items_t = decltype(mp::fold_left(
        // Process items in reverse order to prioritize later definitions
        mp::reverse(mp::concat_all(all_items(Items{})...)), mp::make_tuple(),
        []<typename Acc, typename Elem>(Acc acc, Elem elem) {
          // Skip duplicate items, keeping only the first occurrence
          if constexpr (mp::contains(acc, elem)) {
            return acc;
          }
          else {
            return mp::append(acc, elem);
          };
        }));

    // Flatten all attributes from all items
    using all_attributes_t =
        decltype(mp::flatten(mp::concat_all(all_attributes(Items{})...)));

    // Flatten all properties from all items
    using all_properties_t =
        decltype(mp::flatten(mp::concat_all(all_properties(Items{})...)));

    // Extract attached storage properties (subset of all properties)
    using all_attached_storages_t = decltype(mp::filter(
        all_properties_t{}, mp::derive_from<some::attached_storage>));
  };

  // Base storage map combining declared attributes and attached storages
  using map_t = decltype(mp::unpack(
      mp::concat(typename iinfo::all_attributes_t{},
                 typename iinfo::all_attached_storages_t{}),
      []<typename... Attributes>(Attributes...) {
        return
            typename unit_storage<typename iinfo::env, Attributes...>::type{};
      }));

  // Final storage type with info prepended
  using type = with_info_t<iinfo, map_t>;
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
      return f(item_attribute<attr_t>(typename info_t::all_items_t{}),
               attr_t{}, std::forward<value_t>(value));
    }
    else {
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
        }
        else {
          return std::move(s);
        }
      });
};

/**
 * @brief Factory struct for building complete storage structures
 *
 * @tparam Env The environment type
 * @tparam Items Variadic list of item types to include
 */
template <typename Env, match::item... Items>
struct make {
  /**
   * @brief Internal build method that constructs the complete storage
   * @return The built storage structure
   */
  static constexpr decltype(auto) internal_build()
  {
    using item_storage_t = item_storage<Env, Items...>;
    using info_t = typename item_storage_t::iinfo;
    auto base_storage = typename item_storage_t::type{};

    // Build the storage through a multi-stage transformation pipeline
    return mp::to_database(attribute_storage_transform(
        item_storage_transform<info_t>(
            attribute_storage_transform(
                base_storage,
                // Attribute level: refine attributes based on specifications
                []<match::attribute Attribute, typename Storage>(
                    Attribute, Storage& s) -> decltype(auto) {
                  // Convert attributes using recursive refinement
                  return typename traits::config<typename info_t::env>::
                      template convert<decltype(refine_recursively_attribute(
                          base_storage, Attribute{}))>::type{};
                }),
            // Item level: handle collection wrappers for items
            []<match::item Item, match::attribute Attr>(Item item, Attr attr,
                                                        auto s) {
              using storage_t = std::decay_t<decltype(s)>;

              // Handle wrapped items by applying their wrapper
              if constexpr (match::wrap<Item>) {
                return mp::key_value<
                    Attr,
                    typename traits::config<Env>::template convert<
                        typename Item::template wrapper<storage_t>>::type>{};
              }
              else {
                // Look for wrappers specified in properties
                using all_wrappers_t =
                    decltype(pre_map_all_properties_as<property::wrapped>(
                        base_storage));

                using storage_wrapped =
                    decltype(apply_wrapper<Item, all_wrappers_t>(s));

                // Apply wrapper if different from original storage type
                if constexpr (!std::is_same_v<storage_wrapped, storage_t>) {
                  using storage_wrapped_and_converted =
                      typename traits::config<typename info_t::env>::
                          template convert<storage_wrapped>::type;

                  return mp::key_value<Attr, storage_wrapped_and_converted>{};
                }
                else {
                  // Use environment's default storage if no wrapper
                  return mp::key_value<
                      Attr,
                      typename Env::template default_storage<storage_t>>{};
                }
              }
            }),
        // Memory level: determine memory size based on keeps properties
        []<match::attribute Attribute>(Attribute attr, auto s) {
          using storage_t = std::decay_t<decltype(s)>;

          using all_keeps_t =
              decltype(pre_map_all_properties_as<property::keep>(
                  base_storage));
          return memory_t<storage_t, memory_size<Attribute, all_keeps_t>()>{};
        }));
  };

  // The actual storage member variable
  struct store_t : decltype(internal_build()){};

  store_t _store;

  // Accessor for the storage
  store_t& store() { return _store; }
};
}  // namespace siconos::storage
