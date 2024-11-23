#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/traits/traits.hpp"
#include "siconos/storage/memory.hpp"

namespace siconos::storage {

namespace match = siconos::storage::pattern::match;

template <match::item Item, typename Wrappers, typename Storage>
static constexpr auto apply_wrapper(Storage storage)
{
  auto tpl = ground::filter(Wrappers{},
                            ground::is_a_model<[]<typename T>() consteval {
                              return std::is_same_v<Item, typename T::type>;
                            }>);
  if constexpr (ground::size(tpl) > ground::size_c<0>) {
    return typename std::decay_t<decltype(tpl[0_c])>::template wrapper<
        Storage>{};
  }
  else {
    // without wrapper
    return Storage{};
  }
};

template <typename Env, match::attribute... Attributes>
struct unit_storage {
  using env = Env;

  using type = ground::pre_map<ground::key_value<
      Attributes,
      typename traits::config<env>::template convert<Attributes>::type>...>;
};

template <typename Info, typename M>
using with_info_t =
    decltype(ground::prepend(M{}, ground::key_value<info, Info>{}));

template <typename Env, match::item... Items>
struct item_storage {
  struct iinfo {
    using env = Env;
    using items = gather<Items...>;
    using all_items_t = decltype(ground::fold_left(
        // take last defined items
        ground::reverse(ground::concat_all(all_items(Items{})...)),
        ground::make_tuple(),
        []<typename Acc, typename Elem>(Acc acc, Elem elem) {
          if constexpr (ground::contains(acc, elem)) {
            return acc;
          }
          else {
            return ground::append(acc, elem);
          };
        }));
    using all_attributes_t = decltype(ground::flatten(
        ground::concat_all(all_attributes(Items{})...)));
    using all_properties_t = decltype(ground::flatten(
        ground::concat_all(all_properties(Items{})...)));

    // subset of all_properties_t
    using all_attached_storages_t = decltype(ground::filter(
        all_properties_t{}, ground::derive_from<some::attached_storage>));
  };

  constexpr static auto a =
      ground::concat(typename iinfo::all_attributes_t{},
                     typename iinfo::all_attached_storages_t{});
  // base map with declared attributes
  using map_t = decltype(ground::unpack(
      ground::concat(typename iinfo::all_attributes_t{},
                     typename iinfo::all_attached_storages_t{}),
      []<typename... Attributes>(Attributes...) {
        return
            typename unit_storage<typename iinfo::env, Attributes...>::type{};
      }));

  // final map with info added
  using type = with_info_t<iinfo, map_t>;
};

template <typename Info>
static constexpr auto item_storage_transform =
    []<typename D>(D&& d, auto&& f) constexpr -> decltype(auto) {
  using info_t = Info;

  return ground::transform(d, [&f]<typename P>(P&& key_value) {
    auto&& key = ground::first(key_value);
    auto&& value = ground::second(key_value);
    using key_t = std::decay_t<decltype(key)>;
    using value_t = std::decay_t<decltype(value)>;

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

static constexpr auto attribute_storage_transform =
    []<typename D, typename F>(D&& d, F&& f) constexpr -> decltype(auto) {
  return ground::pre_map_value_transform(
      std::forward<D>(d), [&f]<typename K, typename S>(K, S&& s) {
        if constexpr (match::attribute<typename K::type>) {
          return std::forward<F>(f)(typename K::type{}, std::forward<S>(s));
        }
        else {
          return std::move(s);
        }
      });
};

template <typename Env, match::item... Items>
struct make {
  static constexpr decltype(auto) internal_build()
  {
    using item_storage_t = item_storage<Env, Items...>;
    using info_t = typename item_storage_t::iinfo;
    auto base_storage = typename item_storage_t::type{};
    return ground::to_database(attribute_storage_transform(
        item_storage_transform<info_t>(
            attribute_storage_transform(
                base_storage,
                // attribute level for base storage specifications
                []<match::attribute Attribute, typename Storage>(
                    Attribute, Storage& s) -> decltype(auto) {
                  // if attribute is derived from one of diagonal
                  // specifications, etc.
                  return typename traits::config<typename info_t::env>::
                      template convert<decltype(refine_recursively_attribute(
                          base_storage, Attribute{}))>::type{};
                }),
            // item level: collection depends on item property
            []<match::item Item, match::attribute Attr>(Item item, Attr attr,
                                                        auto s) {
              using storage_t = std::decay_t<decltype(s)>;

              if constexpr (match::wrap<Item>) {
                return ground::key_value<
                    Attr,
                    typename traits::config<Env>::template convert<
                        typename Item::template wrapper<storage_t>>::type>{};
              }
              else {
                // look for wrap specified in properties
                using all_wrappers_t =
                    decltype(pre_map_all_properties_as<property::wrapped>(
                        base_storage));

                using storage_wrapped =
                    decltype(apply_wrapper<Item, all_wrappers_t>(s));

                if constexpr (!std::is_same_v<storage_wrapped, storage_t>) {
                  using storage_wrapped_and_converted =
                      typename traits::config<typename info_t::env>::
                          template convert<storage_wrapped>::type;

                  return ground::key_value<Attr,
                                           storage_wrapped_and_converted>{};
                }
                else {
                  return ground::key_value<
                      Attr,
                      typename Env::template default_storage<storage_t>>{};
                }
              }
            }),
        // attribute level: memory depends on keeps
        []<match::attribute Attribute>(Attribute attr, auto s) {
          using storage_t = std::decay_t<decltype(s)>;

          using all_keeps_t =
              decltype(pre_map_all_properties_as<property::keep>(
                  base_storage));
          return memory_t<storage_t, memory_size<Attribute, all_keeps_t>()>{};
        }));
  };

  struct store_t : decltype(internal_build()) {};

  store_t _store;

  make() : _store(store_t()){};

  store_t& store() { return _store; }
};
}  // namespace siconos::storage
