#pragma once

#include <any>
#include <cstddef>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/refs.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::storage {

namespace remap_detail {
/// the number of columns of the database<Pairs...> base of a storage::make
/// store
template <typename... Pairs>
std::integral_constant<std::size_t, sizeof...(Pairs)> pair_count(
    const mp::database<Pairs...>*);
}  // namespace remap_detail

/**
 * The rows of the items of Item moved as r says (see ref_remap): rewrite the
 * refs to them held in the store, in every column and every kept step.
 * storage::remove, storage::remove_many and storage::permute call it.
 *
 * Registry world: the columns registered as referring to Item. storage::make
 * database: every column of the store, and the existing sparse sets of the
 * dynamic properties.
 */
template <typename Item>
void remap_refs(auto& data, const ref_remap& r) {
  const std::size_t kind = kind_id<Item>::get();
  if constexpr (requires { data.store()._registry; }) {
    data.store()._registry.template remap_refs<Item>(r);
  } else {
    auto& db = data.store();
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      (
          [&] {
            auto& value = boost::hana::second(boost::hana::at_c<I>(db.store));
            if constexpr (refers_to_v<std::decay_t<decltype(value)>, Item>) {
              refs_detail::remap<kind_id>(value, kind, r);
            }
          }(),
          ...);
    }(std::make_index_sequence<decltype(remap_detail::pair_count(&db))::value>{});

    using info_t = get_info_t<std::decay_t<decltype(data)>>;
    mp::for_each(typename info_t::all_properties_t{}, [&]<typename P>(P) {
      if constexpr (requires { typename P::dynamic_storage_t; }) {
        using env_t = typename info_t::template env<typename P::item>;
        using indice = typename env_t::indice;
        using value_t =
            typename traits::config<env_t>::template convert<typename P::type>::type;
        if constexpr (refers_to_v<value_t, Item>) {
          // as prop_dynamic_map, without creating the set
          auto& map = db._dynamic_properties;
          if (auto it = map.find(std::string(typeid(P).name())); it != map.end()) {
            refs_detail::remap<kind_id>(
                std::any_cast<sparse_set<indice, value_t>&>(it->second), kind, r);
          }
        }
      }
    });
  }
}

}  // namespace siconos::storage
