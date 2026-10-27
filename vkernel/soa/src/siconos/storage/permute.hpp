#pragma once

#include <cassert>
#include <utility>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

namespace permute_detail {
/// a storage of one step, reordered: its item k is the item perm[k]
static auto gather = [](const std::vector<std::size_t>& perm, auto& a) {
  if constexpr (match::push_back<std::decay_t<decltype(a)>>) {
    assert(std::size(a) == std::size(perm));
    std::decay_t<decltype(a)> out;
    out.reserve(std::size(perm));
    for (std::size_t k = 0; k < std::size(perm); ++k) {
      out.push_back(std::move(a[perm[k]]));
    }
    a = std::move(out);
  }
  // else a single item (std::array): nothing to reorder
};

/// a sparse set re-keyed: the entry of the item perm[k] is the entry of k
static auto rekey = [](const std::vector<std::size_t>& inverse, auto& set) {
  using set_t = std::decay_t<decltype(set)>;
  using key_t = typename set_t::key_type;
  set_t out;
  out.reserve(set.size());
  for (std::size_t j = 0; j < set.dense_.size(); ++j) {
    out.emplace(static_cast<key_t>(inverse[static_cast<std::size_t>(set.dense_keys_[j])]),
                std::move(set.dense_[j]));
  }
  set = std::move(out);
};
}  // namespace permute_detail

/** Reorder the items of type Item: the item of index k becomes the item of
    index perm[k] before the call, in every attribute, every attached storage
    and every dynamic property (re-keyed), in every kept step.

    perm is a permutation of [0, number of items). The handles and the refs
    on the items of Item held elsewhere are not updated: no item may hold a
    ref on an Item across the call. */
template <typename Item>
void permute(auto& data, const std::vector<std::size_t>& perm) {
  std::vector<std::size_t> inverse(std::size(perm));
  for (std::size_t k = 0; k < std::size(perm); ++k) {
    inverse[perm[k]] = k;
  }

  if constexpr (requires { data.store()._registry; }) {
    data.store()._registry.template permute_item<Item>(perm, inverse);
  } else {
    using info_t = get_info_t<decltype(data)>;
    using all_keeps_t = decltype(all_properties_as<property::keep>(data));
    using indice = typename info_t::template env<Item>::indice;
    using env_t = typename info_t::template env<Item>;

    using all_attrs_t = decltype(mp::tuple_unique(
        concat(attributes(Item{}), attached_storages(Item{}, data))));

    auto attrs = mp::filter(all_attrs_t{}, mp::is_a_model<[]<typename T>() constexpr {
                              return !requires { typename T::dynamic_storage_t; };
                            }>);
    using dynamic_attrs_t =
        decltype(mp::filter(all_attrs_t{}, mp::is_a_model<[]<typename T>() constexpr {
                              return requires { typename T::dynamic_storage_t; };
                            }>));

    mp::for_each(dynamic_attrs_t{}, [&data, &inverse]<typename P>(P) {
      permute_detail::rekey(inverse, prop_dynamic_map<P, env_t, indice>(data));
    });

    mp::for_each(attrs, [&data, &perm]<match::attribute A>(A) {
      return mp::for_each(
          mp::range<memory_size<A, all_keeps_t>()>, [&data, &perm](indice step) {
            permute_detail::gather(perm, memory(step, mp::get<A>(data.store())));
          });
    });
  }
}

}  // namespace siconos::storage
