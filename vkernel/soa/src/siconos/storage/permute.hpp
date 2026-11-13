#pragma once

#include <cassert>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/remap_refs.hpp"

namespace siconos::storage {

namespace permute_detail {
/// a storage of one step, reordered in place: its item k is the item
/// perm[k]. Its buffer is kept: a pointer to it (the coordinates given to
/// the neighbor search) stays valid
static auto gather = [](const std::vector<std::size_t>& perm, auto& a) {
  if constexpr (match::push_back<std::decay_t<decltype(a)>>) {
    assert(std::size(a) == std::size(perm));
    std::decay_t<decltype(a)> tmp;
    tmp.reserve(std::size(perm));
    for (std::size_t k = 0; k < std::size(perm); ++k) {
      tmp.push_back(std::move(a[perm[k]]));
    }
    for (std::size_t k = 0; k < std::size(perm); ++k) {
      a[k] = std::move(tmp[k]);
    }
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

    perm is a permutation of [0, number of items). The refs to the items of
    Item held in the store follow them (see remap_refs); the handles and the
    refs held outside the store do not. */
template <typename Item>
void permute(auto& data, const std::vector<std::size_t>& perm) {
  std::vector<std::size_t> inverse(std::size(perm));
  for (std::size_t k = 0; k < std::size(perm); ++k) {
    if (perm[k] >= std::size(perm)) [[unlikely]] {
      throw std::invalid_argument("storage::permute: index out of the permutation");
    }
    inverse[perm[k]] = k;
  }

  if constexpr (requires { data.store()._registry; }) {
    data.store()._registry.template permute_item<Item>(perm, inverse);
  } else {
    using info_t = get_info_t<decltype(data)>;
    using all_time_levels_t = decltype(all_properties_as<property::time_levels>(data));
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

    if constexpr (mp::size(attrs) > mp::size_c<0>) {
      using first_t = std::decay_t<decltype(attrs[0_c])>;
      const auto& first_storage = memory(0, mp::get<first_t>(data.store()));
      if constexpr (match::push_back<std::decay_t<decltype(first_storage)>>) {
        if (std::size(first_storage) != std::size(perm)) [[unlikely]] {
          throw std::invalid_argument("storage::permute: " + std::to_string(std::size(perm)) +
                                      " indices for " +
                                      std::to_string(std::size(first_storage)) + " items");
        }
      }
    }

    mp::for_each(dynamic_attrs_t{}, [&data, &inverse]<typename P>(P) {
      permute_detail::rekey(inverse, prop_dynamic_map<P, env_t, indice>(data));
    });

    mp::for_each(attrs, [&data, &perm]<match::attribute A>(A) {
      return mp::for_each(
          mp::range<memory_size<A, all_time_levels_t>()>, [&data, &perm](indice step) {
            permute_detail::gather(perm, memory(step, mp::get<A>(data.store())));
          });
    });

    // the refs to the moved items follow them (the registry does it in
    // permute_item)
    remap_refs<Item>(data, ref_remap::table(inverse));
  }
}

}  // namespace siconos::storage
