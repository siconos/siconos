#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/refs.hpp"
#include "siconos/storage/remap_refs.hpp"

namespace siconos::storage {

static auto move_back = [](const auto i, auto& a) constexpr {
  if constexpr (match::push_back<std::decay_t<decltype(a)>>) {
    assert((int)a.size() >= 1);

    a[i] = std::move(a.back());
    a.pop_back();

    assert((int)a.size() >= 0);
  }
  // else...
};

static auto remove = [](auto& data, auto&& h) {
  using item_t = typename std::decay_t<decltype(h)>::type;
  static_assert(!match::not_removable<item_t>,
                "storage::remove: this item lives as long as the store (removable = false)");

  if constexpr (requires { data.store()._registry; }) {
    // registry world: walk the columns listed for item_t at registration,
    // dynamic properties included (their entries follow the swap and pop)
    auto h1 = h;  // clang 19
    if constexpr (has_del(h1)) {
      h1.__del__();
    }
    data.store()._registry.template remove_item<item_t>(h.index().value());
  } else {
    using info_t = get_info_t<decltype(data)>;
    using all_time_levels_t = decltype(all_properties_as<property::time_levels>(data));

    using indice = typename info_t::template env<item_t>::indice;

    using all_attrs_t = decltype(mp::tuple_unique(
        concat(attributes(item_t{}), attached_storages(h.item_type(), data))));

    // the dynamic properties have no static slot: their sparse sets live in
    // the _dynamic_properties map
    auto attrs = mp::filter(all_attrs_t{}, mp::is_a_model<[]<typename T>() constexpr {
                              return !requires { typename T::dynamic_storage_t; };
                            }>);
    using dynamic_attrs_t =
        decltype(mp::filter(all_attrs_t{}, mp::is_a_model<[]<typename T>() constexpr {
                              return requires { typename T::dynamic_storage_t; };
                            }>));

    // the number of items before the removal, for the refs to them
    std::size_t count = ref_remap::invalid;
    if constexpr (mp::size(attrs) > mp::size_c<0>) {
      using first_t = std::decay_t<decltype(attrs[0_c])>;
      const auto& first_storage = memory(0, mp::get<first_t>(data.store()));
      if constexpr (match::push_back<std::decay_t<decltype(first_storage)>>) {
        count = std::size(first_storage);
      }
    }

    if (count != ref_remap::invalid && h.index().value() >= count) [[unlikely]] {
      throw std::out_of_range("storage::remove: index " + std::to_string(h.index().value()) +
                              " past the end (" + std::to_string(count) + " items)");
    }

    // call handle delete function if present
    auto h1 = h;  // clang 19
    if constexpr (has_del(h1)) {
      h1.__del__();
    }

    // the removed item's entries go, the last item's entries are re-keyed to
    // its index, as the storages below are moved back (see registry column)
    if constexpr (mp::size(dynamic_attrs_t{}) > mp::size_c<0>) {
      static_assert(mp::size(attrs) > mp::size_c<0>,
                    "dynamic properties of an item without static storage");
      using first_t = std::decay_t<decltype(attrs[0_c])>;
      const auto& first_storage = memory(0, mp::get<first_t>(data.store()));
      if constexpr (match::push_back<std::decay_t<decltype(first_storage)>>) {
        using env_t = typename info_t::template env<item_t>;
        const indice i = h.index().value();
        const indice last = static_cast<indice>(std::size(first_storage) - 1);
        mp::for_each(dynamic_attrs_t{}, [&data, i, last]<typename P>(P) {
          auto& set = prop_dynamic_map<P, env_t, indice>(data);
          set.erase(i);
          if (last != i && set.contains(last)) {
            auto moved = std::move(set.at(last));
            set.erase(last);
            set.emplace(i, std::move(moved));
          }
        });
      }
    }

    if constexpr (mp::size(attrs) > mp::size_c<0>) {
      mp::for_each(attrs, [&data, &h]<match::attribute A>(A) {
        return mp::for_each(
            mp::range<memory_size<A, all_time_levels_t>()>, [&data, &h](indice step) {
              move_back(h.index().value(), memory(step, mp::get<A>(data.store())));
            });
      });
    }

    // the refs to the removed item become invalid, those to the last one
    // follow it (the registry does it in remove_item)
    if (count != ref_remap::invalid && count > 0) {
      remap_refs<item_t>(data, ref_remap::swap_and_pop(h.index().value(), count - 1));
    }
  }
};

/**
 * Remove the items of Item of the given indices: as many swaps and pops,
 * with the refs to them rewritten once (registry world) instead of once per
 * removal. The indices are those before the call; duplicates are ignored.
 */
template <typename Item>
void remove_many(auto& data, std::vector<std::size_t> indices) {
  static_assert(
      !match::not_removable<Item>,
      "storage::remove_many: this item lives as long as the store (removable = false)");
  std::sort(indices.begin(), indices.end(), std::greater<>{});
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  using indice = typename get_info_t<std::decay_t<decltype(data)>>::template env<Item>::indice;
  if constexpr (requires { data.store()._registry; }) {
    using handle_t = decltype(make_handle(data, index<Item, indice>{}));
    if constexpr (decltype(has_del(std::declval<handle_t&>()))::value) {
      for (auto i : indices) {
        auto h = make_handle(data, index<Item, indice>(static_cast<indice>(i)));
        h.__del__();
      }
    }
    data.store()._registry.template remove_items<Item>(std::move(indices));
  } else {
    // from the largest index down: the last item is never one to remove
    for (auto i : indices) {
      remove(data, make_handle(data, index<Item, indice>(static_cast<indice>(i))));
    }
  }
}

}  // namespace siconos::storage
