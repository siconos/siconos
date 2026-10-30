#pragma once

#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"

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
  }
};

}  // namespace siconos::storage
