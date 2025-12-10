#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/base_concepts.hpp"
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
  using info_t = get_info_t<decltype(data)>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs = ground::tuple_unique(
      concat(attributes(item_t{}), attached_storages(h.item_type(), data)));

  // call handle delete function if present
  auto h1 = h;  // clang 19
  if constexpr (has_del(h1)) {
    h1.__del__();
  }

  if constexpr (ground::size(attrs) > ground::size_c<0>) {
    ground::for_each(attrs, [&data, &h]<match::attribute A>(A) {
      return ground::for_each(
          ground::range<memory_size<A, all_keeps_t>()>,
          [&data, &h](indice step) {
            move_back(h.index().value(),
                      memory(step, ground::get<A>(data.store())));
          });
    });
  }
};

}  // namespace siconos::storage
