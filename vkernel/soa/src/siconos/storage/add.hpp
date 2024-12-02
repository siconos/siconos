#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/memory.hpp"

#include <assert.h>

namespace siconos::storage {

using namespace siconos::storage::pattern;

template <match::item Item>
static auto add = [](auto&& data) constexpr -> decltype(auto) {
  using data_t = std::decay_t<decltype(data)>;
  using info_t = get_info_t<data_t>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;
  constexpr auto attached_storage = ground::filter(
      typename info_t::all_properties_t{}, is_attached_storage<Item>);

  constexpr auto attrs =
      ground::tuple_unique(concat(attributes(Item{}), attached_storage));

  using attrs_t = std::decay_t<decltype(attrs)>;

  if constexpr (ground::size(attrs_t{}) > ground::size_c<0>) {
    indice index = 0;
    ground::for_each(attrs, [&data, &index]<match::attribute A>(A) {
      ground::for_each(
          ground::range<memory_size<A, all_keeps_t>()>,
          [&data, &index](auto step) {
            auto&& storage = memory(
                step, ground::get<A>(static_cast<data_t&&>(data).store()));

            using storage_t = std::decay_t<decltype(storage)>;

            if constexpr (match::push_back<storage_t>) {
              storage.push_back(typename storage_t::value_type{});
              assert(index > 0 ? index == std::size(storage) - 1
                               : index == 0);
              index = std::size(storage) - 1;
            }
            else {
              storage[0] = typename storage_t::value_type{};
              index = 0;
            }
          });
    });
    return make_full_handle<Item>(data, index);
  }
  else {
    return make_full_handle<Item>(data, indice{0});
  }
};

}  // namespace siconos::storage
