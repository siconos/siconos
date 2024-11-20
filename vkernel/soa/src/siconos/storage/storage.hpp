#pragma once

#include "siconos/storage/add.hpp"
#include "siconos/storage/data_holder.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/remove.hpp"

namespace siconos::storage {

static auto apply_fun = []<typename Item, typename SomeFun>(
                            auto& data, Item, SomeFun&& some_fun) {
  using item_t = Item;
  using info_t = get_info_t<decltype(data)>;
  using all_keeps_t = decltype(all_properties_as<property::keep>(data));

  using indice = typename info_t::env::indice;

  auto attrs = ground::tuple_unique(
      concat(attributes(item_t{}), attached_storages(item_t{}, data)));

  if constexpr (ground::size(attrs) > ground::size_c<0>) {
    ground::for_each(attrs, [&data, &some_fun]<match::attribute A>(A) {
      return ground::for_each(ground::range<memory_size<A, all_keeps_t>()>,
                              [&data, &some_fun](indice step) {
                                static_cast<SomeFun&&>(some_fun)(memory(
                                    step, ground::get<A>(data.store())));
                              });
    });
  }
};

template <match::item T>
static constexpr void for_each_attribute(T)
{
  return ground::compose(ground::for_each, attributes)(T{});
};

using pattern::attr_t;
using pattern::wrap;

}  // namespace siconos::storage
