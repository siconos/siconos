#pragma once

#include <assert.h>

#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/properties.hpp"

namespace siconos::storage {

using namespace siconos::storage::pattern;

template <match::item Item>
decltype(auto) add(auto&& data) {
  // the columns of Item were listed at registration, so they are walked at
  // run time, with no compile-time scan of the world
  using indice = typename get_info_t<std::decay_t<decltype(data)>>::template env<Item>::indice;
  const auto index = static_cast<indice>(data.store()._registry.template push_item<Item>());
  auto h = make_full_handle<Item>(data, index);
  // handle initialization is not mandatory
  if constexpr (has_init(h)) {
    h.__init__();
  }
  return h;
}

// Multiple add at once.
template <match::item Item>
decltype(auto) add(auto&& data, std::size_t count) {
  return view::iota(std::size_t{0}, count) |
         view::transform([&data](auto&&) { return add<Item>(data); });
}

}  // namespace siconos::storage
