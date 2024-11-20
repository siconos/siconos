#pragma once

#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

namespace match = pattern::match;
// the storage info key
struct info {};

template <typename D>
auto get_info(D&& data)
{
  if constexpr (match::store<std::decay_t<D>>) {
    // database case
    return ground::get<info>(data.store());
  }
  else {
    // pre-map case
    return ground::second(ground::at(data, ground::size_c<0>));
  }
}

template <typename D>
using get_info_t = decltype(get_info(std::decay_t<D>{}));

}  // namespace siconos::storage
