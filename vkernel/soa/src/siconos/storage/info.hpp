#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

using namespace pattern;
// the storage info key
struct info {};

template <typename D>
auto get_info(D&& data)
{
  if constexpr (match::store<std::decay_t<D>>) {
    // database case
    return mp::get<info>(data.store());
  }
  else {
    // pre-map case
    return mp::second(mp::at(data, mp::size_c<0>));
  }
}

template <typename D>
using get_info_t = decltype(get_info(std::decay_t<D>{}));

}  // namespace siconos::storage
