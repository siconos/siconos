#pragma once

#include <array>

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/pattern/pattern.hpp"

namespace siconos::storage {

using namespace pattern;

template <typename T, std::size_t N>
using memory_t = std::array<T, N>;

static constexpr auto memory = []<typename T>(
                                   typename std::decay_t<T>::size_type step,
                                   T&& mem) constexpr -> decltype(auto) {
  return static_cast<T&&>(mem)[step % std::size(static_cast<T&&>(mem))];
};

template <match::attribute Attr, typename keeps_t>
static constexpr std::size_t memory_size()
{
  auto tpl = mp::filter(keeps_t{},
                            mp::is_a_model<[]<typename T>() consteval {
                              return std::is_same_v<Attr, typename T::type>;
                            }>);
  if constexpr (mp::size(tpl) > mp::size_c<0>) {
    return tpl[0_c].size;
  }
  else {
    // memory size not specified
    return 1;
  }
};

}  // namespace siconos::storage
