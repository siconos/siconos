#pragma once

#include "siconos/storage/add.hpp"
#include "siconos/storage/data_holder.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/info.hpp"
#include "siconos/storage/make.hpp"
#include "siconos/storage/memory.hpp"
#include "siconos/storage/permute.hpp"
#include "siconos/storage/properties.hpp"
#include "siconos/storage/registry.hpp"
#include "siconos/storage/remove.hpp"
#include "siconos/storage/sparse_set.hpp"

namespace siconos::storage {

// template <match::item T>
// static constexpr void for_each_attribute(T)
// {
//   return mp::compose(mp::for_each, attributes)(T{});
// };

using pattern::attr_t;
using pattern::wrap;

}  // namespace siconos::storage
