#pragma once

#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {
using namespace storage;
using namespace storage::pattern;

template <typename...>
constexpr bool always_false = false;

}  // namespace siconos::collision
