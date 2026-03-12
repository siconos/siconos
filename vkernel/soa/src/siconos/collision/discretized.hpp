#pragma once

#include <numbers>

#include "collision_head.hpp"
#include "collision.hpp"

namespace siconos::collision {

struct discretization : item {
  using attributes = gather<attribute<"maxpoints", some::scalar>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) maxpoints() { return attr<"maxpoints">(*self()); }
  };
};

template <match::item Item>
struct discretized : composite_item<Item, discretization> {};
}  // namespace siconos::collision
