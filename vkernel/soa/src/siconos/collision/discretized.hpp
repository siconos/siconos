#pragma once

#include <numbers>

#include "collision.hpp"
#include "collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"

namespace siconos::collision {

template <match::item Item>
struct discretized : Item {
  using items = gather<Item>;
  using item_t = Item;

  using attributes = decltype(storage::pattern::concat(
      typename item_t::attributes{},
      gather<attribute<"maxpoints", some::scalar>>{}));

  template <typename Handle>
  struct interface : item_t::template interface<Handle> {
    using item_t::template interface<Handle>::self;

    struct addons {
      decltype(auto) maxpoints() { return attr<"maxpoints">(*self()); }

      decltype(auto) point_coord(auto point_index)
      {
        return self()->point_coord(maxpoints(), point_index);
      }
    };
    addons _addons;
  };
};
}  // namespace siconos::collision
