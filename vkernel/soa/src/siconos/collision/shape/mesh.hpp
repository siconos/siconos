#pragma once

#include "siconos/collision/collision.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/chained_segment.hpp"

namespace siconos::collision::shape {

struct mesh : item {
  using items = gather<chained_segment>;

  struct attributes {
    some::item_ref<chained_segment> segments;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) segments()
    {
      return storage::handle(self()->data(),
                             storage::attr<"segments">(*self()));
    };

    decltype(auto) segment(auto& indice)
    {
      return storage::handle(self()->data(), segments()[indice]);
    };
  };
};

}  // namespace siconos::collision::shape
