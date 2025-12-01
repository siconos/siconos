#pragma once
// translation only (for symmetric items)

#include "collision.hpp"
#include "collision_head.hpp"

namespace siconos::collision {

template <match::item Item>
struct translated : item<> {
  using item_t = Item;

  using attributes =
      gather<attribute<"translated", some::item_ref<item_t>>,
             attribute<"translation",
                       some::vector<some::scalar, some::indice_value<3>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) translated()
    {
      return storage::handle(self()->data(),
                             storage::attr<"translated">(*self()));
    };
    decltype(auto) translation()
    {
      return storage::attr<"translation">(*self());
    }
  };
};
}  // namespace siconos::collision
