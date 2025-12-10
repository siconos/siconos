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

    // access to the translated shape
    decltype(auto) translated()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"translated">(*self()));
    };

    // translated points coordinates
    decltype(auto) point_coord(auto point_index)
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using vector_t = typename env_t::template vector<scalar, 3>;

      vector_t coord = translation() + translated().point_coord(point_index);

      return coord;
    }

    decltype(auto) points_coords()
    {
      using env_t = decltype(self()->env());
      using indice_t = typename env_t::indice;

      return view::iota((indice_t)0, translated().maxpoints()) |
             view::transform([this](auto i) { return point_coord(i); });
    }

    decltype(auto) translation()
    {
      return storage::attr<"translation">(*self());
    }
  };
};
}  // namespace siconos::collision
