#pragma once

#include <concepts>

#include "siconos/collision/collision_head.hpp"
#include "siconos/model/lagrangian_ds.hpp"

#include "siconos/collision/shape/chained_segment.hpp"
namespace siconos::collision {

// a point linked to an item (a dynamical system or a shape)
template <match::item Item>
struct point : item<> {
  using items = gather<Item>;

  using item_t = Item;
  using attributes = gather<
      // 3D coordinates, 2D => last value = 0.
      attribute<"flag", some::boolean>,
      attribute<"coord", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"item", some::item_ref<item_t>>>;  // reverse link

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item() { return storage::attr<"item">(*self()); };

    void update(auto step)
    {
      auto& data = self()->data();

      if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
        // one body / one point
        auto hbody = storage::handle(data, item());
        storage::attr<"coord">(*self())[0] = storage::attr<"q">(hbody, step)(0);
        storage::attr<"coord">(*self())[1] = storage::attr<"q">(hbody, step)(1);
        storage::attr<"coord">(*self())[2] = 0.; /* 2D */
      }
    }
  };
};

template <>
struct point<collision::shape::chained_segment> {
  using items = gather<collision::shape::chained_segment>;

  using item_t = collision::shape::chained_segment;
  using attributes = gather<
      // 3D coordinates, 2D => last value = 0.
      attribute<"flag", some::boolean>,
      attribute<"coord", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"item", some::item_ref<item_t>>,  // reverse link
      attribute<"index", some::indice>>;  // index of segment in the mesh

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item() { return storage::attr<"item">(*self()); };
    decltype(auto) index() { return storage::attr<"index">(*self()); };

    void update(auto step)
    {
      auto& data = self()->data();

      auto index = self()->index();
      if constexpr (std::derived_from<item_t, model::lagrangian_ds>) {
        // one body / one point
        auto hbody = storage::handle(data, item());
        storage::attr<"coord">(*self())[0] = storage::attr<"q">(hbody, step)(index);
        storage::attr<"coord">(*self())[1] =
          storage::attr<"q">(hbody, step)(index + 1);
        storage::attr<"coord">(*self())[2] = 0.; /* 2D */
      }
    }
  };
};

}  // namespace siconos::collision
