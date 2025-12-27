#pragma once

#include <concepts>
#include <type_traits>  // dependent_false

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/chained_segment.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/collision/translated.hpp"
#include "siconos/model/lagrangian_ds.hpp"

namespace siconos::collision {

// A point linked to an item (a dynamical system or a shape),
// by default is not defined.
template <match::item Item>
struct point : item {
  static_assert(always_false<Item>, "point is not defined for this type");
};

// Direct association to a lagrangian_ds: just a point associated at
// the coordinates of the system.
template <match::item Item>
  requires std::derived_from<Item, model::lagrangian_ds>
struct point<Item> : item {
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
    decltype(auto) item()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"item">(*self()));
    };

    void update(auto step)
    {
      // one body / one point
      storage::attr<"coord">(*self())[0] =
          storage::attr<"q">(item(), step)(0);
      storage::attr<"coord">(*self())[1] =
          storage::attr<"q">(item(), step)(1);
      storage::attr<"coord">(*self())[2] = 0.; /* 2D */
    }
  };
};

// Association to simple shapes.
template <match::item Item>
  requires(std::derived_from<Item, translated<shape::disk>> ||
           std::derived_from<Item, shape::segment>)
struct point<Item> {
  using items = gather<Item>;
  using item_t = Item;

  using attributes = gather<
      // 3D coordinates, 2D => last value = 0.
      attribute<"flag", some::boolean>,
      attribute<"coord", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"item", some::item_ref<item_t>>,  // reverse link
      attribute<"point_index", some::indice>>;    // index of the point on the
                                                  // disk

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"item">(*self()));
    };
    decltype(auto) point_index()
    {
      return storage::attr<"point_index">(*self());
    };

    void update(auto step)
    {
      // associated to static shapes for the moment
      // coord() = item().point_coord(point_index());
    }
  };
};

// association to compound shapes.
template <match::item Item>
  requires std::derived_from<Item, shape::chained_segment>
struct point<Item> {
  using items = gather<Item>;
  using item_t = Item;

  using attributes = gather<
      // 3D coordinates, 2D => last value = 0.
      attribute<"flag", some::boolean>,
      attribute<"coord", some::vector<some::scalar, some::indice_value<3>>>,
      attribute<"item", some::item_ref<item_t>>,  // reverse link
      attribute<"point_index", some::indice>,
      attribute<"seg_index", some::indice>>;  // index of segment in the mesh

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"item">(*self()));
    };
    decltype(auto) point_index()
    {
      return storage::attr<"point_index">(*self());
    };
    decltype(auto) seg_index()
    {
      return storage::attr<"seg_index">(*self());
    };

    void update(auto step)
    {
      coord() = item().point_coord(point_index(), seg_index());
    }
  };
};

}  // namespace siconos::collision
