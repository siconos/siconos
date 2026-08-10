#pragma once

#include <concepts>
#include <type_traits>  // dependent_false

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/shape/mesh.hpp"
#include "siconos/collision/shape/segment.hpp"
#include "siconos/collision/translated.hpp"
#include "siconos/model/lagrangian_ds.hpp"

namespace siconos::collision {

struct empty_shape : empty_item {
  using empty_shape_t = void;
};

// A point linked to an item (a dynamical system or a shape),
// by default is not defined.
template <match::item Item, match::item Shape>
struct point : item {
  static_assert(always_false<Item, Shape>,
                "point is not defined for this association");
};

// Direct association to a lagrangian_ds: just a point associated at
// the coordinates of the system.
template <match::item Item>
  requires(requires { typename Item::lagrangian_dynamical_system_t; })
struct point<Item, empty_shape> : item {
  using items = gather<Item>;
  using item_t = Item;

  struct attributes {
    some::boolean flag;
    some::vector<float, some::indice_value<3>> coord;
    some::item_ref<item_t> item;
  };

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
      coord()[0] = item().q(step)(0);
      coord()[1] = item().q(step)(1);
      coord()[2] = 0; /* 2D */
    }
  };
};

// Association to static shapes.
template <match::item Shape>
  requires(std::derived_from<Shape, translated<shape::disk>> ||
           std::derived_from<Shape, shape::segment>)
struct point<empty_item, Shape> {
  using items = gather<Shape>;
  using item_t = Shape;

  struct attributes {
    some::boolean flag;
    some::vector<float, some::indice_value<3>> coord;
    some::item_ref<item_t> item;
    some::indice point_index;
  };

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

// association to mesh.
template <match::item Item, match::item Shape>
  requires(std::derived_from<Shape, shape::mesh> &&
           (std::derived_from<Item, model::lagrangian_ds> ||
            std::derived_from<Item, model::elastic_lagrangian_ds>))
struct point<Item, Shape> {
  using item_t = Item;
  using shape_t = Shape;
  using items = gather<item_t, shape_t>;

  struct attributes {
    some::boolean flag;
    some::vector<float, some::indice_value<3>> coord;
    some::item_ref<item_t> item;
    some::item_ref<shape_t> shape;
    some::indice point_index;
    some::indice seg_index;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) item()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"item">(*self()));
    }
    decltype(auto) shape()
    {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"shape">(*self()));
    }

    decltype(auto) point_index()
    {
      return storage::attr<"point_index">(*self());
    }
    decltype(auto) seg_index() { return storage::attr<"seg_index">(*self()); }

    void update(auto step)
    {
      // coord() = algebra::cast(
      //     mp::type_c<float>,
      //     storage::make_handle(self()->data(),
      //                          storage::prop<"shape">(self()->item()))
      //         .segments()
      //         .point_coord(point_index(), seg_index()));
    }
  };
};

}  // namespace siconos::collision
