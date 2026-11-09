#pragma once

#include <cmath>
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
  static_assert(always_false<Item, Shape>, "point is not defined for this association");
};

// Direct association to a lagrangian_ds: a point at the coordinates of
// the system, shifted by an offset. A small body has one point with a
// null offset; a body too big for the neighborhood radius gets several
// points on its boundary (see space_filter::make_points). The offset of a
// disk is in the global frame (a circle is invariant by rotation), the
// offset of a polytope in the body frame (body_frame): it turns with it.
template <match::item Item>
  requires(requires { typename Item::lagrangian_dynamical_system_t; })
struct point<Item, empty_shape> : item {
  using items = gather<Item>;
  using item_t = Item;
  // moving lagrangian body: its centre is re-sampled every step
  static constexpr bool is_dynamic_v = true;

  struct attributes {
    some::boolean flag;
    some::vector<float, some::indice_value<3>> coord;
    some::vector<float, some::indice_value<2>> offset;
    some::boolean body_frame;
    // position of the body, and its bounding radius around it: disk radius
    // or polytope circumradius, infinity for a body of unknown shape
    some::vector<float, some::indice_value<2>> center;
    some::scalar bound;
    some::item_ref<item_t> item;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) flags() { return storage::attr<"flag">(*self()); };
    decltype(auto) coord() { return storage::attr<"coord">(*self()); };
    decltype(auto) offset() { return storage::attr<"offset">(*self()); };
    decltype(auto) body_frame() { return storage::attr<"body_frame">(*self()); };
    decltype(auto) center() { return storage::attr<"center">(*self()); };
    decltype(auto) bound() { return storage::attr<"bound">(*self()); };
    decltype(auto) item() {
      return storage::make_ref_handle(self()->data(), storage::attr<"item">(*self()));
    };

    void update(auto step) {
      const auto& q = item().q(step);
      center()[0] = q(0);
      center()[1] = q(1);
      if (body_frame()) {
        const auto c = std::cos(q(2));
        const auto s = std::sin(q(2));
        coord()[0] = q(0) + c * offset()[0] - s * offset()[1];
        coord()[1] = q(1) + s * offset()[0] + c * offset()[1];
      } else {
        coord()[0] = q(0) + offset()[0];
        coord()[1] = q(1) + offset()[1];
      }
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
  // fixed shape: probes never move, so CompactNSearch skips position updates
  static constexpr bool is_dynamic_v = false;

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
    decltype(auto) item() {
      return storage::make_ref_handle(self()->data(), storage::attr<"item">(*self()));
    };
    decltype(auto) point_index() { return storage::attr<"point_index">(*self()); };

    void update(auto step) {
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
  // mesh of a moving (FEM) body: re-sampled every step
  static constexpr bool is_dynamic_v = true;

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
    decltype(auto) item() {
      return storage::make_ref_handle(self()->data(), storage::attr<"item">(*self()));
    }
    decltype(auto) shape() {
      return storage::make_ref_handle(self()->data(), storage::attr<"shape">(*self()));
    }

    decltype(auto) point_index() { return storage::attr<"point_index">(*self()); }
    decltype(auto) seg_index() { return storage::attr<"seg_index">(*self()); }

    void update(auto step) {
      // coord() = algebra::cast(
      //     mp::type_c<float>,
      //     storage::make_handle(self()->data(),
      //                          storage::prop<"shape",
      //                          storage::config_declared>(self()->item()))
      //         .segments()
      //         .point_coord(point_index(), seg_index()));
    }
  };
};

}  // namespace siconos::collision
