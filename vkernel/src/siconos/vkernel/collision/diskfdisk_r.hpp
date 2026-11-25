// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/algebra/math.hpp"
#include "siconos/vkernel/collision/collision_head.hpp"
#include "siconos/vkernel/collision/shape/disk.hpp"
#include "siconos/vkernel/collision/shape_of.hpp"
#include "siconos/vkernel/collision/translated.hpp"
// #include "siconos/vkernel/model/model.hpp"

namespace siconos::vkernel::collision {

// disk fixed disk

struct diskfdisk_r : item, model::relation1, model::any_lagrangian_relation {
  struct attributes {
    some::item_ref<collision::translated<shape::disk>> translated_disk_shape;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) translated_disk_shape() {
      return make_ref_handle(self()->data(), attr<"translated_disk_shape">(*self()));
    };

    decltype(auto) shape() { return self()->translated_disk_shape(); }

    [[gnu::always_inline]] decltype(auto) compute_h(auto step, auto& ds) {
      auto& q = storage::attr<"q">(ds, step);
      return collision::distance(q, translated_disk_shape().translation()) -
             collision::disk_radius(ds) - translated_disk_shape().translated().radius();
    }

    [[gnu::always_inline]] void compute_jachq(auto step, auto& ds, auto& h_matrix1) {
      using scalar = typename decltype(self()->env())::scalar;

      const auto& q1 = storage::attr<"q">(ds, step);
      const auto& q2 = translated_disk_shape().translation();
      const scalar r1 = collision::disk_radius(ds);

      //      const scalar& r2 =
      //        translated_disk_shape().item().radius();

      auto x1 = q1(0);
      auto y1 = q1(1);
      auto x2 = q2(0);
      auto y2 = q2(1);

      auto dx = x2 - x1;
      auto dy = y2 - y1;

      auto d = algebra::hypot(dx, dy);

      auto dxsd = dx / d;
      auto dysd = dy / d;

      auto& g1 = h_matrix1;

      g1(0, 0) = -dxsd;
      g1(1, 0) = dysd;
      g1(0, 1) = -dysd;
      g1(1, 1) = -dxsd;
      g1(0, 2) = 0.;
      g1(1, 2) = -r1;

      // the rolling row of a law with rolling resistance: omega1, the
      // fixed disk does not turn
      if constexpr (std::decay_t<decltype(g1)>::RowsAtCompileTime > 2) {
        g1(2, 2) = 1.;
      }
    }
  };
};
}  // namespace siconos::vkernel::collision
