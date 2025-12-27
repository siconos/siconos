#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/disk.hpp"
#include "siconos/collision/translated.hpp"
// #include "siconos/model/model.hpp"

namespace siconos::collision {

// disk fixed disk

struct diskfdisk_r : item,
                     model::relation1,
                     model::any_lagrangian_relation {
  using attributes =
      gather<attribute<"translated_disk_shape",
                       some::item_ref<collision::translated<shape::disk>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) translated_disk_shape()
    {
      return make_ref_handle(self()->data(),
                             attr<"translated_disk_shape">(*self()));
    };

    decltype(auto) shape() { return self()->translated_disk_shape(); }

    decltype(auto) compute_h(auto& ds)
    {
      auto& data = self()->data();
      auto& q = storage::attr<"q">(ds);
      return collision::distance(q, translated_disk_shape().translation()) -
             make_handle(data, prop<"shape">(ds)).radius() -
             translated_disk_shape().translated().radius();
    }

    void compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      auto& data = self()->data();
      using scalar = typename decltype(self()->env())::scalar;

      const auto& q1 = storage::attr<"q">(ds);
      const auto& q2 = translated_disk_shape().translation();
      const scalar& r1 =
          storage::make_handle(data, storage::prop<"shape">(ds)).radius();

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
    }
  };
};
}  // namespace siconos::collision
