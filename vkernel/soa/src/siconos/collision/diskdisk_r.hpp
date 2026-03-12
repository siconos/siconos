#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

struct diskdisk_r : item, model::relation2, model::any_lagrangian_relation {
  using dof = some::indice_parameter<"dof">;

  struct attributes {}; // empty

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    template <match::handle<model::lagrangian_ds> DS1,
              match::handle<model::lagrangian_ds> DS2>
    decltype(auto) compute_h(DS1& ds1, DS2& ds2)
    {
      auto& q1 = storage::attr<"q">(ds1);
      auto& q2 = storage::attr<"q">(ds2);

      auto& r1 = storage::make_handle(self()->data(), storage::prop<"shape">(ds1))
                     .radius();
      auto& r2 = storage::make_handle(self()->data(), storage::prop<"shape">(ds2))
                     .radius();

      auto dx = q2[0] - q1[0];
      auto dy = q2[1] - q1[1];

      return algebra::hypot(dx, dy) - (r1 + r2);
    }

    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::lagrangian_ds> DS2, typename M>
    void compute_jachq(S step, DS1& ds1, DS2& ds2, M& h_matrix1, M& h_matrix2)
    {
      auto& data = self()->data();

      auto& q1 = storage::attr<"q">(ds1);
      auto& q2 = storage::attr<"q">(ds2);

      auto& r1 = storage::make_handle(data, storage::prop<"shape">(ds1)).radius();
      auto& r2 = storage::make_handle(data, storage::prop<"shape">(ds2)).radius();

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
      auto& g2 = h_matrix2;

      g1(0, 0) = -dxsd;
      g1(1, 0) = dysd;
      g1(0, 1) = -dysd;
      g1(1, 1) = -dxsd;
      g1(0, 2) = 0.;
      g1(1, 2) = -r1;

      g2(0, 0) = dxsd;
      g2(1, 0) = -dysd;
      g2(0, 1) = dysd;
      g2(1, 1) = dxsd;
      g2(0, 2) = 0.;
      g2(1, 2) = -r2;
    }
  };

  //             method<"compute_jachq", &interface<H>::compute_jachq>>;
  // cf
  // https://stackoverflow.com/questions/47906426/type-of-member-functions-arguments
};
}  // namespace siconos::collision
