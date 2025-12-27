#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/segment.hpp"
// #include "siconos/model/model.hpp"

namespace siconos::collision {

// disk fixed segment

struct diskfsegment_r : item,
                        model::relation1,
                        model::any_lagrangian_relation {
  using attributes =
      gather<attribute<"segment", some::item_ref<shape::segment>>>;

  using properties =
      gather<storage::attached<shape::segment, symbol<"relation">,
                               some::item_ref<diskfsegment_r>>>;
  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) segment()
    {
      return make_ref_handle(self()->data(), attr<"segment">(*self()));
    }

    decltype(auto) shape() { return self()->segment(); }

    decltype(auto) compute_h(auto& ds)
    {
      auto& q = storage::attr<"q">(ds);
      return segment().distance(q) -
             make_handle(self()->data(), prop<"shape">(ds)).radius();
    }

    void compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      auto& data = self()->data();
      using scalar = typename decltype(self()->env())::scalar;

      const auto& q = storage::attr<"q">(ds);
      const scalar& r =
          storage::make_handle(data, storage::prop<"shape">(ds)).radius();

      const scalar& x = q(0);
      const scalar& y = q(1);

      const scalar& x1 = segment().x1();
      const scalar& x2 = segment().x2();
      const scalar& y1 = segment().y1();
      const scalar& y2 = segment().y2();

      const scalar tmp0 = x1 - x2;
      const scalar tmp1 = y1 - y2;
      const scalar tmp2 = 1.0 / (pow(tmp0, 2) + pow(tmp1, 2));
      const scalar tmp3 = x - x1;
      const scalar tmp4 = y - y1;
      const scalar tmp5 = tmp2 * (-tmp0 * tmp3 - tmp1 * tmp4);
      const scalar tmp6 = fmin(1, tmp5);
      const scalar tmp7 = fmax(0, tmp6);
      const scalar tmp8 = tmp0 * tmp7 + tmp3;
      const scalar tmp9 = tmp1 * tmp7 + tmp4;
      const scalar tmp10 = pow(pow(tmp8, 2) + pow(tmp9, 2), -1.0 / 2.0);
      const scalar tmp11 =
          ((tmp5 > 1) ? (0) : ((tmp5 == 1) ? (1.0 / 2.0) : (1)));
      const scalar tmp12 =
          ((tmp6 < 0) ? (0) : ((tmp6 == 0) ? (1.0 / 2.0) : (1)));
      const scalar tmp13 = -tmp0 * tmp11 * tmp12 * tmp2;
      const scalar tmp14 = -tmp1 * tmp11 * tmp12 * tmp2;

      h_matrix1(0, 0) = tmp10 * (tmp1 * tmp13 * tmp9 +
                                 (1.0 / 2.0) * tmp8 * (2 * tmp0 * tmp13 + 2));
      h_matrix1(1, 0) =
          -tmp10 *
          (tmp0 * tmp14 * tmp8 + (1.0 / 2.0) * tmp9 * (2 * tmp1 * tmp14 + 2));
      h_matrix1(0, 1) = -h_matrix1(1, 0);
      h_matrix1(1, 1) = h_matrix1(0, 0);
      h_matrix1(0, 2) = 0;
      h_matrix1(1, 2) = -r;
    }
  };
};
}  // namespace siconos::collision
