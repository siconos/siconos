#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/chained_segment.hpp"

namespace siconos::collision {

struct diskmesh_r : item, model::relation2, model::any_lagrangian_relation {
  using attributes =
      gather<attribute<"mesh", some::item_ref<shape::chained_segment>>,
             attribute<"contact_index", some::indice>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mesh()
    {
      return storage::make_ref_handle(self()->data(), attr<"mesh">(*self()));
    }
    decltype(auto) shape() { return self()->mesh(); }
    decltype(auto) contact_index()
    {
      return storage::attr<"contact_index">(*self());
    }
    decltype(auto) compute_h(auto& ds1, auto& ds2)
    {
      auto& q1 = storage::attr<"q">(ds1); /* disk */

      // auto& mesh = storage::prop<"shape">(ds2);
      // assert ds2 mesh == mesh()

      auto& cindex = self()->contact_index();

      return mesh().distance(q1, cindex) -
             make_handle(self()->data(), prop<"shape">(ds1)).radius();
    }

    template <typename I, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::rt_lagrangian_ds> DS2, typename M1,
              typename M2>
    void compute_jachq(I step, DS1& ds1, DS2& ds2, M1& h_matrix1,
                       M2& h_matrix2)
    {
      auto& data = self()->data();
      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = storage::attr<"q">(ds1);
      auto& q2 = storage::attr<"q">(ds2);

      auto& r =
          storage::make_handle(data, storage::prop<"shape">(ds1)).radius();

      /* disk coordinates */
      const scalar& x = q1(0);
      const scalar& y = q1(1);

      /* segment end points */
      const scalar& x1 = q2(self()->contact_index());
      const scalar& y1 = q2(self()->contact_index() + 1);

      const scalar& x2 = q2(self()->contact_index() + 2);
      const scalar& y2 = q2(self()->contact_index() + 3);

      // auto& g1 = h_matrix1;
      // auto& g2 = h_matrix2;

      double tmp0 = x1 - x2;
      double tmp1 = y1 - y2;
      double tmp2 = tmp0 * tmp0 + tmp1 * tmp1;
      double tmp3 = 1.0 / tmp2;
      double tmp4 = -tmp0 * tmp3;
      double tmp5 = x - x1;
      double tmp6 = y - y1;
      double tmp7 = -tmp0 * tmp5 - tmp1 * tmp6;
      double tmp8 = tmp3 * tmp7;
      double tmp9 = fmin(1, tmp8);
      double tmp10 = fmax(0, tmp9);
      double tmp11 = tmp1 * tmp10 + tmp6;
      double tmp12 = 2 *
                     ((tmp9 < 0) ? (0) : ((tmp9 == 0) ? (1.0 / 2.0) : (1))) *
                     ((tmp8 > 1) ? (0) : ((tmp8 == 1) ? (1.0 / 2.0) : (1)));
      double tmp13 = tmp1 * tmp12;
      double tmp14 = tmp11 * tmp13;
      double tmp15 = tmp0 * tmp10 + tmp5;
      double tmp16 = tmp0 * tmp12;
      double tmp17 =
          (1.0 / 2.0) * pow(tmp11 * tmp11 + tmp15 * tmp15, -1.0 / 2.0);
      double tmp18 = -tmp1 * tmp3;
      double tmp19 = tmp15 * tmp16;
      double tmp20 = 2 * x1;
      double tmp21 = tmp20 - 2 * x2;
      double tmp22 = tmp7 * 1.0 / (tmp2 * tmp2);
      double tmp23 = -tmp21 * tmp22 + tmp3 * (tmp20 - x - x2);
      double tmp24 = 2 * tmp10;
      double tmp25 = tmp24 - 2;
      double tmp26 = 2 * y1;
      double tmp27 = tmp26 - 2 * y2;
      double tmp28 = -tmp22 * tmp27 + tmp3 * (tmp26 - y - y2);
      double tmp29 = tmp21 * tmp22 + tmp3 * tmp5;
      double tmp30 = -tmp24;
      double tmp31 = tmp22 * tmp27 + tmp3 * tmp6;
      h_matrix1(0, 0) = -tmp17 * (tmp14 * tmp4 + tmp15 * (tmp16 * tmp4 + 2));
      h_matrix1(0, 1) =
          -tmp17 * (tmp11 * (tmp13 * tmp18 + 2) + tmp18 * tmp19);
      h_matrix1(0, 2) = 0;

      h_matrix1(1, 0) = -h_matrix1(0, 1);
      h_matrix1(1, 1) = h_matrix1(0, 0);

      h_matrix1(1, 2) = -r;

      h_matrix2(0, 0) =
          -tmp17 * (tmp14 * tmp23 + tmp15 * (tmp16 * tmp23 + tmp25));
      h_matrix2(0, 1) =
          -tmp17 * (tmp11 * (tmp13 * tmp28 + tmp25) + tmp19 * tmp28);
      h_matrix2(0, 2) =
          -tmp17 * (tmp14 * tmp29 + tmp15 * (tmp16 * tmp29 + tmp30));
      h_matrix2(0, 3) =
          -tmp17 * (tmp11 * (tmp13 * tmp31 + tmp30) + tmp19 * tmp31);
    }
  };
};
}  // namespace siconos::collision
