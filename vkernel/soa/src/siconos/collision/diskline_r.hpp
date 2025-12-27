#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/line.hpp"
// #include "siconos/model/model.hpp"

namespace siconos::collision {

struct diskline_r : item, model::relation1, model::any_lagrangian_relation {
  using attributes = gather<attribute<"line", some::item_ref<shape::line>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) line()
    {
      return make_ref_handle(self()->data(), attr<"line">(*self()));
    }

    decltype(auto) shape() { return self()->line(); }

    decltype(auto) compute_h(auto& ds)
    {
      auto& q = storage::attr<"q">(ds);
      return line().distance(q) -
             make_handle(self()->data(), prop<"shape">(ds)).radius();
    }

    void compute_jachq(auto step, auto& ds, auto& h_matrix1)
    {
      auto& data = self()->data();

      auto& q = storage::attr<"q">(ds);
      auto& r = storage::make_handle(data, storage::prop<"shape">(ds)).radius();

      auto& x = q(0);
      auto& y = q(1);
      auto& a = attr<"a">(line());
      auto& b = attr<"b">(line());
      auto& c = attr<"c">(line());
      auto& invsqrta2pb2 = attr<"invsqrta2pb2">(line());

      auto d1 = a * x + b * y + c;
      auto sd1 = copysign(1, d1);

      h_matrix1(0, 0) = a * sd1 * invsqrta2pb2;
      h_matrix1(1, 0) = -b * sd1 * invsqrta2pb2;
      h_matrix1(0, 1) = b * sd1 * invsqrta2pb2;
      h_matrix1(1, 1) = a * sd1 * invsqrta2pb2;
      h_matrix1(0, 2) = 0;
      h_matrix1(1, 2) = -r;
    }
  };
};
}  // namespace siconos::collision
