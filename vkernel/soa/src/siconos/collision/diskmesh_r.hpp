#pragma once

#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/shape/mesh.hpp"

#include<print>

namespace siconos::collision {

struct mesh_relation {};

struct diskmesh_r : item,
                    mesh_relation,
                    model::relation2,
                    model::any_lagrangian_relation {
  struct attributes {
    some::item_ref<shape::mesh> mesh;
    some::indice contact_index;
  };

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
    decltype(auto) compute_h(auto step, auto& ds1, auto& ds2)
    {
      auto& q1 = ds1.q(step); /* disk */

      // auto& mesh = storage::prop<"shape">(ds2);
      // assert ds2 mesh == mesh()

      auto& cindex = self()->contact_index();

      return mesh().segments().distance(q1, cindex) -
             make_handle(self()->data(), prop<"shape">(ds1)).radius();
    }

    template <typename I, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2, typename M1,
              typename M2>
    void compute_jachq(I step, DS1& hds1, DS2& hds2, M1& h_matrix1,
                       M2& h_matrix2)
    {
      auto& data = self()->data();
      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = hds1.q(step);

      auto& r =
          storage::make_handle(data, storage::prop<"shape">(hds1)).radius();

      /* disk coordinates */
      const scalar& xb = q1(0);
      const scalar& yb = q1(1);

      /* segment end points */
      const scalar& xw1 = mesh().segments().x1(contact_index());
      const scalar& yw1 = mesh().segments().y1(contact_index());

      const scalar& xw2 = mesh().segments().x2(contact_index());
      const scalar& yw2 = mesh().segments().y2(contact_index());

      // auto& g1 = h_matrix1;
      // auto& g2 = h_matrix2;

      double dxs = xw2 - xw1;
      double dys = yw2 - yw1;

      double dx1 = xb - xw1;
      double dy1 = yb - yw1;

      double ds2 = dxs * dxs + dys * dys;

      double t = (dx1 * dxs + dy1 * dys) / ds2;
      t = std::min(t, 1.0);
      t = std::max(t, 0.0);

      double dx = dx1 - t * dxs;
      double dy = dy1 - t * dys;

      double nnorm = pow(dx * dx + dy * dy, 1.0 / 2.0);

      double dtx = dx / nnorm;
      double dty = dy / nnorm;

      h_matrix1(0, 0) = dtx;
      h_matrix1(0, 1) = dty;
      h_matrix1(0, 2) = 0;

      h_matrix1(1, 0) = -h_matrix1(0, 1);
      h_matrix1(1, 1) = h_matrix1(0, 0);
      h_matrix1(1, 2) = -r;

      h_matrix2(0, 0) = -(1 - t) * dtx;
      h_matrix2(0, 1) = -(1 - t) * dty;
      h_matrix2(1, 0) = -h_matrix2(0, 1);
      h_matrix2(1, 1) = h_matrix2(0, 0);

      h_matrix2(0, 2) = -t * dtx;
      h_matrix2(0, 3) = -t * dty;
      h_matrix2(1, 2) = -h_matrix2(0, 3);
      h_matrix2(1, 3) = h_matrix2(0, 2);

      std::println("rrel.contact_index()={}", contact_index());
      std::println("dtx={}", dtx);
      std::println("dty={}", dty);

      std::println("nnorm={}", nnorm);

      std::println("xw1={},yw1={},xw2={},yw2={}",xw1,yw1,xw2,yw2);
      std::println("q1={},{}", q1(0),q1(1));

    }
  };
};

}  // namespace siconos::collision

namespace siconos::storage::pattern::match
{
  template <typename T>
  concept mesh_relation = std::derived_from<T, collision::mesh_relation>;
}
