#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/cxd_narrowphase.hpp"
#include "siconos/collision/shape/mesh.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Polytope x mesh-edge contact relation (2D convex polygon vs one mesh
   segment, per (mesh, contact_index)). The mesh edge is treated as a 2-vertex
   polytope whose world vertices are the segment endpoints (re-packed each
   call, since meshes may deform for FEM). Analytic-ish: geometry is evaluated
   per pair via a one-pair cxd batch, so this relation does not depend on the
   cxd_narrowphase stage. Consistent with diskmesh_r / polytopefsegment_r.
   The contact normal (from cxd) points body1 (polytope) -> body2 (mesh edge). */

template <typename Topo>
struct polytopemesh_r : item, mesh_relation, model::relation2, model::any_lagrangian_relation {
  struct attributes {
    some::item_ref<shape::mesh> mesh;
    some::indice contact_index;
    some::indice status;
    some::scalar h;
    some::vector<some::scalar, some::indice_value<2>> normal;
    some::vector<some::scalar, some::indice_value<2>> cp1;
    some::vector<some::scalar, some::indice_value<2>> cp2;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mesh() {
      return storage::make_ref_handle(self()->data(), storage::attr<"mesh">(*self()));
    }
    decltype(auto) shape() { return self()->mesh(); }
    decltype(auto) contact_index() { return storage::attr<"contact_index">(*self()); }
    decltype(auto) status() { return storage::attr<"status">(*self()); }
    decltype(auto) h() { return storage::attr<"h">(*self()); }
    decltype(auto) normal() { return storage::attr<"normal">(*self()); }
    decltype(auto) cp1() { return storage::attr<"cp1">(*self()); }
    decltype(auto) cp2() { return storage::attr<"cp2">(*self()); }

    /* run cxd on this single polytope x segment pair and cache the result */
    template <match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2>
    void evaluate(auto step, DS1& ds1, DS2& ds2) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();

      auto& q1 = storage::attr<"q">(ds1, step);
      auto sh =
          storage::make_handle(data, storage::prop<"shape", storage::config_declared>(ds1));
      const scalar c = std::cos(q1(2)), s = std::sin(q1(2));
      const auto n = sh.nvertices();

      narrowphase_input<scalar> in;
      // entry 0: polytope (world vertices)
      const std::uint32_t off0 = 0;
      for (std::size_t i = 0; i < n; ++i) {
        const auto& v = sh.vertices()[i];
        in.vertices.push_back(
            std::array<scalar, 2>{q1(0) + c * v[0] - s * v[1], q1(1) + s * v[0] + c * v[1]});
      }
      in.entries.push_back({off0, static_cast<std::uint32_t>(n)});
      const std::uint32_t slot0 = 0;
      // entry 1: mesh edge endpoints (static world coords within this step)
      auto mesh = self()->mesh();
      auto ci = self()->contact_index();
      const std::uint32_t off1 = static_cast<std::uint32_t>(in.vertices.size());
      in.vertices.push_back(
          std::array<scalar, 2>{mesh.segments().x1(ci), mesh.segments().y1(ci)});
      in.vertices.push_back(
          std::array<scalar, 2>{mesh.segments().x2(ci), mesh.segments().y2(ci)});
      in.entries.push_back({off1, 2});
      const std::uint32_t slot1 = 1;

      in.ids_a.push_back(slot0);
      in.ids_b.push_back(slot1);
      in.init_dirs.push_back({1.0, 0.0});

      auto out = narrowphase_batch(std::move(in));
      const auto& res = out[0];
      storage::attr<"status">(*self()) = res.status;
      storage::attr<"h">(*self()) = res.h;
      storage::attr<"normal">(*self()) = {res.normal[0], res.normal[1]};
      storage::attr<"cp1">(*self()) = {res.cp1[0], res.cp1[1]};
      storage::attr<"cp2">(*self()) = {res.cp2[0], res.cp2[1]};
    }

    template <match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2>
    [[gnu::always_inline]] decltype(auto) compute_h(auto step, DS1& ds1, DS2& ds2) {
      evaluate(step, ds1, ds2);
      return storage::attr<"h">(*self());
    }

    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2, typename M1, typename M2>
    void compute_jachq(S step, DS1& ds1, DS2& ds2, M1& g1, M2& g2) {
      evaluate(step, ds1, ds2);

      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = storage::attr<"q">(ds1, step);
      auto& n = storage::attr<"normal">(*self());
      auto& c1 = storage::attr<"cp1">(*self());
      auto& c2 = storage::attr<"cp2">(*self());

      const scalar nx = n[0];
      const scalar ny = n[1];
      const scalar tx = -ny;
      const scalar ty = nx;

      const scalar a1x = c1[0] - q1(0);
      const scalar a1y = c1[1] - q1(1);
      const scalar za1x = -a1y;
      const scalar za1y = a1x;

      /* polytope side (body 1): normal/tangent rows, theta column via lever arm */
      g1(0, 0) = -nx;
      g1(0, 1) = -ny;
      g1(0, 2) = -(nx * za1x + ny * za1y);
      g1(1, 0) = -tx;
      g1(1, 1) = -ty;
      g1(1, 2) = -(tx * za1x + ty * za1y);

      /* mesh edge side (body 2): endpoints x1,y1,x2,y2, weighted by (1-t), t */
      auto mesh = self()->mesh();
      auto ci = self()->contact_index();
      const scalar xw1 = mesh.segments().x1(ci);
      const scalar yw1 = mesh.segments().y1(ci);
      const scalar xw2 = mesh.segments().x2(ci);
      const scalar yw2 = mesh.segments().y2(ci);

      const scalar dxs = xw2 - xw1;
      const scalar dys = yw2 - yw1;
      const scalar seglen2 = dxs * dxs + dys * dys;
      scalar t = 0;
      if (seglen2 > scalar{0}) {
        t = ((c2[0] - xw1) * dxs + (c2[1] - yw1) * dys) / seglen2;
        t = std::min(std::max(t, scalar{0}), scalar{1});
      }

      /* normal row (1 -> 2 direction): mesh pushed toward polytope = +normal */
      g2(0, 0) = (1 - t) * nx;
      g2(0, 1) = (1 - t) * ny;
      g2(0, 2) = t * nx;
      g2(0, 3) = t * ny;
      /* tangent row: CCW perpendicular */
      g2(1, 0) = -(1 - t) * ny;
      g2(1, 1) = (1 - t) * nx;
      g2(1, 2) = -t * ny;
      g2(1, 3) = t * nx;
    }
  };
};

}  // namespace siconos::collision
