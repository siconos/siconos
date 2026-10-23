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
   The contact normal (from cxd) points body1 (polytope) -> body2 (mesh edge).
   One relation item serves every polytope near this edge, so it holds no
   per-pair state: contact() evaluates the pair and returns the result. */

template <typename Topo>
struct polytopemesh_r : item, mesh_relation, model::relation2, model::any_lagrangian_relation {
  struct attributes {
    some::item_ref<shape::mesh> mesh;
    some::indice contact_index;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mesh() {
      return storage::make_ref_handle(self()->data(), storage::attr<"mesh">(*self()));
    }
    decltype(auto) shape() { return self()->mesh(); }
    decltype(auto) contact_index() { return storage::attr<"contact_index">(*self()); }

    /* run cxd on this single polytope x mesh edge pair */
    template <match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2>
    auto contact(auto step, DS1& ds1, DS2& ds2) {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      auto& data = self()->data();

      auto& q1 = storage::attr<"q">(ds1, step);
      auto sh = collision::shape_of<shape::polytope>(ds1);
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

      return narrowphase_batch(std::move(in))[0];
    }

    template <match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2>
    [[gnu::always_inline]] auto compute_h(auto step, DS1& ds1, DS2& ds2) {
      return contact(step, ds1, ds2).h;
    }

    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::elastic_lagrangian_ds> DS2, typename M1, typename M2>
    void compute_jachq(S step, DS1& ds1, DS2& ds2, M1& g1, M2& g2) {
      const auto ct = contact(step, ds1, ds2);

      using scalar = typename decltype(self()->env())::scalar;

      auto& q1 = storage::attr<"q">(ds1, step);
      const auto& n = ct.normal;
      const auto& c1 = ct.cp1;
      const auto& c2 = ct.cp2;

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
