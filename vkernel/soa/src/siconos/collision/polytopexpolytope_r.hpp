#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/collision/collision_head.hpp"
#include "siconos/collision/cxd_narrowphase.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/storage/storage.hpp"

namespace siconos::collision {

/* Polytope x polytope contact relation (2D convex polygons, dof 3).

   Structural reference: diskdisk_r (item + model::relation2 +
   any_lagrangian_relation, G-matrix writing into M& args). The relation is
   batch-BACKED: it carries no compute of its own, the cxd_narrowphase
   stage item (Node template below; its own stage in Task 5) fills the
   cached status/h/normal/cp1/cp2 in one batch per step, and the two
   compute methods below read the cache/make sure the cache is fresh. An
   unbound narrowphase ref (standalone hand-cached world) no-ops.

   Templated on Topo only, so the relation <-> narrowphase mutual refs
   stay acyclic (narrowphase = cxd_narrowphase<Topo>, same Topo). */
template <typename Topo>
struct polytopexpolytope_r : item, model::relation2, model::any_lagrangian_relation {
  using dof = some::indice_parameter<"dof">;
  using narrowphase_t = cxd_narrowphase<Topo>;

  struct attributes {
    // loaded once at rebind time by the stage; unbound = stage-less
    // standalone (manually cached) relation
    some::item_ref<narrowphase_t> narrowphase;
    some::indice status;  // 0 separated, 1 enclosed (cxd::status value)
    some::scalar h;       // +distance when separated / -depth when enclosed
    some::vector<some::scalar, some::indice_value<2>> normal;  // 1 -> 2
    some::vector<some::scalar, some::indice_value<2>> cp1;
    some::vector<some::scalar, some::indice_value<2>> cp2;
    some::vector<some::scalar, some::indice_value<2>> init_dir;  // batch warm start
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) status() { return storage::attr<"status">(*self()); }
    decltype(auto) h() { return storage::attr<"h">(*self()); }
    decltype(auto) normal() { return storage::attr<"normal">(*self()); }
    decltype(auto) cp1() { return storage::attr<"cp1">(*self()); }
    decltype(auto) cp2() { return storage::attr<"cp2">(*self()); }
    decltype(auto) init_dir() { return storage::attr<"init_dir">(*self()); }

    /* ensure the stage cache is fresh for this step. The stage ref is
       bound once at rebind (Task 5); an unset ref (= the attribute's
       defaulted index 0 while the world holds no stage item to bind) is
       the standalone hand-cached case and must no-op (ledger ruling). */
    void ensure(auto step) {
      auto& ref = storage::attr<"narrowphase">(*self());
      if (ref.value() == std::remove_cvref_t<decltype(ref.value())>{0}) return;
      storage::make_ref_handle(self()->data(), ref).ensure(step);
    }

    template <match::handle<model::lagrangian_ds> DS1, match::handle<model::lagrangian_ds> DS2>
    [[gnu::always_inline]] decltype(auto) compute_h(auto step, DS1& ds1, DS2& ds2) {
      ensure(step);
      return storage::attr<"h">(*self());
    }

    /* G layout mirrors diskdisk_r (rows: normal / tangent, columns x y
       theta), hand-derived in the plan Task 4 pins:

       t = z^n (ccw perpendicular); arm_i = cp_i - c(q_i); z^arm = (-ay, ax)

       G1 = [[-n_x, -n_y, -n.(z^arm1)], [-t_x, -t_y, -t.(z^arm1)]]
       G2 = [[+n_x, +n_y, +n.(z^arm2)], [+t_x, +t_y, +t.(z^arm2)]]  */
    template <typename S, match::handle<model::lagrangian_ds> DS1,
              match::handle<model::lagrangian_ds> DS2, typename M>
    [[gnu::always_inline]] void compute_jachq(S step, DS1& ds1, DS2& ds2, M& g1, M& g2) {
      ensure(step);

      auto& q1 = storage::attr<"q">(ds1, step);
      auto& q2 = storage::attr<"q">(ds2, step);
      auto& n = storage::attr<"normal">(*self());
      auto& c1 = storage::attr<"cp1">(*self());
      auto& c2 = storage::attr<"cp2">(*self());

      const double nx = n[0];
      const double ny = n[1];
      const double tx = -ny;
      const double ty = nx;

      const double a1x = c1[0] - q1(0);
      const double a1y = c1[1] - q1(1);
      const double a2x = c2[0] - q2(0);
      const double a2y = c2[1] - q2(1);

      const double za1x = -a1y;
      const double za1y = a1x;
      const double za2x = -a2y;
      const double za2y = a2x;

      g1(0, 0) = -nx;
      g1(0, 1) = -ny;
      g1(0, 2) = -(nx * za1x + ny * za1y);
      g1(1, 0) = -tx;
      g1(1, 1) = -ty;
      g1(1, 2) = -(tx * za1x + ty * za1y);

      g2(0, 0) = nx;
      g2(0, 1) = ny;
      g2(0, 2) = nx * za2x + ny * za2y;
      g2(1, 0) = tx;
      g2(1, 1) = ty;
      g2(1, 2) = tx * za2x + ty * za2y;
    }
  };
};

}  // namespace siconos::collision
