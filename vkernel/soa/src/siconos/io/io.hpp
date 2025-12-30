#pragma once

#include <stdexcept>

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/collision.hpp"
#include "siconos/collision/diskfdisk_r.hpp"
#include "siconos/collision/diskfsegment_r.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/storage/pattern/base.hpp"
#include "siconos/storage/pattern/pattern.hpp"
#include "siconos/storage/some/some.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::io {
using namespace storage;
using namespace storage::pattern;

template <typename Osi, typename... ContactShapes>
struct io : item {
  using osi = Osi;
  /* /!\ only system sizes defined at compile time and 2D for the moment */
  using system = nth_t<0, typename osi::systems_t>;
  using interaction = nth_t<0, typename osi::interactions_t>;
  using relations_t = typename interaction::relations;
  using contact_shapes = gather<ContactShapes...>;

  using attributes = gather<
      attribute<"osi", some::item_ref<osi>>,
      attribute<"p0_info", some::unbounded_collection<some::vector<
                               some::scalar, some::indice_value<4>>>>,
      attribute<"radii_info", some::unbounded_collection<some::vector<
                                  some::scalar, some::indice_value<2>>>>,
      attribute<"pos_info", some::unbounded_collection<some::vector<
                                some::scalar, some::indice_value<4>>>>,
      attribute<"vel_info", some::unbounded_collection<some::vector<
                                some::scalar, some::indice_value<4>>>>,
      attribute<"cp_info", some::unbounded_collection<some::vector<
                               some::scalar, some::indice_value<25>>>>,
      attribute<"co_info", some::unbounded_collection<some::vector<
                               some::scalar, some::indice_value<4>>>>,
      attribute<"work_info", some::unbounded_collection<some::vector<
                                 some::scalar, some::indice_value<5>>>>>;

  using properties =
      // an attached global ident for contact shapes
      gather<storage::attached<ContactShapes, symbol<"ident">,
                               some::integer>...>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) osi() { return storage::attr<"osi">(*self()); }

    decltype(auto) p0s(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& involveds = storage::prop_values<system, "involved">(data, step);
      auto& indices = storage::prop_values<system, "index">(data, step);

      /* /!\ first osi */
      const auto& p0_v = self()->make_handle(osi()).template visit_element<0>(
          [](auto elem) { return elem.p0_vector_assembled(); });

      attr<"p0_info">(*self()).clear();

      for (auto [id, involved, index] : view::zip(ids, involveds, indices)) {
        if (involved) {
          /* contact */
          auto&& p0 = algebra::get_vector(p0_v, index);
          attr<"p0_info">(*self()).push_back(
              {(scalar)id, p0[0], p0[1], p0[2]});
        }
        else {
          /* without contact */
          attr<"p0_info">(*self()).push_back({(scalar)id, 0., 0., 0.});
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"p0_info">(*self()).data()->data(),
          attr<"p0_info">(*self()).size(),
          attr<"p0_info">(*self()).data()->size());
    }

    decltype(auto) radii(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& shapes = storage::prop_values<system, "shape">(data, step);

      attr<"radii_info">(*self()).clear();

      for (auto [id, shape] : view::zip(ids, shapes)) {
        attr<"radii_info">(*self()).push_back(
            {(scalar)id, storage::make_handle(data, shape).radius()});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 2>>(
          attr<"radii_info">(*self()).data()->data(),
          attr<"radii_info">(*self()).size(),
          attr<"radii_info">(*self()).data()->size());
    }

    decltype(auto) positions(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& qs = storage::attr_values<system, "q">(data, step);

      attr<"pos_info">(*self()).clear();

      for (auto [id, q] : view::zip(ids, qs)) {
        attr<"pos_info">(*self()).push_back({(scalar)id, q[0], q[1], q[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"pos_info">(*self()).data()->data(),
          attr<"pos_info">(*self()).size(),
          attr<"pos_info">(*self()).data()->size());
    }

    decltype(auto) velocities(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id">(data, step);
      auto& velos = storage::attr_values<system, "velocity">(data, step);

      attr<"vel_info">(*self()).clear();

      for (auto [id, velo] : view::zip(ids, velos)) {
        attr<"vel_info">(*self()).push_back(
            {(scalar)id, velo[0], velo[1], velo[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"vel_info">(*self()).data()->data(),
          attr<"vel_info">(*self()).size(),
          attr<"vel_info">(*self()).data()->size());
    }

    decltype(auto) contact_points(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      /* /!\ first osi */
      const auto& p0_v = self()->make_handle(osi()).template visit_element<0>(
          [](auto elem) { return elem.p0_vector_assembled(); });

      auto& ys = storage::attr_values<interaction, "y">(data, step);
      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);
      auto& relations =
          storage::attr_values<interaction, "relation">(data, step);

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto& activations =
          storage::prop_values<interaction, "activation">(data, step);

      attr<"cp_info">(*self()).clear();

      indice k = 0;
      for (auto [relation, nslaw, y, ydot, lambda, ds1, ds2, activation] :
           view::zip(relations, nslaws, ys, ydots, lambdas, ds1s, ds2s,
                     activations)) {
        auto hds1 = storage::make_handle(data, ds1);
        auto hds2 = storage::make_handle(data, ds2);

        using vect = std::decay_t<decltype(hds1.q())>; /* in 2D, 3 components:
                                                          translation 2 +
                                                          orientation 1 */

        if (activation) {
          auto index_ds1 =
              prop<"index">(hds1); /* cf one_step_integrator.hpp,
                                    * assemble_h_matrix_for_involved_ds =>
                                    * row of p0_vector_assembled */
          auto index_ds2 = prop<"index">(hds2);
          auto p0 =
              algebra::get_vector(p0_v, index_ds1); /* in 2D, 2 components */

          vect c1 = {hds1.q()[0], hds1.q()[1], 0.};
          vect c2 = {hds2.q()[0], hds2.q()[1], 0.};

          vect cn;

          vect ca;

          vect cb;

          if (ds1 != ds2) {
            /* 2 disks */
            scalar dc2c1 = collision::distance(c2, c1);

            cn = (c2 - c1) / dc2c1;
            ca = c1 + storage::make_handle(data, storage::prop<"shape">(hds1))
                              .radius() *
                          cn;
            cb = c2 - storage::make_handle(data, storage::prop<"shape">(hds2))
                              .radius() *
                          cn;
          }
          else {
            variant::visit(
                data, relation,
                mp::overload(
                    /* disk / segment */
                    [&](storage::index<collision::diskfsegment_r, indice>
                            rrel) {
                      auto hrel = storage::make_handle(data, rrel);
                      /* cb is the proj point on the segment, computed a
                       * second time!
                       */
                      auto segment = hrel.segment();
                      const scalar t =
                          fmax(0, fmin(1, algebra::dot(c1 - segment.p1(),
                                                       segment.dp2p1()) /
                                              segment.length_sq()));
                      cb = segment.p1() + t * segment.dp2p1();
                      scalar dcbc1 = collision::distance(cb, c1);
                      cn = (cb - c1) / dcbc1;
                      ca = c1 + cn * storage::make_handle(
                                         data, storage::prop<"shape">(hds1))
                                         .radius();
                    },
                    /* disk / fixed disk */
                    [&](storage::index<collision::diskfdisk_r, indice> rrel) {
                      auto hrel = storage::make_handle(data, rrel);
                      auto tds = hrel.translated_disk_shape();

                      c2 = tds.translation();
                      scalar radius2 = tds.translated().radius();

                      scalar dc2c1 = collision::distance(c2, c1);

                      cn = (c2 - c1) / dc2c1;
                      ca = c1 + storage::make_handle(
                                    data, storage::prop<"shape">(hds1))
                                        .radius() *
                                    cn;
                      cb = c2 - radius2 * cn;
                    },
                    [&](auto) {
                      throw(
                          std::runtime_error("Error: not a disk-segment or "
                                             "disk-disk relation"));
                    }));
          };

          /* disk / segment */
          attr<"cp_info">(*self()).push_back(
              {storage::make_handle(data, nslaw).mu(),
               ca[0],
               ca[1],
               0. /* 2D */,
               cb[0],
               cb[1],
               0. /* 2D */,
               cn[0],
               cn[1],
               0. /* 2D */,
               p0[0],
               p0[1],
               0. /* 2D */,
               y[0],
               y[1],
               0. /* 2D */,
               ydot[0],
               ydot[1],
               0. /* 2D */,
               lambda[0],
               lambda[1],
               0. /* 2D */,
               (scalar)k,
               (scalar)index_ds1,
               (scalar)index_ds2});
          k++;
        }
      }
      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 25>>(
          attr<"cp_info">(*self()).data()->data(),
          attr<"cp_info">(*self()).size(),
          attr<"cp_info">(*self()).data()->size());
    }

    decltype(auto) contact_info(auto step)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;
      using integer = typename env_t::integer;

      auto& relations =
          storage::attr_values<interaction, "relation">(data, step);
      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto& activations =
          storage::prop_values<interaction, "activation">(data, step);

      attr<"co_info">(*self()).clear();

      indice k = 0;
      for (auto [relation, ds1, ds2, activation] :
           view::zip(relations, ds1s, ds2s, activations)) {
        if (activation) {
          auto hds1 = storage::make_handle(data, ds1);
          auto hds2 = storage::make_handle(data, ds2);
          auto index_ds1 = prop<"index">(hds1);
          auto index_ds2 = prop<"index">(hds2);

          indice inter_index = k++; /* index of interaction in indexset 1 */
          // a pair type of shape (unsigned int) + index
          auto static_shape_info = variant::visit(
              data, relation,
              mp::overload(
                  // relation1 with static shape : the concept is missing
                  [&]<match::relation1 Relation>(Relation rrel) {
                    // must provide shape method
                    return storage::prop<"ident">(rrel.shape());
                  },
                  [&](auto) {
                    // another kind of relation
                    return (integer)0;
                  }));

          attr<"co_info">(*self()).push_back(
              {(scalar)inter_index, (scalar)index_ds1, (scalar)index_ds2,
               (scalar)static_shape_info});
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"co_info">(*self()).data()->data(),
          attr<"co_info">(*self()).size(),
          attr<"co_info">(*self()).data()->size());
    }

    decltype(auto) contact_work(auto step, auto omega, auto tol)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& ydot_ks =
          storage::attr_values<interaction, "ydot">(data, step - 1);
      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);

      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);

      attr<"work_info">(*self()).clear();

      for (auto [ydot, ydot_k, lambda, index_nslaw] :
           view::zip(ydots, ydot_ks, lambdas, nslaws)) {

        auto nslaw = storage::make_handle(data, index_nslaw);

        scalar e = nslaw.e();
        scalar mu = nslaw.mu();

        // Compute normal contact work
        scalar vn_minus = ydot_k(0);
        scalar vn_plus = ydot(0);
        scalar pn = lambda(0);

        scalar vn_average = omega * vn_plus + (1. - omega) * vn_minus;
        scalar normal_contact_work = vn_average * pn;

        // Compute tangent contact work of impulse
        scalar vt_1_minus = ydot_k(1);
        scalar vt_1_plus = ydot(1);

        scalar vt_1_average = omega * vt_1_plus + (1. - omega) * vt_1_minus;
        scalar pt_1 = lambda(1);

        scalar tangent_contact_work = vt_1_average * pt_1;

        // Compute work dissipated by friction impulse
        scalar norm_vt_average = sqrt(vt_1_average * vt_1_average);
        scalar friction_dissipation = mu * norm_vt_average * pn;

        // Compute contact status
        // Warning the status are consistent for the sticking and sliding
        // only with fully implicit discretization o NewtonImpact law
        // and not wih Fremond impact law
        scalar norm_pt = sqrt(pt_1 * pt_1);
        scalar norm_vt_plus = sqrt(vt_1_plus * vt_1_plus);

        scalar answer_4;
        scalar answer_5;
        if ((pn < tol) and (vn_plus + e * vn_minus > tol)) {
          answer_4 = 0;  // take-off = 0
        }
        else if ((pn > tol) and (vn_plus + e * vn_minus > tol)) {
          if ((norm_pt - mu * pn > tol)) {
            answer_4 = -3;  // outside the cone
          }
          else if (norm_pt - mu * pn < -tol) {
            if (norm_vt_plus > tol) {
              answer_4 =
                  -2;  // sticking with a non zero slifing velocity = -2
            }
            else {
              answer_4 = 1;  // sticking
            }
          }
          else {
            answer_4 = 2;
          }  // sliding
        }
        else {
          answer_4 = -1;
        }  // undetermined
        if ((pn > tol) and (vn_minus > tol)) {
          answer_5 = normal_contact_work;
        }

        attr<"work_info">(*self()).push_back(
            {normal_contact_work, tangent_contact_work, friction_dissipation,
             answer_4, answer_5});
      }
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("p0s", &interface<Handle>::p0s<indice>),
          method("radii", &interface<Handle>::radii<indice>),
          method("positions", &interface<Handle>::positions<indice>),
          method("velocities", &interface<Handle>::velocities<indice>),
          method("contact_points",
                 &interface<Handle>::contact_points<indice>),
          method("contact_info", &interface<Handle>::contact_info<indice>),
          method("contact_work",
                 &interface<Handle>::contact_work<indice, scalar, scalar>));
    }
  };
};
}  // namespace siconos::io
