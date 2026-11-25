// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include <array>
#include <functional>
#include <map>
#include <stdexcept>

#include "siconos/vkernel/algebra/numerics_backend.hpp"
#include "siconos/vkernel/collision/collision.hpp"
#include "siconos/vkernel/collision/diskfdisk_r.hpp"
#include "siconos/vkernel/collision/diskfsegment_r.hpp"
#include "siconos/vkernel/collision/diskmesh_r.hpp"
#include "siconos/vkernel/collision/polytopedisk_r.hpp"
#include "siconos/vkernel/collision/polytopefdisk_r.hpp"
#include "siconos/vkernel/collision/polytopefsegment_r.hpp"
#include "siconos/vkernel/collision/polytopemesh_r.hpp"
#include "siconos/vkernel/collision/polytopepolytope_r.hpp"
#include "siconos/vkernel/collision/shape_of.hpp"
#include "siconos/vkernel/model/lagrangian_r.hpp"
#include "siconos/vkernel/model/nslaws.hpp"
#include "siconos/vkernel/storage/handle.hpp"
#include "siconos/vkernel/storage/pattern/base.hpp"
#include "siconos/vkernel/storage/pattern/pattern.hpp"
#include "siconos/vkernel/storage/some/some.hpp"
#include "siconos/vkernel/storage/storage.hpp"
#include "siconos/vkernel/utils/variant.hpp"

namespace siconos::vkernel::io {
using namespace storage;
using namespace storage::pattern;

template <typename Osi, typename... ContactShapes>
struct io : item {
  using osi = Osi;
  // declares the attached storages of systems and interactions (id, ds1, ...)
  using topology = typename osi::topology;
  /* /!\ only system sizes defined at compile time and 2D for the moment */
  using system = typename topology::fixed_dof_system;     // fixed dof system
  using dsystem = typename topology::dynamic_dof_system;  // runtime dof system
                                                          // maybe empty_item
  using interaction = typename topology::fixed_dof_interaction;
  // disk / fem interaction (maybe empty_item)
  using dinteraction = typename topology::dynamic_dof_fixed_dof_interaction;
  using relations_t = typename interaction::relations;
  using contact_shapes = gather<ContactShapes...>;

  // the columns of a contact (contact_points): 25, and with a law with
  // rolling resistance two 3-vectors more, the rolling rows of ydot and
  // lambda (the relative angular velocity and the rolling moment), the
  // second and third components zero in 2D
  static constexpr bool with_rolling_rows =
      std::derived_from<typename interaction::nslaw, model::newton_impact_rolling_friction>;
  static constexpr int contact_columns = with_rolling_rows ? 31 : 25;

  using attributes = gather<
      attribute<"osi", some::item_ref<osi>>,
      attribute<"disp_info", some::map<some::indice, some::unbounded_vector<some::scalar>>>,
      attribute<"p0_info",
                some::unbounded_collection<some::vector<some::scalar, some::indice_value<4>>>>,
      attribute<"radii_info",
                some::unbounded_collection<some::vector<some::scalar, some::indice_value<2>>>>,
      attribute<"pos_info",
                some::unbounded_collection<some::vector<some::scalar, some::indice_value<4>>>>,
      attribute<"vel_info",
                some::unbounded_collection<some::vector<some::scalar, some::indice_value<4>>>>,
      attribute<"cp_info", some::unbounded_collection<some::vector<
                               some::scalar, some::indice_value<contact_columns>>>>,
      attribute<"co_info",
                some::unbounded_collection<some::vector<some::scalar, some::indice_value<4>>>>,
      attribute<"work_info", some::unbounded_collection<
                                 some::vector<some::scalar, some::indice_value<7>>>>>;

  using properties =
      // an attached global ident for contact shapes
      gather<storage::attached<ContactShapes, symbol<"ident">, some::integer>...>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) osi() { return storage::attr<"osi">(*self()); }

    decltype(auto) p0s(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);
      auto& involveds = storage::prop_values<system, "involved", topology>(data, step);
      auto& indices = storage::prop_values<system, "index", topology>(data, step);

      /* /!\ the rigid systems, the first kind */
      auto hosi = self()->make_handle(osi());
      auto&& p0h = hosi.p0_vector_assembled();
      // p0 is assembled only when there are active contacts (the solve block):
      // a predicted (involved) contact can be output at a step where the solve
      // did not run, and then p0 has no storage yet.  Do not dereference it.
      const bool p0_ready = (p0h._v != nullptr);
      const auto& p0_v = std::get<0>(hosi.systems()).p0_view(p0h);

      attr<"p0_info">(*self()).clear();

      for (auto [id, involved, index] : view::zip(ids, involveds, indices)) {
        // p0 is sized by the solve (impulse) to the involved systems: an
        // involved system is read only when p0 exists and its block is in
        // range, so a scheme that reaches the output before the solve has
        // sized p0 (or an empty p0) does not read past the assembled vector
        if (involved and p0_ready and std::size_t(index) < algebra::size0(p0_v)) {
          /* contact */
          auto&& p0 = algebra::get_vector(p0_v, index);
          attr<"p0_info">(*self()).push_back({(scalar)id, p0[0], p0[1], p0[2]});
        } else {
          /* without contact */
          attr<"p0_info">(*self()).push_back({(scalar)id, 0., 0., 0.});
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"p0_info">(*self()).data()->data(), attr<"p0_info">(*self()).size(),
          attr<"p0_info">(*self()).data()->size());
    }

    decltype(auto) radii(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);

      attr<"radii_info">(*self()).clear();

      // a disk's radius; a polytope's circumradius (vview scales the disk
      // glyphs with it, polytopes are drawn from their vertices)
      for (auto [id, h] : view::zip(ids, storage::handles<system>(data, step))) {
        using h_t = std::decay_t<decltype(h)>;
        scalar r = 0;
        if constexpr (collision::may_have_shape<collision::shape::disk, h_t>()) {
          if (collision::has_shape<collision::shape::disk>(h))
            r = collision::shape_of<collision::shape::disk>(h).radius();
        }
        if constexpr (collision::may_have_shape<collision::shape::polytope, h_t>()) {
          if (collision::has_shape<collision::shape::polytope>(h))
            r = collision::shape_of<collision::shape::polytope>(h).rmax();
        }
        attr<"radii_info">(*self()).push_back({(scalar)id, r});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 2>>(
          attr<"radii_info">(*self()).data()->data(), attr<"radii_info">(*self()).size(),
          attr<"radii_info">(*self()).data()->size());
    }

    decltype(auto) displacements(auto step) {
      auto& data = self()->data();
      if constexpr (!match::empty_item<dsystem>) {
        auto& ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& qs = storage::attr_values<dsystem, "q">(data, step);
        auto& q0s = storage::prop_values<dsystem, "q0", topology>(data, step);

        attr<"disp_info">(*self()).clear();

        for (auto [id, q, q0] : view::zip(ids, qs, q0s)) {
          attr<"disp_info">(*self())[id] = q - q0;
        }

        return attr<"disp_info">(*self());
      }
    }

    auto stiffness_matrices(auto step) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;
      using stiffness_matrices_t =
          std::map<indice, std::reference_wrapper<algebra::sparse_matrix<scalar>>>;
      stiffness_matrices_t result;

      if constexpr (!match::empty_item<dsystem>) {
        auto& data = self()->data();
        auto& ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& matrices = storage::attr_values<dsystem, "k_matrix">(data, step);
        for (auto [id, matrix] : view::zip(ids, matrices))
          result.emplace(id, std::ref(matrix));
      }

      return result;
    }

    // the mass matrix of each FEM system, by reference, as stiffness_matrices
    // (the consistent mass of the FEM model, constant)
    auto mass_matrices(auto step) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;
      using mass_matrices_t =
          std::map<indice, std::reference_wrapper<algebra::sparse_matrix<scalar>>>;
      mass_matrices_t result;

      if constexpr (!match::empty_item<dsystem>) {
        auto& data = self()->data();
        auto& ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& matrices = storage::attr_values<dsystem, "mass_matrix">(data, step);
        for (auto [id, matrix] : view::zip(ids, matrices))
          result.emplace(id, std::ref(matrix));
      }

      return result;
    }

    // the velocity of each FEM system at step (all its dofs), keyed by id as
    // displacements
    auto fem_velocities(auto step) {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;
      std::map<indice, algebra::unbounded_vector<scalar>> result;

      if constexpr (!match::empty_item<dsystem>) {
        auto& data = self()->data();
        auto& ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& velocities = storage::attr_values<dsystem, "velocity">(data, step);
        for (auto [id, v] : view::zip(ids, velocities)) result.emplace(id, v);
      }

      return result;
    }

    // the kinetic energy 1/2 v^T M v of each system at step, (id, K): the
    // rigid systems (M diagonal, m, m and the inertia: the rotation is
    // included), then the FEM systems (the sparse M)
    auto kinetic_energies(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using rows_t = algebra::unbounded_col_matrix<scalar, 2>;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);
      auto& masses = storage::attr_values<system, "mass_matrix">(data, step);
      auto& velocities = storage::attr_values<system, "velocity">(data, step);
      std::size_t n = ids.size();
      if constexpr (!match::empty_item<dsystem>) {
        n += storage::prop_values<dsystem, "id", topology>(data, step).size();
      }
      rows_t result(n, 2);

      std::size_t r = 0;
      for (auto [id, m, v] : view::zip(ids, masses, velocities)) {
        result(r, 0) = id;
        result(r++, 1) = 0.5 * v.dot(m * v);
      }

      if constexpr (!match::empty_item<dsystem>) {
        auto& d_ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& d_masses = storage::attr_values<dsystem, "mass_matrix">(data, step);
        auto& d_velocities = storage::attr_values<dsystem, "velocity">(data, step);
        for (auto [id, m, v] : view::zip(d_ids, d_masses, d_velocities)) {
          algebra::unbounded_vector<scalar> mv = m * v;
          result(r, 0) = id;
          result(r++, 1) = 0.5 * v.dot(mv);
        }
      }

      return result;
    }

    // the work of the forces over the step that ends at step, (id, W), as
    // MoreauJeanOSI::computeWorkForces of the AOS kernel:
    //   W = h f_theta . v_theta,  x_theta = theta x(step) + (1 - theta) x(step - 1)
    // f the total forces: fext for a rigid system (constant over the step),
    // fext - K q for a FEM system (the elastic force included). theta and h
    // are those of the integrator. The rigid systems, then the FEM systems.
    auto work_of_forces(auto step, auto theta, auto h) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;
      using rows_t = algebra::unbounded_col_matrix<scalar, 2>;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);
      auto& fexts = storage::attr_values<system, "fext">(data, step);
      auto& velocities = storage::attr_values<system, "velocity">(data, step);
      auto& velocities_k = storage::attr_values<system, "velocity">(data, step - 1);
      std::size_t n = ids.size();
      if constexpr (!match::empty_item<dsystem>) {
        n += storage::prop_values<dsystem, "id", topology>(data, step).size();
      }
      rows_t result(n, 2);

      std::size_t r = 0;
      for (auto [id, f, v, v_k] : view::zip(ids, fexts, velocities, velocities_k)) {
        result(r, 0) = id;
        result(r++, 1) = h * f.dot(theta * v + (1 - theta) * v_k);
      }

      if constexpr (!match::empty_item<dsystem>) {
        auto& d_ids = storage::prop_values<dsystem, "id", topology>(data, step);
        auto& d_fexts = storage::attr_values<dsystem, "fext">(data, step);
        auto& d_ks = storage::attr_values<dsystem, "k_matrix">(data, step);
        auto& d_qs = storage::attr_values<dsystem, "q">(data, step);
        auto& d_qs_k = storage::attr_values<dsystem, "q">(data, step - 1);
        auto& d_velocities = storage::attr_values<dsystem, "velocity">(data, step);
        auto& d_velocities_k = storage::attr_values<dsystem, "velocity">(data, step - 1);
        for (auto [id, f, k, q, q_k, v, v_k] :
             view::zip(d_ids, d_fexts, d_ks, d_qs, d_qs_k, d_velocities, d_velocities_k)) {
          algebra::unbounded_vector<scalar> q_theta = theta * q + (1 - theta) * q_k;
          algebra::unbounded_vector<scalar> f_theta = f - k * q_theta;
          algebra::unbounded_vector<scalar> v_theta = theta * v + (1 - theta) * v_k;
          result(r, 0) = id;
          result(r++, 1) = h * f_theta.dot(v_theta);
        }
      }

      return result;
    }

    decltype(auto) positions(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);
      auto& qs = storage::attr_values<system, "q">(data, step);

      attr<"pos_info">(*self()).clear();

      for (auto [id, q] : view::zip(ids, qs)) {
        attr<"pos_info">(*self()).push_back({(scalar)id, q[0], q[1], q[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"pos_info">(*self()).data()->data(), attr<"pos_info">(*self()).size(),
          attr<"pos_info">(*self()).data()->size());
    }

    decltype(auto) velocities(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      //      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ids = storage::prop_values<system, "id", topology>(data, step);
      auto& velos = storage::attr_values<system, "velocity">(data, step);

      attr<"vel_info">(*self()).clear();

      for (auto [id, velo] : view::zip(ids, velos)) {
        attr<"vel_info">(*self()).push_back({(scalar)id, velo[0], velo[1], velo[2]});
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"vel_info">(*self()).data()->data(), attr<"vel_info">(*self()).size(),
          attr<"vel_info">(*self()).data()->size());
    }

    // a row of contact_points: its 25 columns, then the rolling rows of ydot
    // and lambda with a law with rolling resistance
    void add_contact_point(const auto& columns, const auto& ydot, const auto& lambda) {
      auto& rows = attr<"cp_info">(*self());
      rows.emplace_back();
      auto& row = rows.back();
      for (std::size_t i = 0; i < columns.size(); ++i) row[i] = columns[i];
      if constexpr (with_rolling_rows) {
        row[25] = ydot[2];
        row[26] = 0.;
        row[27] = 0.;
        row[28] = lambda[2];
        row[29] = 0.;
        row[30] = 0.;
      }
    }

    decltype(auto) contact_points(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      auto& ys = storage::attr_values<interaction, "y">(data, step);
      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);
      auto& relations = storage::attr_values<interaction, "relation">(data, step);

      auto& ds1s = storage::prop_values<interaction, "ds1", topology>(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2", topology>(data, step);

      auto& h_matrices1 = storage::attr_values<interaction, "h_matrix1">(data, step);

      auto& activations =
          storage::prop_values<interaction, "activation", topology>(data, step);

      attr<"cp_info">(*self()).clear();

      indice k = 0;
      for (auto [relation, h_matrix1, nslaw, y, ydot, lambda, ds1, ds2, activation] :
           view::zip(relations, h_matrices1, nslaws, ys, ydots, lambdas, ds1s, ds2s,
                     activations)) {
        auto hds1 = storage::make_handle(data, ds1);
        auto hds2 = storage::make_handle(data, ds2);

        using vect = std::decay_t<decltype(hds1.q(step))>; /* in 2D, 3 components:
                                                          translation 2 +
                                                          orientation 1 */

        if (activation) {
          // computation at the interaction level
          auto p0 = h_matrix1.transpose() * lambda;

          vect c1 = {hds1.q(step)[0], hds1.q(step)[1], 0.};
          vect c2 = {hds2.q(step)[0], hds2.q(step)[1], 0.};

          vect cn;

          vect ca;

          vect cb;

          // a contact given by its two points and its normal (1 -> 2)
          auto from_contact = [&](const auto& ct) {
            ca = {ct.cp1[0], ct.cp1[1], 0.};
            cb = {ct.cp2[0], ct.cp2[1], 0.};
            cn = {ct.normal[0], ct.normal[1], 0.};
          };

          if (ds1 != ds2) {
            variant::visit(
                data, relation,
                mp::overload(
                    /* 2 disks */
                    [&](match::handle<collision::diskdisk_r> auto) {
                      scalar dc2c1 = collision::distance(c2, c1);
                      cn = (c2 - c1) / dc2c1;
                      ca = c1 + collision::disk_radius(hds1) * cn;
                      cb = c2 - collision::disk_radius(hds2) * cn;
                    },
                    /* 2 polytopes: the stage's cached manifold point */
                    [&](match::handle<collision::polytopepolytope_r<topology>> auto rrel) {
                      ca = {rrel.cp1()[0], rrel.cp1()[1], 0.};
                      cb = {rrel.cp2()[0], rrel.cp2()[1], 0.};
                      cn = {rrel.normal()[0], rrel.normal()[1], 0.};
                    },
                    /* disk / polytope */
                    [&](match::handle<collision::polytopedisk_r<topology>> auto rrel) {
                      from_contact(rrel.contact(step, hds1, hds2));
                    },
                    [&](auto) {
                      throw(std::runtime_error("io::contact_points: unknown ds-ds relation"));
                    }));
          } else {
            variant::visit(
                data, relation,
                mp::overload(
                    /* disk / segment */
                    [&](storage::index<collision::diskfsegment_r, indice> rrel) {
                      auto hrel = storage::make_handle(data, rrel);
                      /* cb is the proj point on the segment, computed a
                       * second time!
                       */
                      auto segment = hrel.segment();
                      const scalar t =
                          fmax(0, fmin(1, algebra::dot(c1 - segment.p1(), segment.dp2p1()) /
                                              segment.length_sq()));
                      cb = segment.p1() + t * segment.dp2p1();
                      scalar dcbc1 = collision::distance(cb, c1);
                      cn = (cb - c1) / dcbc1;
                      ca = c1 + cn * collision::disk_radius(hds1);
                    },
                    /* polytope / fixed segment, polytope / fixed disk */
                    [&](match::handle<collision::polytopefsegment_r<topology>> auto rrel) {
                      from_contact(rrel.contact(step, hds1));
                    },
                    [&](match::handle<collision::polytopefdisk_r<topology>> auto rrel) {
                      from_contact(rrel.contact(step, hds1));
                    },
                    /* disk / fixed disk */
                    [&](storage::index<collision::diskfdisk_r, indice> rrel) {
                      auto hrel = storage::make_handle(data, rrel);
                      auto tds = hrel.translated_disk_shape();

                      c2 = tds.translation();
                      scalar radius2 = tds.translated().radius();

                      scalar dc2c1 = collision::distance(c2, c1);

                      cn = (c2 - c1) / dc2c1;
                      ca = c1 + collision::disk_radius(hds1) * cn;
                      cb = c2 - radius2 * cn;
                    },
                    [&](auto) {
                      throw(
                          std::runtime_error("Error: not a disk-segment or "
                                             "disk-disk relation"));
                    }));
          };

          auto id_ds1 = storage::prop<"id", topology>(hds1);
          auto id_ds2 = storage::prop<"id", topology>(hds2);

          /* disk / segment */
          add_contact_point(std::array<scalar, 25>{storage::make_handle(data, nslaw).mu(),
                                                   cb[0],
                                                   cb[1],
                                                   0. /* 2D */,
                                                   ca[0],
                                                   ca[1],
                                                   0. /* 2D */,
                                                   -cn[0],
                                                   -cn[1],
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
                                                   (scalar)id_ds1,
                                                   (scalar)id_ds2},
                            ydot, lambda);
          k++;
        }
      }

      // disk / fem contacts, numbered after the others as in the assembled
      // problem
      if constexpr (!match::empty_item<dinteraction>) {
        auto& d_ys = storage::attr_values<dinteraction, "y">(data, step);
        auto& d_ydots = storage::attr_values<dinteraction, "ydot">(data, step);
        auto& d_lambdas = storage::attr_values<dinteraction, "lambda">(data, step);
        auto& d_nslaws = storage::attr_values<dinteraction, "nslaw">(data, step);
        auto& d_relations = storage::attr_values<dinteraction, "relation">(data, step);
        auto& d_ds1s = storage::prop_values<dinteraction, "ds1", topology>(data, step);
        auto& d_ds2s = storage::prop_values<dinteraction, "ds2", topology>(data, step);
        auto& d_h_matrices1 = storage::attr_values<dinteraction, "h_matrix1">(data, step);
        auto& d_activations =
            storage::prop_values<dinteraction, "activation", topology>(data, step);

        for (auto [relation, h_matrix1, nslaw, y, ydot, lambda, ds1, ds2, activation] :
             view::zip(d_relations, d_h_matrices1, d_nslaws, d_ys, d_ydots, d_lambdas, d_ds1s,
                       d_ds2s, d_activations)) {
          if (activation) {
            auto hds1 = storage::make_handle(data, ds1); /* disk */
            auto hds2 = storage::make_handle(data, ds2); /* fem */

            using vect = std::decay_t<decltype(hds1.q(step))>;

            // impulse on the disk
            auto p0 = h_matrix1.transpose() * lambda;

            vect c1 = {hds1.q(step)[0], hds1.q(step)[1], 0.};
            vect ca;
            vect cb;
            vect cn;

            variant::visit(
                data, relation,
                mp::overload(
                    /* disk / mesh segment: the segment of the mesh in its
                     * reference configuration, as in diskmesh_r */
                    [&](match::handle<collision::diskmesh_r> auto rrel) {
                      auto segments = rrel.mesh().segments();
                      auto ci = rrel.contact_index();
                      vect p1 = {segments.x1(ci), segments.y1(ci), 0.};
                      vect p2 = {segments.x2(ci), segments.y2(ci), 0.};
                      vect dp = p2 - p1;
                      vect d1 = c1 - p1;
                      const scalar t =
                          fmax(0, fmin(1, algebra::dot(d1, dp) / algebra::dot(dp, dp)));
                      cb = p1 + t * dp;
                      cn = (cb - c1) / collision::distance(cb, c1);
                      ca = c1 + cn * collision::disk_radius(hds1);
                    },
                    /* polytope / mesh segment */
                    [&](match::handle<collision::polytopemesh_r<topology>> auto rrel) {
                      const auto ct = rrel.contact(step, hds1, hds2);
                      ca = {ct.cp1[0], ct.cp1[1], 0.};
                      cb = {ct.cp2[0], ct.cp2[1], 0.};
                      cn = {ct.normal[0], ct.normal[1], 0.};
                    },
                    [&](auto) {
                      throw(std::runtime_error("Error: not a disk-mesh relation"));
                    }));

            auto id_ds1 = storage::prop<"id", topology>(hds1);
            auto id_ds2 = storage::prop<"id", topology>(hds2);

            add_contact_point(std::array<scalar, 25>{storage::make_handle(data, nslaw).mu(),
                                                     cb[0],
                                                     cb[1],
                                                     0. /* 2D */,
                                                     ca[0],
                                                     ca[1],
                                                     0. /* 2D */,
                                                     -cn[0],
                                                     -cn[1],
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
                                                     (scalar)id_ds1,
                                                     (scalar)id_ds2},
                              ydot, lambda);
            k++;
          }
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, contact_columns>>(
          attr<"cp_info">(*self()).data()->data(), attr<"cp_info">(*self()).size(),
          attr<"cp_info">(*self()).data()->size());
    }

    decltype(auto) contact_info(auto step) {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;
      using integer = typename env_t::integer;

      auto& relations = storage::attr_values<interaction, "relation">(data, step);
      auto& ds1s = storage::prop_values<interaction, "ds1", topology>(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2", topology>(data, step);
      auto& activations =
          storage::prop_values<interaction, "activation", topology>(data, step);

      attr<"co_info">(*self()).clear();

      indice k = 0;
      for (auto [relation, ds1, ds2, activation] :
           view::zip(relations, ds1s, ds2s, activations)) {
        if (activation) {
          auto hds1 = storage::make_handle(data, ds1);
          auto hds2 = storage::make_handle(data, ds2);
          auto id_ds1 = prop<"id", topology>(hds1);
          auto id_ds2 = prop<"id", topology>(hds2);

          indice inter_index = k++; /* index of interaction in indexset 1 */
          // a pair type of shape (unsigned int) + index
          auto static_shape_info =
              variant::visit(data, relation,
                             mp::overload(
                                 // relation1 with static shape : the concept is missing
                                 [&]<match::relation1 Relation>(Relation rrel) {
                                   // must provide shape method
                                   return storage::prop<"ident", io>(rrel.shape());
                                 },
                                 [&](auto) {
                                   // another kind of relation
                                   return (integer)0;
                                 }));

          attr<"co_info">(*self()).push_back({(scalar)inter_index, (scalar)id_ds1,
                                              (scalar)id_ds2, (scalar)static_shape_info});
        }
      }

      // disk / fem contacts, numbered after the others as in contact_points
      if constexpr (!match::empty_item<dinteraction>) {
        auto& d_ds1s = storage::prop_values<dinteraction, "ds1", topology>(data, step);
        auto& d_ds2s = storage::prop_values<dinteraction, "ds2", topology>(data, step);
        auto& d_activations =
            storage::prop_values<dinteraction, "activation", topology>(data, step);

        for (auto [ds1, ds2, activation] : view::zip(d_ds1s, d_ds2s, d_activations)) {
          if (activation) {
            auto id_ds1 = prop<"id", topology>(storage::make_handle(data, ds1));
            auto id_ds2 = prop<"id", topology>(storage::make_handle(data, ds2));
            attr<"co_info">(*self()).push_back(
                {(scalar)k++, (scalar)id_ds1, (scalar)id_ds2, (scalar)0});
          }
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 4>>(
          attr<"co_info">(*self()).data()->data(), attr<"co_info">(*self()).size(),
          attr<"co_info">(*self()).data()->size());
    }

    decltype(auto) contact_work(auto step, auto theta, auto tol, auto gap_law, auto h) {
      // fix: theta is omega in mechanics_run.py & MechanicsIO.hpp
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      attr<"work_info">(*self()).clear();

      // active interactions only (index set 1), numbered as in contact_info:
      // the rigid contacts, then the disk / fem contacts
      scalar k = 0;
      // gap: the inelastic gap laws apply to the rigid contacts only, the disk /
      // fem contacts keep the Newton impact law
      auto add_contact_work = [&](const auto& ydot, const auto& ydot_k, const auto& y,
                                  const auto& lambda, const auto& index_nslaw, bool gap) {
        handle nslaw = storage::make_handle(data, index_nslaw);

        scalar e = nslaw.e();
        scalar mu = nslaw.mu();

        // index of the interaction in index set 1: column 0 of contact_info
        // gives the body ids
        scalar answer_0 = k++;

        // Compute normal contact work
        scalar vn_minus = ydot_k(0);
        scalar vn_plus = ydot(0);
        scalar pn = lambda(0);

        scalar vn_average = 1 / 2. * (vn_plus + vn_minus);
        scalar normal_contact_work = vn_average * pn;
        scalar answer_1 = normal_contact_work;

        scalar vn_average_theta = theta * vn_plus + (1. - theta) * vn_minus;
        scalar normal_contact_work_theta = vn_average_theta * pn;
        scalar answer_3 = normal_contact_work_theta;

        // Compute tangent contact work of impulse
        scalar vt_1_minus = ydot_k(1);
        scalar vt_1_plus = ydot(1);
        scalar pt_1 = lambda(1);
        ;

        scalar vt_1_average = 1 / 2. * (vt_1_plus + vt_1_minus);
        scalar tangent_contact_work = vt_1_average * pt_1;
        scalar answer_2 = tangent_contact_work;

        scalar vt_1_average_theta = theta * vt_1_plus + (1. - theta) * vt_1_minus;
        scalar tangent_contact_work_theta = vt_1_average_theta * pt_1;
        scalar answer_4 = tangent_contact_work_theta;

        // // Compute directly work dissipated by friction impulse
        // scalar norm_vt_average = sqrt(vt_1_average * vt_1_average);
        // scalar friction_dissipation = mu * norm_vt_average * pn;
        // answer(3) = friction_dissipation;

        /* Compute contact status
         * Warning the status are consistent for the sticking and sliding
         * only with fully implicit discretization o NewtonImpact law;
         * with the Fremond impact law friction acts on the theta-weighted
         * tangential velocity
         */

        scalar answer_5;
        scalar answer_6 = 0.;  // work of the impact law, 0 without impact

        // the normal complementarity residual at step + 1: the Newton impact law
        // vn_plus + e vn_minus, or, for the inelastic gap laws (quasi inelastic and
        // the Anitescu-Tasora stabilisation), vn_plus + g_pos / h
        const scalar normal_residual =
            vn_plus + (gap ? std::max(scalar(0), y(0)) / h : e * vn_minus);

        scalar norm_pt = sqrt(pt_1 * pt_1);
        scalar norm_vt_plus = sqrt(vt_1_plus * vt_1_plus);
        if constexpr (std::derived_from<typename interaction::nslaw,
                                        model::fremond_impact_friction>) {
          norm_vt_plus = sqrt(vt_1_average_theta * vt_1_average_theta);
        }
        // scalar norm_vt_minus = sqrt(vt_1_minus*vt_1_minus);
        if ((pn < tol) and (normal_residual > tol))
          answer_5 = 0;  // take-off = 0
        else if ((pn > tol) and (normal_residual < tol)) {
          if ((norm_pt - mu * pn > tol)) {
            // std::cout << "WARNING: the impulse is outside the Coulomb cone
            // " << std::endl;
            answer_5 = -3;  // outside the cone = -3
          } else if ((norm_pt - mu * pn < -tol)) {
            // std::cout << "the impulse is in the *interior* of  the Coulomb
            // cone  " << std::endl; std::cout << "norm_vt_plus  " <<
            // norm_vt_plus << std::endl;
            if (norm_vt_plus > tol) {
              // std::cout << "WARNING: but the norm of vt is not zero  " <<
              // std::endl;
              answer_5 = -2;  // sticking with a non zero slifing velocity = -2
            }
            // ?? answer_5 = -2 is always overwritten
            answer_5 = 1;  // sticking = 1
          } else {
            // std::cout << "the impulse is on the *boundary* of the Coulomb
            // cone  " << std::endl; std::cout << "norm_vt_plus  " <<
            // norm_vt_plus << std::endl;
            answer_5 = 2;  // sliding = 2
          }
        } else
          answer_5 = -1;  // undetermined = -1
        if ((pn > tol) and (vn_minus > tol)) {
          // std::cout << "WARNING: we apply the impact law of positive
          // velocity " << std::endl; std::cout << "pn " << pn << " vn minus "
          // << vn_minus << " vn plus " << vn_plus
          // 		<< " normal_contact_work " << normal_contact_work
          // 		<< " -e * vn_minus   " << -e*vn_minus
          // 		<< std::endl;
          answer_6 = normal_contact_work;
        }

        attr<"work_info">(*self()).push_back(
            {answer_0, answer_1, answer_2, answer_3, answer_4, answer_5, answer_6});
      };

      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& ydot_ks = storage::attr_values<interaction, "ydot">(data, step - 1);
      auto& ys = storage::attr_values<interaction, "y">(data, step);
      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);
      auto& activations =
          storage::prop_values<interaction, "activation", topology>(data, step);
      for (auto [ydot, ydot_k, y, lambda, index_nslaw, activation] :
           view::zip(ydots, ydot_ks, ys, lambdas, nslaws, activations)) {
        if (activation) add_contact_work(ydot, ydot_k, y, lambda, index_nslaw, gap_law);
      }

      if constexpr (!match::empty_item<dinteraction>) {
        auto& d_ydots = storage::attr_values<dinteraction, "ydot">(data, step);
        auto& d_ydot_ks = storage::attr_values<dinteraction, "ydot">(data, step - 1);
        auto& d_ys = storage::attr_values<dinteraction, "y">(data, step);
        auto& d_lambdas = storage::attr_values<dinteraction, "lambda">(data, step);
        auto& d_nslaws = storage::attr_values<dinteraction, "nslaw">(data, step);
        auto& d_activations =
            storage::prop_values<dinteraction, "activation", topology>(data, step);
        for (auto [ydot, ydot_k, y, lambda, index_nslaw, activation] :
             view::zip(d_ydots, d_ydot_ks, d_ys, d_lambdas, d_nslaws, d_activations)) {
          if (activation) add_contact_work(ydot, ydot_k, y, lambda, index_nslaw, false);
        }
      }

      return algebra::matrix_view<algebra::unbounded_col_matrix<scalar, 7>>(
          attr<"work_info">(*self()).data()->data(), attr<"work_info">(*self()).size(),
          attr<"work_info">(*self()).data()->size());
    }

    auto methods() {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("p0s", &interface<Handle>::p0s<indice>),
          method("radii", &interface<Handle>::radii<indice>),
          method("positions", &interface<Handle>::positions<indice>),
          method("displacements", &interface<Handle>::displacements<indice>),
          method("stiffness_matrices", &interface<Handle>::stiffness_matrices<indice>),
          method("mass_matrices", &interface<Handle>::mass_matrices<indice>),
          method("fem_velocities", &interface<Handle>::fem_velocities<indice>),
          method("kinetic_energies", &interface<Handle>::kinetic_energies<indice>),
          method("work_of_forces", &interface<Handle>::work_of_forces<indice, scalar, scalar>),
          method("velocities", &interface<Handle>::velocities<indice>),
          method("contact_points", &interface<Handle>::contact_points<indice>),
          method("contact_info", &interface<Handle>::contact_info<indice>),
          method("contact_work",
                 &interface<Handle>::contact_work<indice, scalar, scalar, int, scalar>));
    }
  };
};
}  // namespace siconos::vkernel::io
