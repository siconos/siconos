#pragma once

#include <fstream>
#include <print>

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

template <typename System, typename Inter, typename OsiAssembled>
struct moreau_jean_element : item {
  using system = System;
  using interaction = Inter;
  using osi_assembled_t = OsiAssembled;
  using items = gather<System, Inter, OsiAssembled>;

  using nonsmooth_law = typename interaction::nslaw;
  using nslaw_size_t = typename interaction::nslaw_size;
  using nslaw = typename interaction::nslaw;
  using y = attr_t<interaction, "y">;
  using ydot = attr_t<interaction, "ydot">;
  using lambda = attr_t<interaction, "lambda">;
  using relation = attr_t<interaction, "relation">;
  using h_matrix1 = attr_t<interaction, "h_matrix1">;
  using h_matrix2 = attr_t<interaction, "h_matrix2">;

  using q = attr_t<system, "q">;
  using velocity = attr_t<system, "velocity">;
  using fext = attr_t<system, "fext">;

  struct attributes {
    some::item_ref<osi_assembled_t> assembled_osi;
    some::unbounded_collection<some::indice> sum_dofs;
    some::indice ds_offset;
    some::indice inter_offset;
    some::indice number_of_involved_ds;
    some::indice number_of_interactions;
  };

  using minv_storage_attachment =
      storage::attached<system, symbol<"minv_f">, attr_t<system, "fext">>;

  using properties = gather<minv_storage_attachment,
                            storage::keep<minv_storage_attachment, 1>,
                            storage::keep<attr_t<system, "fext">, 1>,
                            storage::keep<attr_t<system, "q">, 2>,
                            storage::keep<attr_t<system, "velocity">, 2>,
                            storage::keep<y, 2>, storage::keep<ydot, 2>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    static constexpr bool runtime_dof()
    {
      using data_t = typename Handle::data_t;
      return model::runtime_dof<decltype(storage::make_handle<system>(
          std::declval<data_t>(), storage::index<system, int>{0}))>();
    }

    decltype(auto) sum_dofs() { return storage::attr<"sum_dofs">(*self()); }

    void compute_total_dofs(auto step)
    {
      if constexpr (runtime_dof()) {
        // runtime degrees of freedom: we must compute the sum
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;

        auto& data = self()->data();
        auto& mass_matrices =
            storage::attr_values<system, "mass_matrix">(data, step);
        auto& involveds =
            storage::prop_values<system, "involved">(data, step);

        indice sum_cols = 0;
        sum_dofs().clear();
        sum_dofs().push_back(0.);
        for (auto [mass_matrix, involved] :
             view::zip(mass_matrices, involveds)) {
          if (involved) {
            sum_cols += algebra::ncols(mass_matrix);
            sum_dofs().push_back(sum_cols);
          }
        }
      }
      else {
        // compile time systems => sum dofs useless
      }
    }

    decltype(auto) total_dofs()
    {
      if constexpr (runtime_dof()) {
        return sum_dofs()[sum_dofs().size() - 1];
      }
      else {
        // compile time degree of freedom: same dof for all systems of this
        // element
        using env_t = decltype(self()->env());
        return (typename traits::config<env_t>::template convert<
                    typename interaction::dof>::type{}
                    .value) *
               number_of_involved_ds();
      }
    }

    decltype(auto) assembled_osi()
    {
      return storage::make_ref_handle(
          self()->data(), storage::attr<"assembled_osi">(*self()));
    }

    decltype(auto) ds_offset() { return storage::attr<"ds_offset">(*self()); }
    decltype(auto) inter_offset()
    {
      return storage::attr<"inter_offset">(*self());
    }

    decltype(auto) nslaw_size() { return nslaw_size_t{}.value; }

    static constexpr auto nslaw_with_friction()
    {
      return std::derived_from<nslaw, model::newton_impact_friction>;
    }

    static constexpr auto system_with_k_matrix()
    {
      return model::has_k_matrix(system{});
    }

    decltype(auto) number_of_involved_ds()
    {
      return attr<"number_of_involved_ds">(*self());
    };
    decltype(auto) number_of_interactions()
    {
      return attr<"number_of_interactions">(*self());
    };
    decltype(auto) theta() { return assembled_osi().theta(); }
    decltype(auto) gamma() { return assembled_osi().gamma(); }
    decltype(auto) constraint_activation_threshold()
    {
      return assembled_osi().constraint_activation_threshold();
    }
    decltype(auto) mass_matrix_assembled()
    {
      auto mass_matrix_storage =
          get_storage_type(self()->data(), attr_t<system, "mass_matrix">{});

      return algebra::mat_view(mass_matrix_storage,
                               assembled_osi().mass_matrix_assembled(),
                               ds_offset(), ds_offset());
    }

    decltype(auto) k_matrix_assembled()
    {
      auto k_matrix_storage =
          get_storage_type(self()->data(), attr_t<system, "k_matrix">{});

      return algebra::mat_view(k_matrix_storage,
                               assembled_osi().k_matrix_assembled(),
                               ds_offset(), ds_offset());
    }

    decltype(auto) h_matrix_assembled()
    {
      auto h_matrix2_storage = get_storage_type(self()->data(), h_matrix2{});

      return algebra::mat_view(h_matrix2_storage,
                               assembled_osi().h_matrix_assembled(),
                               inter_offset(), ds_offset());
    }

    decltype(auto) w_matrix_assembled()
    {
      auto w_matrix_storage = convert_storage_type(
          system{}, self()->data(),
          some::matrix<some::scalar, nth_t<0, typename h_matrix1::sizes>,
                       nth_t<0, typename h_matrix1::sizes>>{});

      return algebra::mat_view(w_matrix_storage,
                               assembled_osi().w_matrix_assembled(),
                               inter_offset(), inter_offset());
    }

    decltype(auto) q_nsp_vector_assembled()
    {
      auto lambda_storage = get_storage_type(self()->data(), lambda{});

      return algebra::vec_view(lambda_storage,
                               assembled_osi().q_nsp_vector_assembled(),
                               inter_offset());
    }
    decltype(auto) velocity_vector_assembled()
    {
      auto velocity_storage = get_storage_type(self()->data(), velocity{});

      return algebra::vec_view(velocity_storage,
                               assembled_osi().velocity_vector_assembled(),
                               ds_offset());
    }
    decltype(auto) lambda_vector_assembled()
    {
      auto lambda_storage = get_storage_type(self()->data(), lambda{});

      return algebra::vec_view(lambda_storage,
                               assembled_osi().lambda_vector_assembled(),
                               inter_offset());
    }
    decltype(auto) p0_vector_assembled()
    {
      using env_t = decltype(self()->env());
      if constexpr (!runtime_dof()) {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, typename interaction::dof>>::type;

        return algebra::vec_view<vec_t>(assembled_osi().p0_vector_assembled(),
                                        ds_offset());
      }
      else {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<1>>>::type;

        return algebra::vec_view<vec_t>(assembled_osi().p0_vector_assembled(),
                                        ds_offset());
      }
    }
    decltype(auto) y_vector_assembled()
    {
      auto y_storage = get_storage_type(self()->data(), y{});

      return algebra::vec_view(
          y_storage, assembled_osi().y_vector_assembled(), inter_offset());
    }

    decltype(auto) ydot_vector_assembled()
    {
      auto ydot_storage = get_storage_type(self()->data(), ydot{});

      return algebra::vec_view(ydot_storage,
                               assembled_osi().ydot_vector_assembled(),
                               inter_offset());
    }

    decltype(auto) mu_vector_assembled()
    {
      if constexpr (nslaw_with_friction()) {
        using env_t = decltype(self()->env());
        using vec_mu_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<1>>>::type;
        return algebra::vec_view<vec_mu_t>(
            assembled_osi().mu_vector_assembled(), inter_offset());
      }
      else {
        // cf
        // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>() {
          static_assert(flag, "no mu vector with this nslaw");
        }();
      }
    }

    void initialize(auto step)
    {
      auto& data = self()->data();

      auto& vs_next =
          storage::attr_values<system, "velocity">(data, step + 1);

      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);

      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto& ydots_next =
          storage::attr_values<interaction, "ydot">(data, step + 1);

      for (auto [v_next] : view::zip(vs_next)) {
        algebra::set_zero(v_next);
      }

      // useless at initialisation ?
      for (auto [lambda, ydot, ydot_next] :
           view::zip(lambdas, ydots, ydots_next)) {
        algebra::set_zero(lambda);
        algebra::set_zero(ydot);
        algebra::set_zero(ydot_next);
      };

      compute_iteration_matrix(step);
      compute_h_matrices(step);
    }

    void compute_iteration_matrix(auto step)
    {
      auto& data = self()->data();
      using env = decltype(self()->env());
      using scalar = typename env::scalar;

      auto& mats = storage::attr_values<system, "mass_matrix">(data, step);
      auto& fs = storage::attr_values<system, "fext">(data, step);
      auto& fs_next = storage::attr_values<system, "fext">(data, step + 1);
      auto& minv_fs = storage::prop_values<system, "minv_f">(data, step);
      auto& minv_fs_next =
          storage::prop_values<system, "minv_f">(data, step + 1);

      if constexpr (system_with_k_matrix()) {
        auto& ks = storage::attr_values<system, "k_matrix">(data, step);
        auto& qs = storage::attr_values<system, "q">(data, step);
        for (auto [mat, k, q, f, f_next, minv_f, minv_f_next] :
             view::zip(mats, ks, qs, fs, fs_next, minv_fs, minv_fs_next)) {
          algebra::unbounded_vector<scalar> fmkq = f - k * q;
          algebra::solve_linear_system(mat, fmkq, minv_f);

          // constant fext
          minv_f_next = minv_f;
        }
      }
      else {
        for (auto [mat, f, minv_f, minv_f_next] :
             view::zip(mats, fs, minv_fs, minv_fs_next)) {
          minv_f = f;
          algebra::solve_in_place(mat, minv_f);
          minv_f_next = minv_f;
        }
      }
    }

    void compute_h_matrices(auto step)
    {
      auto& data = self()->data();

      auto& h_matrices1 = storage::attr_values<h_matrix1>(data, step);
      auto& h_matrices2 = storage::attr_values<h_matrix2>(data, step);
      auto& ndss = storage::prop_values<interaction, "nds">(data, step);

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);
      auto& relations = storage::attr_values<relation>(data, step);

      for (auto [rel, hm1, hm2, nds, ds1, ds2] :
           view::zip(relations, h_matrices1, h_matrices2, ndss, ds1s, ds2s)) {
        // local binding not enough to be passed to lambda...
        auto& hhm1 = hm1;
        auto& hhm2 = hm2;
        auto hds1 = storage::make_handle(data, ds1);
        auto hds2 = storage::make_handle(data, ds2);
        auto rnds = nds;

        siconos::variant::visit(
            data, rel,
            mp::overload(
                [&step, &hhm1, &hhm2, &hds1, &hds2,
                 &rnds](match::linear_relation auto&& rrel) {
                  if (rnds == 1) {
                    rrel.compute_jachq(step, hds1, hhm1);
                  }
                  else {
                    assert(rnds == 2);
                    rrel.compute_jachq(step, hds1, hds2, hhm1, hhm2);
                  }
                },
                [&step, &hhm1, &hds1](match::relation1 auto& rrel) {
                  rrel.compute_jachq(step, hds1, hhm1);
                },
                [&step, &hhm1, &hhm2, &hds1,
                 &hds2](match::relation2 auto& rrel) {
                  rrel.compute_jachq(step, hds1, hds2, hhm1, hhm2);
                },
                [](auto rrel) { assert(false); }));
      }
    }

    void update_iteration_matrix(auto step)
    {
      using data_t = const std::decay_t<decltype(self()->data())>;
      if constexpr (!storage::has_property_t<attr_t<system, "fext">,
                                             property::time_invariant,
                                             data_t>() ||
                    system_with_k_matrix()) {
        // constant fext => constant iteration matrix
        compute_iteration_matrix(step);
      }
    }

    void update_h_matrices(auto step)
    {
      auto& data = self()->data();
      using data_t = const std::decay_t<decltype(data)>;
      if constexpr (!storage::has_property_t<
                        interaction, property::time_invariant, data_t>()) {
        compute_h_matrices(step);
      }
    };

    // update v(step + 1)
    void compute_free_state(auto step, auto h)
    {
      auto& data = self()->data();
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      scalar theta = self()->theta();

      auto& vs = storage::attr_values<system, "velocity">(data, step);
      auto& vs_next =
          storage::attr_values<system, "velocity">(data, step + 1);

      auto& qs = storage::attr_values<system, "q">(data, step);
      auto& qs_next = storage::attr_values<system, "q">(data, step + 1);

      auto& bc_indices =
          storage::prop_values<system, "bc_velocities_0">(data, step);

      if constexpr (system_with_k_matrix()) {
        // W = M + h^2*theta^2*K
        // RHS = M*v - h^2*theta*(1-theta)*K*v + h*(f_ext - K*q)
        // Note: for constant fext, F_k = f_ext

        auto& fexts = storage::attr_values<system, "fext">(data, step);

        auto& mass_matrices =
            storage::attr_values<system, "mass_matrix">(data, step);
        auto& k_matrices =
            storage::attr_values<system, "k_matrix">(data, step);
        auto& vs = storage::attr_values<system, "velocity">(data, step);
        auto& vs_next =
            storage::attr_values<system, "velocity">(data, step + 1);
        auto& bc_indices =
            storage::prop_values<system, "bc_velocities_0">(data, step);

        for (auto [m_mat, k_mat, q, q_next, v, v_next, f, bc_indice] :
             view::zip(mass_matrices, k_matrices, qs, qs_next, vs, vs_next,
                       fexts, bc_indices)) {
          using vector = typename env_t::template unbounded_vector<scalar>;

          // Force at step k: F_k = f_ext - K*q_k
          vector tf = f - k_mat * q;
          vector tf_next = f - k_mat * (q + h * v);

          // W = M + h^2*theta^2*K
          auto w_mat = m_mat;
          w_mat += (h * h * theta * theta) * k_mat;

          // rhs
          vector rhs = h * theta * tf_next + h * (1 - theta) * tf;

          // W * dv = RHS
          vector dv(v_next.size());
          algebra::solve_linear_system(w_mat, rhs, dv);
          v_next = v + dv;

          // boundary conditions
          for (auto bc_idx : bc_indice) {
            v_next(bc_idx) = 0.0;
          }

          // update q next:
          q_next = q + h * theta * v_next + h * (1 - theta) * v;
        }
      }

      else {
        // K = 0
        auto& minv_fs = storage::prop_values<system, "minv_f">(data, step);
        auto& minv_fs_next =
            storage::prop_values<system, "minv_f">(data, step + 1);

        for (auto [q, q_next, v, v_next, minv_f, minv_f_next, bc_indice] :
             view::zip(qs, qs_next, vs, vs_next, minv_fs, minv_fs_next,
                       bc_indices)) {
          v_next = v + h * theta * minv_f_next + h * (1 - theta) * minv_f;

          // apply boundary conditions
          for (auto bc_idx : bc_indice) {
            v_next(bc_idx) = 0.0;
          }

          // update q next:
          q_next = q + h * theta * v_next + h * (1 - theta) * v;
        }
      }
    }

    void compute_output(auto step)
    {
      auto& data = self()->data();
      using env = decltype(self()->env());
      using scalar = typename env::scalar;
      using indice = typename env::indice;
      using nslaw_size_vector =
          typename env::template vector<scalar, nslaw_size_t{}.value>;

      auto& ys = storage::attr_values<y>(data, step);
      auto& ydots = storage::attr_values<ydot>(data, step);
      auto& h_matrices1 = storage::attr_values<h_matrix1>(data, step);
      auto& h_matrices2 = storage::attr_values<h_matrix2>(data, step);

      auto& qs = storage::attr_values<system, "q">(data, step);
      auto& velocities =
          storage::attr_values<attr_t<system, "velocity">>(data, step);

      auto& ds1s = storage::prop_values<interaction, "ds1">(data, step);
      auto& ds2s = storage::prop_values<interaction, "ds2">(data, step);

      auto& ndss = storage::prop_values<interaction, "nds">(data, step);

      const auto& inters = storage::handles<interaction>(data, step);

      // global h_matrix is not assembled at this stage
      for (auto [y, ydot, hm1, hm2, ds1, ds2, nds, inter] :
           view::zip(ys, ydots, h_matrices1, h_matrices2, ds1s, ds2s, ndss,
                     inters)) {
        bool linear_case = siconos::variant::visit(
            data, inter.relation(),
            mp::overload(
                [&](match::linear_relation auto rrel) { return true; },
                [&](auto) { return false; }));

        if (linear_case) {
          auto b = siconos::variant::visit(
              data, inter.relation(),
              mp::overload(
                  [](match::linear_relation auto& real_relation) {
                    return real_relation.b();
                  },
                  // no b() present
                  [](auto) {
                    nslaw_size_vector b;
                    b(0) = 0.;
                    return b;
                  }));

          y = hm1 * qs[ds1.value()] + b;
          ydot = hm1 * velocities[ds1.value()];

          if (nds == 2) {
            y += hm2 * qs[ds2.value()];
            ydot += hm2 * velocities[ds2.value()];
          }
        }
        else {
          auto hds1 = storage::make_handle(data, ds1);
          auto hds2 = storage::make_handle(data, ds2);

          siconos::variant::visit(
              data, inter.relation(),
              mp::overload(
                  [&](match::linear_relation auto rrel) {
                    assert(false);
                    return 0.;
                  },
                  [&](match::handle<collision::diskmesh_r> auto rrel) {
                    using vector_4d =
                        typename env::template vector<scalar, 4>;
                    using ds1_t = typename decltype(hds1)::type;

                    auto& ds1_velocities =
                        storage::attr_values<ds1_t, "velocity">(
                            self()->data(), step);

                    y[0] = rrel.compute_h(step, hds1, hds2);

                    indice i = rrel.contact_index();
                    indice idx1 = rrel.mesh().global_indices()[i];
                    indice idy1 = rrel.mesh().global_indices()[i + 1];
                    indice idx2 = rrel.mesh().global_indices()[i + 2];
                    indice idy2 = rrel.mesh().global_indices()[i + 3];

                    ydot = hm1 * ds1_velocities[ds1.value()];

                    /* segment nodes */
                    vector_4d local_velocity = {
                        velocities[ds2.value()][idx1],
                        velocities[ds2.value()][idy1],
                        velocities[ds2.value()][idx2],
                        velocities[ds2.value()][idy2]};

                    ydot += hm2 * local_velocity;
                  },
                  [&](match::relation1 auto rrel) {
                    y[0] = rrel.compute_h(step, hds1);
                    ydot = hm1 * velocities[ds1.value()];
                  },
                  [&](match::relation2 auto rrel) {
                    y[0] = rrel.compute_h(step, hds1, hds2);
                    ydot = hm1 * velocities[ds1.value()];
                    ydot += hm2 * velocities[ds2.value()];
                  }));
        }
      }
    }

    void keep_lambdas(auto step)
    {
      auto& data = self()->data();

      auto&& lambda_assembled = lambda_vector_assembled();
      auto&& ydot_assembled = ydot_vector_assembled();

      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& ydots_bck =
          storage::prop_values<interaction, "ydot_backup">(data, step);

      auto activations =
          storage::prop_values<interaction, "activation">(data, step);

      size_t k = 0;
      for (auto [lambda, ydot_bck, activation] :
           view::zip(lambdas, ydots_bck, activations)) {
        if (activation) {
          if constexpr (match::fixed_size_vector<velocity>) {
            lambda = get_vector(lambda_assembled, k);
            ydot_bck = get_vector(ydot_assembled, k);
            k++;
          }
          else {
            lambda = get_vector(lambda_assembled, k, lambda.size());
            ydot_bck = get_vector(ydot_assembled, k, ydot_bck.size());
            k++;
          }
        }
      }
    }

    void update_velocities(auto step, auto h)
    {
      auto& data = self()->data();
      auto&& velo = velocity_vector_assembled();

      auto& vs_next = storage::attr_values<velocity>(data, step + 1);

      auto& involveds = storage::prop_values<system, "involved">(data, step);

      auto& indices = storage::prop_values<system, "index">(data, step);

      auto& bc_indices =
          storage::prop_values<system, "bc_velocities_0">(data, step);

      // involved ds velocities -> ds velocities
      for (auto [v_next, involved, index, bc_indice] :
           view::zip(vs_next, involveds, indices, bc_indices)) {
        if (involved) {
          if constexpr (match::fixed_size_vector<velocity>) {
            v_next += get_vector(velo, index);
          }
          else {
            v_next += get_vector(velo, index, v_next.size());
          }

          for (auto bc_idx : bc_indice) {
            v_next(bc_idx) = 0.0;
          }
        }
      }
    }

    void update_positions(auto step, auto h)
    {
      auto& data = self()->data();

      auto& xs = storage::attr_values<system, "q">(data, step);
      auto& xs_next = storage::attr_values<system, "q">(data, step + 1);
      auto& vs = storage::attr_values<velocity>(data, step);
      auto& vs_next = storage::attr_values<velocity>(data, step + 1);

      auto& involveds = storage::prop_values<system, "involved">(data, step);
      for (auto [x, x_next, v, v_next, involved] :
           view::zip(xs, xs_next, vs, vs_next, involveds)) {
        x_next = x + h * (theta() * v_next + (1.0 - theta()) * v);
      }
    }

    auto assemble_vectors(auto step)
    {
      auto& data = self()->data();

      auto& lambdas = storage::attr_values<interaction, "lambda">(data, step);
      auto& nslaws = storage::attr_values<interaction, "nslaw">(data, step);

      auto& ydots = storage::attr_values<interaction, "ydot">(data, step);
      auto activations =
          storage::prop_values<interaction, "activation">(data, step);

      size_t k = 0;
      for (auto [lambda, nslaw, ydot, activation] :
           view::zip(lambdas, nslaws, ydots, activations)) {
        if (activation) {
          set_value(lambda_vector_assembled(), k, lambda);

          if constexpr (nslaw_with_friction()) {
            auto hnslaw = storage::make_handle(data, nslaw);
            set_value(mu_vector_assembled(), k, hnslaw.mu());
          }

          set_value(ydot_vector_assembled(), k, ydot);
          k++;
        }
      }
    }

    void assemble_mass_matrix_for_involved_ds(auto step)
    {
      auto& data = self()->data();

      auto& mass_matrices =
          storage::attr_values<system, "mass_matrix">(data, step);
      auto& involveds = storage::prop_values<system, "involved">(data, step);
      auto& indices = storage::prop_values<system, "index">(data, step);

      // size may be 0
      for (auto [mass_matrix, involved, index] :
           view::zip(mass_matrices, involveds, indices)) {
        if (involved) {
          set_value(mass_matrix_assembled(), index, index, mass_matrix);
        }
      }
    }

    void assemble_k_matrix_for_involved_ds(auto step)
    {
      if constexpr (system_with_k_matrix()) {
        auto& data = self()->data();

        auto& k_matrices =
            storage::attr_values<system, "k_matrix">(data, step);
        auto& involveds =
            storage::prop_values<system, "involved">(data, step);
        auto& indices = storage::prop_values<system, "index">(data, step);
        // size may be 0

        for (auto [k_matrix, involved, index] :
             view::zip(k_matrices, involveds, indices)) {
          if (involved) {
            set_value(k_matrix_assembled(), index, index, k_matrix);
          }
        }
      }
      // else the system is rigid and the sparse k_matrix is not filled.
    }

    // nonsmooth law effect
    void nsl_effect_on_free_output(auto step)
    {
      auto& data = self()->data();
      auto& ydots = storage::attr_values<ydot>(data, step);
      auto& ydots_next = storage::attr_values<ydot>(data, step + 1);

      auto& es = storage::attr_values<nslaw, "e">(data, step);

      auto& inslaws =
          storage::attr_values<attr_t<interaction, "nslaw">>(data, step);

      for (auto [ydot, ydot_next, inslaw] :
           view::zip(ydots, ydots_next, inslaws)) {
        ydot_next += es[inslaw.value()] * ydot;
      }
    }
  };
};

template <typename System, typename OsiAssembled>
struct moreau_jean_element<System, empty_item, OsiAssembled> : empty_item {
  using attributes = gather<>;
};

template <typename T>
struct moreau_jean_element<empty_item, empty_item, T> : empty_item {
  using attributes = gather<>;
};

}  // namespace siconos::simul
