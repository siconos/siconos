#pragma once

#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/storage/storage.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

template <typename DynamicalSystem, typename Interaction>
struct one_step_integrator {
  using items = gather<DynamicalSystem, Interaction>;
  using interaction_i = Interaction;
  using nonsmooth_law = typename interaction_i::nslaw;
  using system_i = DynamicalSystem;
  using dof = typename interaction_i::dof;
  using nslaw_size = typename interaction_i::nslaw_size;
  using nslaw = typename interaction_i::nslaw;
  using y = attr_t<interaction_i, "y">;
  using ydot = attr_t<interaction_i, "ydot">;
  using lambda = attr_t<interaction_i, "lambda">;
  using relation = attr_t<interaction_i, "relation">;
  using h_matrix1 = attr_t<interaction_i, "h_matrix1">;
  using h_matrix2 = attr_t<interaction_i, "h_matrix2">;

  using q = attr_t<system_i, "q">;
  using velocity = attr_t<system_i, "velocity">;
  using fext = attr_t<system_i, "fext">;

  struct moreau_jean : item<> {
    using system = DynamicalSystem;
    using interaction = Interaction;

    using attributes = gather<
        attribute<"theta", some::scalar>, attribute<"gamma", some::scalar>,
        attribute<"constraint_activation_threshold", some::scalar>,
        attribute<"h_matrix_assembled", some::unbounded_matrix<h_matrix1>>,
        attribute<"mass_matrix_assembled",
                  some::unbounded_matrix<attr_t<system, "mass_matrix">>>,
        attribute<"w_matrix",
                  some::unbounded_matrix<some::matrix<
                      some::scalar, nth_t<0, typename h_matrix1::sizes>,
                      nth_t<0, typename h_matrix1::sizes>>>>,
        attribute<
            "q_nsp_vector_assembled",
            some::unbounded_vector<some::vector<some::scalar, nslaw_size>>>,
        attribute<"velocity_vector_assembled",
                  some::unbounded_vector<velocity>>,
        attribute<"y_vector_assembled", some::unbounded_vector<ydot>>,
        attribute<"ydot_vector_assembled", some::unbounded_vector<ydot>>,
        attribute<"lambda_vector_assembled", some::unbounded_vector<lambda>>,
        attribute<"p0_vector_assembled",
                  some::unbounded_vector<some::vector<some::scalar, dof>>>>;

    using properties = gather<storage::keep<attr_t<system, "q">, 2>,
                              storage::keep<attr_t<system, "velocity">, 2>,
                              storage::keep<y, 2>, storage::keep<ydot, 2>>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      decltype(auto) theta() { return attr<"theta">(*self()); }
      decltype(auto) gamma() { return attr<"gamma">(*self()); }
      decltype(auto) constraint_activation_threshold()
      {
        return attr<"constraint_activation_threshold">(*self());
      }
      decltype(auto) h_matrix_assembled()
      {
        return attr<"h_matrix_assembled">(*self());
      }
      decltype(auto) q_nsp_vector_assembled()
      {
        return attr<"q_nsp_vector_assembled">(*self());
      }
      decltype(auto) velocity_vector_assembled()
      {
        return attr<"velocity_vector_assembled">(*self());
      }
      decltype(auto) lambda_vector_assembled()
      {
        return attr<"lambda_vector_assembled">(*self());
      }
      decltype(auto) p0_vector_assembled()
      {
        return attr<"p0_vector_assembled">(*self());
      }
      decltype(auto) y_vector_assembled()
      {
        return attr<"y_vector_assembled">(*self());
      }
      decltype(auto) ydot_vector_assembled()
      {
        return attr<"ydot_vector_assembled">(*self());
      }
      decltype(auto) mass_matrix_assembled()
      {
        return attr<"mass_matrix_assembled">(*self());
      }
      decltype(auto) w_matrix() { return attr<"w_matrix">(*self()); }

      auto compute_iteration_matrix(auto step)
      {
        auto &data = self()->data();
        auto &mass_matrices =
            storage::attr_memory<system, "mass_matrix">(data);
        auto &external_forces = storage::attr_memory<fext>(data);

        auto &mats = storage::memory(step, mass_matrices);
        auto &fs = storage::memory(step, external_forces);

        for (auto [mat, f] : view::zip(mats, fs)) {
          f = mat.inverse() * f;
        }
      }

      void initialize(auto step)
      {
        auto &data = self()->data();
        auto &vs_next =
            storage::attr_values<system, "velocity">(data, step + 1);
        auto &lambdas =
            storage::attr_values<interaction, "lambda">(data, step);

        auto &ydots = storage::attr_values<interaction, "ydot">(data, step);
        auto &ydots_next =
            storage::attr_values<interaction, "ydot">(data, step + 1);

        for (auto [v_next, lambda, ydot, ydot_next] :
             view::zip(vs_next, lambdas, ydots, ydots_next)) {
          algebra::set_zero(v_next);
          algebra::set_zero(lambda);
          algebra::set_zero(ydot);
          algebra::set_zero(ydot_next);
        };

        compute_iteration_matrix(step);
        compute_h_matrices(step);
      }

      void compute_h_matrices(auto step)
      {
        auto &data = self()->data();

        auto &h_matrices1 = storage::attr_values<h_matrix1>(data, step);
        auto &h_matrices2 = storage::attr_values<h_matrix2>(data, step);
        auto &ndss = storage::prop_values<interaction, "nds">(data, step);

        auto &ds1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ds2s = storage::prop_values<interaction, "ds2">(data, step);
        auto &relations = storage::attr_values<relation>(data, step);

        for (auto [rel, hm1, hm2, nds, ds1, ds2] : view::zip(
                 relations, h_matrices1, h_matrices2, ndss, ds1s, ds2s)) {
          // local binding not enough to be passed to lambda...
          auto &hhm1 = hm1;
          auto &hhm2 = hm2;
          auto hds1 = storage::handle(data, ds1);
          auto hds2 = storage::handle(data, ds2);
          auto rnds = nds;

          siconos::variant::visit(
              data, rel,
              ground::overload(
                  [&step, &hhm1, &hhm2, &hds1, &hds2,
                   &rnds](match::linear_relation auto &&rrel) {
                    if (rnds == 1) {
                      rrel.compute_jachq(step, hds1, hhm1);
                    }
                    else {
                      assert(rnds == 2);
                      rrel.compute_jachq(step, hds1, hds2, hhm1, hhm2);
                    }
                  },
                  [&step, &hhm1, &hds1](match::relation1 auto &rrel) {
                    rrel.compute_jachq(step, hds1, hhm1);
                  },
                  [&step, &hhm1, &hhm2, &hds1,
                   &hds2](match::relation2 auto &rrel) {
                    rrel.compute_jachq(step, hds1, hds2, hhm1, hhm2);
                  },
                  [](auto rrel) { assert(false); }));

          // std::cout << "HM1:" << hm1 << std::endl;
        }
      }

      void update_h_matrices(auto step)
      {
        auto &data = self()->data();
        using data_t = const std::decay_t<decltype(data)>;
        if constexpr (!storage::has_property_t<
                          interaction, property::time_invariant, data_t>()) {
          compute_h_matrices(step);
        }
      };

      void update_iteration_matrix(auto current_step)
      {
        using data_t = const std::decay_t<decltype(self()->data())>;
        if constexpr (!storage::has_property_t<attr_t<system, "fext">,
                                               property::time_invariant,
                                               data_t>()) {
          // constant fext => constant iteration matrix
          // FIX ISSUE HERE compute_iteration_matrix(current_step);
        }
      }

      // update v(step + 1)
      auto compute_free_state(auto step, auto h)
      {
        auto &data = self()->data();

        auto &vs = storage::attr_values<system, "velocity">(data, step);
        auto &vs_next =
            storage::attr_values<system, "velocity">(data, step + 1);

        auto &minv_fs = storage::attr_values<system, "fext">(data, step);
        auto &minv_fs_next =
            storage::attr_values<system, "fext">(data, step + 1);

        // for all ds
        for (auto [v, v_next, minv_f, minv_f_next] :
             view::zip(vs, vs_next, minv_fs, minv_fs_next)) {
          // note: theta useless if fext is constant
          v_next = v + h * theta() * minv_f + h * (1 - theta()) * minv_f_next;

          // std::cout << "v" << v << std::endl;
          // std::cout << "v_next:" << v_next << std::endl;
        }
      }

      void compute_output(auto step)
      {
        auto &data = self()->data();

        auto &ys = storage::attr_values<y>(data, step);
        auto &ydots = storage::attr_values<ydot>(data, step);
        auto &h_matrices1 = storage::attr_values<h_matrix1>(data, step);
        auto &h_matrices2 = storage::attr_values<h_matrix2>(data, step);

        auto &qs = storage::attr_values<system, "q">(data, step);
        auto &velocities =
            storage::attr_values<attr_t<system, "velocity">>(data, step);

        auto &ds1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ds2s = storage::prop_values<interaction, "ds2">(data, step);

        auto &ndss = storage::prop_values<interaction, "nds">(data, step);

        const auto &inters = storage::handles<interaction>(data, step);

        // global h_matrix is not assembled at this stage
        for (auto [y, ydot, hm1, hm2, ds1, ds2, nds, inter] :
             view::zip(ys, ydots, h_matrices1, h_matrices2, ds1s, ds2s, ndss,
                       inters)) {
          bool linear_case = siconos::variant::visit(
              data, inter.relation(),
              ground::overload(
                  [&](match::linear_relation auto rrel) { return true; },
                  [&](auto) { return false; }));

          if (linear_case) {
            y = hm1 * qs[ds1.get()];
            ydot = hm1 * velocities[ds1.get()];

            if (nds == 2) {
              y += hm2 * qs[ds2.get()];
              ydot += hm2 * velocities[ds2.get()];
            }
          }
          else {
            auto hds1 = storage::handle(data, ds1);
            auto hds2 = storage::handle(data, ds2);

            y[0] = siconos::variant::visit(
                data, inter.relation(),
                ground::overload(
                    [&](match::linear_relation auto rrel) {
                      assert(false);
                      return 0.;
                    },
                    [&](match::relation1 auto rrel) {
                      return rrel.compute_h(hds1);
                    },
                    [&](match::relation2 auto rrel) {
                      return rrel.compute_h(hds1, hds2);
                    }));
            ydot = hm1 * velocities[ds1.get()];
            if (nds == 2) {
              ydot += hm2 * velocities[ds2.get()];
            }
            // std::cout << "interaction:" << inter.get() << ","
            //           << "y:" << y[0] << ",ydot:" << ydot << std::endl;
          }
        }
      }

      void keep_lambdas(auto step)
      {
        auto &data = self()->data();

        auto &lambda_assembled = lambda_vector_assembled();
        auto &ydot_assembled = ydot_vector_assembled();

        auto &lambdas =
            storage::attr_values<interaction, "lambda">(data, step);
        auto &ydots_bck =
            storage::prop_values<interaction, "ydot_backup">(data, step);

        auto activations =
            storage::prop_values<interaction, "activation">(data, step);

        size_t k = 0;
        for (auto [lambda, ydot_bck, activation] :
             view::zip(lambdas, ydots_bck, activations)) {
          if (activation) {
            lambda = get_vector(lambda_assembled, k);
            ydot_bck = get_vector(ydot_assembled, k);
            k++;
          }
        }
      }

      void compute_input()
      {
        auto &h_matrix = h_matrix_assembled();
        auto &lambda = lambda_vector_assembled();
        //        auto &ydot = ydot_vector_assembled();
        auto &p0 = p0_vector_assembled();
        auto &velo = velocity_vector_assembled();
        auto &mass_matrix = mass_matrix_assembled();

        resize(p0, size1(h_matrix));
        resize(velo, size1(h_matrix));

        transpose(h_matrix);

        // p0 += h_matrix^t * lambda
        prodt1(h_matrix, lambda, p0);

        // velo <- mass_matrix^-1 * p0
        solve(mass_matrix, p0, velo);

        // ??
        // velo += h_matrix^t * ydot
        // prodt1(h_matrix, ydot, velo);

        /*      print("ydot assembled:\n");
                numerics::display(ydot);
                print("lambda assembled:\n");
                numerics::display(lambda);
                print("p0 assembled:\n");
                numerics::display(p0);
                print("velo assembled:\n");
                numerics::display(velo);*/
      }

      void update_velocities(auto step, auto h)
      {
        auto &data = self()->data();
        auto &velo = velocity_vector_assembled();

        auto &vs_next = storage::attr_values<velocity>(data, step + 1);

        auto &involveds =
            storage::prop_values<system, "involved">(data, step);

        auto &indices = storage::prop_values<system, "index">(data, step);

        // involved ds velocities -> ds velocities
        for (auto [v_next, involved, index] :
             view::zip(vs_next, involveds, indices)) {
          if (involved) {
            v_next += get_vector(velo, index);
          }
        }
      }

      void update_positions(auto step, auto h)
      {
        auto &data = self()->data();

        auto &xs = storage::attr_values<system, "q">(data, step);
        auto &xs_next = storage::attr_values<system, "q">(data, step + 1);
        auto &vs = storage::attr_values<velocity>(data, step);
        auto &vs_next = storage::attr_values<velocity>(data, step + 1);

        auto &involveds =
            storage::prop_values<system, "involved">(data, step);
        for (auto [x, x_next, v, v_next, involved] :
             view::zip(xs, xs_next, vs, vs_next, involveds)) {
          x_next = x + h * (theta() * v + (1.0 - theta()) * v_next);
        }
      }

      auto compute_active_interactions(auto step, auto h)
      {
        auto &data = self()->data();

        using info_t = storage::get_info_t<decltype(data)>;

        using env = typename info_t::env;
        using indice = typename env::indice;

        auto &ys = storage::attr_values<y>(data, step + 1);
        auto &ydots = storage::attr_values<ydot>(data, step + 1);

        auto &ids1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ids2s = storage::prop_values<interaction, "ds2">(data, step);

        auto &activations =
            storage::prop_values<interaction, "activation">(data, step);

        auto &involveds =
            storage::prop_values<system, "involved">(data, step);

        const auto &interactions = storage::handles<interaction>(data, step);

        // all ds -> not involved
        // without zip : involved is a copy not a ref!!
        for (auto [involved] : view::zip(involveds)) {
          involved = false;
        };

        auto gamma_v = 0.5;

        indice ds_counter = 0;
        indice inter_counter = 0;
        for (auto [y, ydot, activation, ids1, ids2, inter] :
             view::zip(ys, ydots, activations, ids1s, ids2s, interactions)) {
          auto b = siconos::variant::visit(
              data, inter.relation(),
              ground::overload(
                  [](match::linear_relation auto &real_relation) {
                    return real_relation.b();
                  },
                  // no b() present
                  [](auto) { return 0.; }));
          // on normal component
          // std::cout << "y:" << y[0] << " ydot:" << ydot[0]
          //           << "ACTIVATION:" << (y + gamma_v * h * ydot)(0)
          //           << std::endl;
          activation = ((y + gamma_v * h * ydot)(0) + b <=
                        self()->constraint_activation_threshold());

          if (activation) {
            inter_counter++;

            auto ds1 = storage::handle(data, ids1);
            auto ds2 = storage::handle(data, ids2);

            if (!prop<"involved">(ds1)) {
              prop<"involved">(ds1) = true;
              prop<"index">(ds1) = ds_counter++;
            };

            if (!prop<"involved">(ds2)) {
              prop<"involved">(ds2) = true;
              prop<"index">(ds2) = ds_counter++;
            }

            //          print(
            //              "\nstep: {}, time: {} => ACTIVATION {}<->{}
            //              !\ny:{}, " "ydot:{}\n",
            //             step, step * h, ids1.get(), ids2.get(), y,
            //             ydot);
          }
        }


        print(
          "  [compute_active_interactions] total number of ds: {}, total number of "
          "interactions: {}\n",
          std::size(involveds), std::size(activations));

        print(
          "  [compute_active_interactions] number of involved ds:{}, number of "
          "activated interactions: {}\n",
          ds_counter, inter_counter);
        return std::pair{inter_counter, ds_counter};
      }

      // strategy 1 : assemble the matrix for involved ds only
      auto assemble_h_matrix_for_involved_ds(auto step, auto ninter, auto nds)
      {
        auto &data = self()->data();
        auto &h_matrix = self()->h_matrix_assembled();
        resize(h_matrix, ninter, nds);

        auto &activations =
            storage::prop_values<interaction, "activation">(data, step);

        auto &h_mat1s =
            storage::attr_values<interaction, "h_matrix1">(data, step);
        auto &h_mat2s =
            storage::attr_values<interaction, "h_matrix2">(data, step);

        auto &ids1s = storage::prop_values<interaction, "ds1">(data, step);
        auto &ids2s = storage::prop_values<interaction, "ds2">(data, step);

        // auto &involveds =
        //     storage::prop_values<system, "involved">(data, step);
        auto &indices = storage::prop_values<system, "index">(data, step);

        size_t i = 0;
        for (auto [activation, h_mat1, h_mat2, ids1, ids2] :
             view::zip(activations, h_mat1s, h_mat2s, ids1s, ids2s)) {
          if (activation) {
            assert(storage::prop<"involved">(storage::handle(data, ids1)));
            assert(storage::prop<"involved">(storage::handle(data, ids2)));

            auto j1 = indices[ids1.get()];
            auto j2 = indices[ids2.get()];
            if (j1 == j2) {
              // one block
              set_value(h_matrix, i, j1,
                        h_mat1);  // sparse block matrix
            }
            else {
              // i!=j blocks
              set_value(h_matrix, i, j1,
                        h_mat1);  // sparse block matrix
              set_value(h_matrix, i, j2,
                        h_mat2);  // sparse block matrix
            }
            i++;
          }
        }
      }
      /*      print("h_matrix:\n");
              numerics::display(h_matrix_assembled());
              print("================\n");*/

      auto resize_assembled_vectors(auto step, auto ninter, auto nds)
      {
        auto &data = self()->data();

        auto &lambdas =
            storage::attr_values<interaction, "lambda">(data, step);
        auto &ydots_bck =
            storage::prop_values<interaction, "ydot_backup">(data, step);
        auto activations =
            storage::prop_values<interaction, "activation">(data, step);

        resize(self()->ydot_vector_assembled(), ninter);
        resize(self()->lambda_vector_assembled(), ninter);

        size_t k = 0;
        for (auto [lambda, ydot_bck, activation] :
             view::zip(lambdas, ydots_bck, activations)) {
          if (activation) {
            set_value(lambda_vector_assembled(), k, lambda);
            set_value(ydot_vector_assembled(), k, ydot_bck);
            k++;
          }
        }
      }

      void assemble_mass_matrix_for_involved_ds(auto step, auto size)
      {
        auto &data = self()->data();

        auto &mass_matrices =
            storage::attr_values<system, "mass_matrix">(data, step);
        auto &involveds =
            storage::prop_values<system, "involved">(data, step);
        auto &indices = storage::prop_values<system, "index">(data, step);
        // size may be 0
        resize(mass_matrix_assembled(), size, size);

        for (auto [mass_matrix, involved, index] :
             view::zip(mass_matrices, involveds, indices)) {
          if (involved) {
            set_value(mass_matrix_assembled(), index, index, mass_matrix);
          }

          /*      print("mass_matrix:\n");
                  numerics::display(mass_matrix_assembled());
                  print("================\n");
                  assert(size0(mass_matrix_assembled()) ==
                  size1(mass_matrix_assembled()));*/
        }
      }

      // compute H M^-1 H^t
      void compute_w_matrix(auto step)
      {
        auto &data = self()->data();
        using info_t = storage::get_info_t<decltype(data)>;

        using env = typename info_t::env;
        auto tmp_matrix = typename traits::config<env>::template convert<
            some::unbounded_matrix<some::transposed_matrix<h_matrix1>>>::
            type{};

        resize(tmp_matrix, size1(h_matrix_assembled()),
               size0(h_matrix_assembled()));
        // M^-1 H^t
        //        if constexpr (has_property<mass_matrix, some::diagonal>)
        //        {
        //          prod(inv(mass_matrix_assembled()),
        //          trans(h_matrix_assembled(), tmp_matrix));
        //        }
        //        else  // general case
        //        {
        solvet(mass_matrix_assembled(), h_matrix_assembled(), tmp_matrix);
        //        }

        // print("h_matrix_assembled:\n");
        // algebra::display(h_matrix_assembled());
        // print("================\n");

        // print("mass_matrix_assembled:\n");
        // algebra::display(mass_matrix_assembled());
        // print("================\n");

        // print("tmp_matrix:\n");
        // algebra::display(tmp_matrix);
        // print("================\n");

        // aliasing ?
        resize(w_matrix(), size0(h_matrix_assembled()), size1(tmp_matrix));

        prod(h_matrix_assembled(), tmp_matrix, w_matrix());

        // print("w_matrix:\n");
        // algebra::display(w_matrix());
        // print("================\n");
      }

      void nsl_effect_on_free_output(auto step)
      {
        auto &data = self()->data();
        auto &ydots = storage::attr_values<ydot>(data, step);
        auto &ydots_next = storage::attr_values<ydot>(data, step + 1);

        auto &es = storage::attr_values<nslaw, "e">(data, step);

        auto &inslaws =
            storage::attr_values<attr_t<interaction, "nslaw">>(data, step);

        for (auto [ydot, ydot_next, inslaw] :
             view::zip(ydots, ydots_next, inslaws)) {
          ydot_next += es[inslaw.get()] * ydot;
        }
      }

      // compute H vfree
      void compute_q_nsp_vector_assembled(auto step, auto ninter)
      {
        auto &data = self()->data();

        resize(q_nsp_vector_assembled(), ninter);

        auto &ydots_next = storage::attr_values<ydot>(data, step + 1);
        auto &activations =
            storage::prop_values<interaction, "activation">(data, step);

        auto k = 0;
        for (auto [ydot_next, activation] :
             view::zip(ydots_next, activations)) {
          if (activation) {
            set_value(q_nsp_vector_assembled(), k++, ydot_next);
          }
        }
      }

      auto methods()
      {
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;
        using scalar = typename env_t::scalar;

        return collect(
            method("compute_output",
                   &interface<Handle>::compute_output<indice>),
            method("compute_input", &interface<Handle>::compute_input),
            method("compute_active_interactions",
                   &interface<Handle>::compute_active_interactions<indice,
                                                                   scalar>),
            method("assemble_mass_matrix_for_involved_ds",
                   &interface<Handle>::assemble_h_matrix_for_involved_ds<
                       indice, indice, indice>),
            method(
                "resize_assembled_vectors",
                &interface<Handle>::resize_assembled_vectors<indice, indice,
                                                             indice>),
            method("assemble_mass_matrix_for_involved_ds",
                   &interface<Handle>::assemble_mass_matrix_for_involved_ds<
                       indice, indice>),
            method(
                "compute_q_nsp_vector_assembled",
                &interface<Handle>::compute_q_nsp_vector_assembled<indice,
                                                                   indice>),
            method("compute_w_matrix",
                   &interface<Handle>::compute_w_matrix<indice>),
            method("nsl_effect_on_free_output",
                   &interface<Handle>::nsl_effect_on_free_output<indice>),
            method("update_velocities",
                   &interface<Handle>::update_velocities<indice, indice>),
            method("update_positions",
                   &interface<Handle>::update_positions<indice, indice>),
            method("compute_iteration_matrix",
                   &interface<Handle>::compute_iteration_matrix<indice>),
            method("compute_free_state",
                   &interface<Handle>::compute_free_state<indice, scalar>),
            method("initialize", &interface<Handle>::initialize<indice>),
            method("compute_h_matrices",
                   &interface<Handle>::compute_h_matrices<indice>),
            method("update_iteration_matrix",
                   &interface<Handle>::update_iteration_matrix<indice>));
      }
    };
  };
};
}  // namespace siconos::simul
