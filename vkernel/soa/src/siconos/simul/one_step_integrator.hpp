#pragma once

#include "interaction.hpp"
#include "moreau_jean_element.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/moreau_jean_assembled.hpp"
#include "siconos/simul/moreau_jean_element.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

template <typename Topology>
struct one_step_integrator {
  using indice_t = std::size_t;

  using ct_interaction = typename Topology::fixed_dof_interaction;

  using rt_ct_interaction =
      typename Topology::dynamic_dof_fixed_dof_interaction;

  using rt_rt_interaction =
      typename Topology::dynamic_dof_dynamic_dof_interaction;

  using ct_system = typename Topology::fixed_dof_system;
  using rt_system = typename Topology::dynamic_dof_system;

  struct moreau_jean : item {
    using topology = Topology;

    using systems_t = gather<ct_system, rt_system>;
    using interactions_t =
        gather<ct_interaction, rt_ct_interaction, rt_rt_interaction>;

    using ct_osi_t =
        moreau_jean_element<ct_system, ct_interaction, moreau_jean_assembled>;
    using rt_ct_osi_t = moreau_jean_element<rt_system, rt_ct_interaction,
                                            moreau_jean_assembled>;
    using rt_rt_osi_t = moreau_jean_element<rt_system, rt_rt_interaction,
                                            moreau_jean_assembled>;

    using all_elements_t =
        mp::tuple<some::item_ref<ct_osi_t>, some::item_ref<rt_ct_osi_t>,
                  some::item_ref<rt_rt_osi_t>>;

    using elements_t = decltype(mp::unpack(
        mp::filter(all_elements_t{},
                   [](const auto& h) constexpr {
                     using t = typename std::decay_t<decltype(h)>::type;
                     return mp::bool_c<!std::derived_from<t, empty_item>>;
                   }),
        []<typename... Elems>(Elems...) { return some::tuple<Elems...>{}; }));

    using assembled_osi_t = moreau_jean_assembled;

    using items = gather<topology, ct_osi_t, rt_ct_osi_t, rt_rt_osi_t,
                         moreau_jean_assembled>;

    struct attributes {
      elements_t elements;
      some::item_ref<assembled_osi_t> assembled_osi;
    };

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      decltype(auto) assembled_osi()
      {
        return storage::make_ref_handle(self()->data(),
                                        attr<"assembled_osi">(*self()));
      }
      void __init__()
      {
        auto osi = storage::add<assembled_osi_t>(self()->data());
        assembled_osi() = osi;

        mp::for_each(elements(), [&]<typename Elem>(Elem elem) {
          elem = storage::add<typename Elem::type>(self()->data());
          elem.assembled_osi() = osi;
        });
      }

      void __del__()
      {
        // not implemented
        throw std::runtime_error("Stable pointers need to be implemented");
        storage::remove(assembled_osi());

        mp::for_each(elements(), [&]<typename Elem>(Elem elem) {
          storage::remove(elem);
        });
      }

      void initialize(auto step)
      {
        mp::for_each(elements(), [&](auto elem) { elem.initialize(step); });
      }

      decltype(auto) theta() { return assembled_osi().theta(); }

      decltype(auto) gamma() { return assembled_osi().gamma(); }

      decltype(auto) constraint_activation_threshold()
      {
        return assembled_osi().constraint_activation_threshold();
      }

      static constexpr auto index_elements() { return elements_t{}; }

      decltype(auto) elements()
      {
        return mp::transform(
            storage::attr<"elements">(*self()), [&](auto elem) {
              return storage::make_handle(self()->data(), elem);
            });
      }

      template <size_t N, typename Func>
      decltype(auto) visit_element(Func&& func)
      {
        return func(std::get<N>(elements()));
      }

      decltype(auto) mass_matrix_assembled()
      {
        return assembled_osi().mass_matrix_assembled();
      };

      decltype(auto) iteration_matrix_assembled()
      {
        return assembled_osi().iteration_matrix_assembled();
      }

      decltype(auto) k_matrix_assembled()
      {
        return assembled_osi().k_matrix_assembled();
      };

      decltype(auto) h_matrix_assembled()
      {
        return assembled_osi().h_matrix_assembled();
      };

      decltype(auto) w_matrix_assembled()
      {
        return assembled_osi().w_matrix_assembled();
      };

      decltype(auto) lambda_vector_assembled()
      {
        return assembled_osi().lambda_vector_assembled();
      };

      decltype(auto) velocity_vector_assembled()
      {
        return assembled_osi().velocity_vector_assembled();
      };

      decltype(auto) p0_vector_assembled()
      {
        return assembled_osi().p0_vector_assembled();
      };

      decltype(auto) q_nsp_vector_assembled()
      {
        return assembled_osi().q_nsp_vector_assembled();
      };

      decltype(auto) y_vector_assembled()
      {
        return assembled_osi().y_vector_assembled();
      };

      decltype(auto) ydot_vector_assembled()
      {
        return assembled_osi().ydot_vector_assembled();
      };

      decltype(auto) mu_vector_assembled()
      {
        return assembled_osi().mu_vector_assembled();
      };

      static constexpr auto with_k_matrix()
      {
        return mp::any_of(
            typename elements_t::types{}, []<typename Elem>(Elem) {
              return typename Elem::type::template interface<Handle>{}
                  .system_with_k_matrix();
            });
      }

      static constexpr auto with_friction()
      {
        return mp::any_of(
            typename elements_t::types{}, []<typename Elem>(Elem) {
              return typename Elem::type::template interface<Handle>{}
                  .nslaw_with_friction();
            });
      }

      void assemble_setup(auto step)
      {
        using env_t = decltype(self()->env());
        using indice_t = typename env_t::indice;

        indice_t ods = 0;     // offset ds
        indice_t ointer = 0;  // offset inter
        indice_t num_fric =
            0;  // total number of friction nslaws for mu vector size.

        mp::for_each(
            elements(), [&ods, &ointer, &num_fric, &step](auto elem) {
              elem.compute_total_dofs(step);

              elem.ds_offset() = ods;
              elem.inter_offset() = ointer;
              ods += elem.total_dofs();

              ointer += elem.number_of_interactions() * elem.nslaw_size();

              if constexpr (elem.nslaw_with_friction()) {
                num_fric += elem.number_of_interactions();
              }
            });

        indice_t raw_ds_size = ods;
        indice_t raw_inter_size = ointer;

        if constexpr (with_k_matrix()) {
          algebra::resize(k_matrix_assembled(), raw_ds_size, raw_ds_size);
          algebra::resize(iteration_matrix_assembled(), raw_ds_size,
                          raw_ds_size);
        }
        algebra::resize(mass_matrix_assembled(), raw_ds_size, raw_ds_size);
        algebra::resize(h_matrix_assembled(), raw_inter_size, raw_ds_size);
        algebra::resize(w_matrix_assembled(), raw_inter_size, raw_inter_size);
        algebra::resize(lambda_vector_assembled(), raw_inter_size);
        algebra::resize(velocity_vector_assembled(), raw_ds_size);
        algebra::resize(y_vector_assembled(), raw_inter_size);
        algebra::resize(ydot_vector_assembled(), raw_inter_size);
        algebra::resize(q_nsp_vector_assembled(), raw_inter_size);
        algebra::resize(p0_vector_assembled(), raw_ds_size);

        if constexpr (with_friction()) {
          if (num_fric > 0) {
            algebra::resize(mu_vector_assembled(), num_fric);
          }
        }
      }

      auto compute_active_interactions(auto step, auto h)
      {
        auto& data = self()->data();
        using env = decltype(self()->env());
        using indice = typename env::indice;

        auto& ct_ys = storage::attr_values<ct_interaction, "y">(data, step);
        auto& ct_ydots =
            storage::attr_values<ct_interaction, "ydot">(data, step);

        auto& rt_ct_ys =
            storage::attr_values<rt_ct_interaction, "y">(data, step);
        auto& rt_ct_ydots =
            storage::attr_values<rt_ct_interaction, "ydot">(data, step);

        auto& ct_ids1s =
            storage::prop_values<ct_interaction, "ds1">(data, step);
        auto& ct_ids2s =
            storage::prop_values<ct_interaction, "ds2">(data, step);
        auto& ct_ndss =
            storage::prop_values<ct_interaction, "nds">(data, step);

        auto& rt_ct_ids1s =
            storage::prop_values<rt_ct_interaction, "ds1">(data, step);
        auto& rt_ct_ids2s =
            storage::prop_values<rt_ct_interaction, "ds2">(data, step);
        auto& rt_ct_ndss =
            storage::prop_values<rt_ct_interaction, "nds">(data, step);

        auto& ct_activations =
            storage::prop_values<ct_interaction, "activation">(data, step);

        auto& rt_ct_activations =
            storage::prop_values<rt_ct_interaction, "activation">(data, step);

        auto& ct_involveds =
            storage::prop_values<ct_system, "involved">(data, step);

        auto& rt_involveds =
            storage::prop_values<rt_system, "involved">(data, step);

        const auto& ct_interactions =
            storage::handles<ct_interaction>(data, step);
        const auto& rt_ct_interactions =
            storage::handles<rt_ct_interaction>(data, step);

        // all ds -> not involved
        // view::concat: c++26
        for (auto [involved] : view::zip(ct_involveds)) {
          involved = false;
        };

        for (auto [involved] : view::zip(rt_involveds)) {
          involved = false;
        };

        auto gamma_v = 0.5;

        indice ct_ds_counter = 0;
        indice ct_inter_counter = 0;
        for (auto [y, ydot, activation, nds, ids1, ids2, inter] :
             view::zip(ct_ys, ct_ydots, ct_activations, ct_ndss, ct_ids1s,
                       ct_ids2s, ct_interactions)) {
          activation = ((y + gamma_v * h * ydot)(0) <=
                        self()->constraint_activation_threshold());

          if (activation) {
            ct_inter_counter++;

            auto ds2 = storage::make_handle(data, ids2);

            if (!prop<"involved">(ds2)) {
              prop<"involved">(ds2) = true;
              prop<"index">(ds2) = ct_ds_counter++;
            }

            if (nds == 2) {
              auto ds1 = storage::make_handle(data, ids1);

              if (!prop<"involved">(ds1)) {
                prop<"involved">(ds1) = true;
                prop<"index">(ds1) = ct_ds_counter++;
              };
            }
          }
        }

        indice rt_ds_counter = 0;
        indice rt_ct_inter_counter = 0;
        for (auto [y, ydot, activation, nds, ids1, ids2, inter] :
             view::zip(rt_ct_ys, rt_ct_ydots, rt_ct_activations, rt_ct_ndss,
                       rt_ct_ids1s, rt_ct_ids2s, rt_ct_interactions)) {
          activation = ((y + gamma_v * h * ydot)(0) <=
                        self()->constraint_activation_threshold());

          if (activation) {
            rt_ct_inter_counter++;

            auto ds2 = storage::make_handle(data, ids2);

            if (!prop<"involved">(ds2)) {
              prop<"involved">(ds2) = true;
              prop<"index">(ds2) = rt_ds_counter++;
            }

            // /!\ a rt system can only be in interaction with a ct system
            // i.e no rt rt interaction yet
            assert(nds == 2);
            {
              auto ds1 = storage::make_handle(data, ids1);

              if (!prop<"involved">(ds1)) {
                prop<"involved">(ds1) = true;
                prop<"index">(ds1) = ct_ds_counter++;
              };
            }
          }
        }

        std::print(
            "  [compute_active_interactions] total number of ct ds: "
            "{}, "
            "total "
            "number of "
            "ct interactions: {}\n",
            std::size(ct_involveds), std::size(ct_activations));

        std::print(
            "  [compute_active_interactions] total number of rt ds: "
            "{}, "
            "total "
            "number of "
            "rt ct interactions: {}\n",
            std::size(rt_involveds), std::size(rt_ct_activations));

        std::print(
            "  [compute_active_interactions] number of involved ct ds:{}, "
            "number of "
            "activated ct interactions: {}\n",
            ct_ds_counter, ct_inter_counter);

        std::print(
            "  [compute_active_interactions] number of involved rt ds:{}, "
            "number of "
            "activated ct rt interactions: {}\n",
            rt_ds_counter, rt_ct_inter_counter);

        auto ct_elem = std::get<0>(elements());
        auto rt_elem = std::get<1>(elements());

        ct_elem.number_of_interactions() = ct_inter_counter;
        ct_elem.number_of_involved_ds() = ct_ds_counter;

        rt_elem.number_of_interactions() = rt_ct_inter_counter;
        rt_elem.number_of_involved_ds() = rt_ds_counter;

        return ct_ds_counter + rt_ds_counter;
      }

      auto assemble_h_matrix_for_involved_ds(auto step)
      {
        using env_t = decltype(self()->env());
        using scalar = typename env_t::scalar;
        using matrix_1x1_t = typename env_t::template matrix<scalar, 1, 1>;

        auto& data = self()->data();

        // assumption: compile time element is the first one.
        auto ct_elem = std::get<0>(elements());
        auto rt_elem = std::get<1>(elements());

        auto&& ct_h_matrix = ct_elem.h_matrix_assembled();

        auto&& rt_h_matrix = algebra::mat_view<matrix_1x1_t>(
            assembled_osi().h_matrix_assembled(), rt_elem.inter_offset(),
            rt_elem.ds_offset());

        auto& ct_activations =
            storage::prop_values<ct_interaction, "activation">(data, step);

        auto& rt_ct_activations =
            storage::prop_values<rt_ct_interaction, "activation">(data, step);

        // auto& rt_rt_activations =
        //     storage::prop_values<rt_rt_interaction, "activation">(data,
        //     step);

        auto& ct_h_mat1s =
            storage::attr_values<ct_interaction, "h_matrix1">(data, step);
        auto& ct_h_mat2s =
            storage::attr_values<ct_interaction, "h_matrix2">(data, step);

        auto& rt_ct_h_mat1s =
            storage::attr_values<rt_ct_interaction, "h_matrix1">(data, step);
        auto& rt_ct_h_mat2s =
            storage::attr_values<rt_ct_interaction, "h_matrix2">(data, step);

        // auto& rt_rt_h_mat1s =
        //     storage::attr_values<rt_rt_interaction, "h_matrix1">(data,
        //     step);
        // auto& rt_rt_h_mat2s =
        //     storage::attr_values<rt_rt_interaction, "h_matrix2">(data,
        //     step);

        auto& ct_ids1s =
            storage::prop_values<ct_interaction, "ds1">(data, step);
        auto& ct_ids2s =
            storage::prop_values<ct_interaction, "ds2">(data, step);

        auto& rt_ct_ids1s =
            storage::prop_values<rt_ct_interaction, "ds1">(data, step);
        auto& rt_ct_ids2s =
            storage::prop_values<rt_ct_interaction, "ds2">(data, step);

        auto& rt_ct_rels =
            storage::attr_values<rt_ct_interaction, "relation">(data, step);

        // auto& rt_rt_ids1s =
        //     storage::prop_values<rt_rt_interaction, "ds1">(data, step);
        // auto& rt_rt_ids2s =
        //     storage::prop_values<rt_rt_interaction, "ds2">(data, step);

        auto& ct_indices =
            storage::prop_values<ct_system, "index">(data, step);
        auto& rt_indices =
            storage::prop_values<rt_system, "index">(data, step);

        size_t i_ct = 0;
        for (auto [activation, h_mat1, h_mat2, ids1, ids2] :
             view::zip(ct_activations, ct_h_mat1s, ct_h_mat2s, ct_ids1s,
                       ct_ids2s)) {
          // ct / ct activation (i.e disk/disk)
          if (activation) {
            auto j1 = ct_indices[ids1.value()];
            auto j2 = ct_indices[ids2.value()];

            // BC velocities for ds1
            auto handle_ds1 = storage::make_handle(data, ids1);
            auto& bc_vel_1 = storage::prop<"bc_velocities_0">(handle_ds1);

            // modification on a copy
            auto h_mat1_mod = h_mat1;

            // zero columns in h_mat1_mod / BC DOFs in ds1
            for (auto bc_local_idx : bc_vel_1) {
              h_mat1_mod.col(bc_local_idx).setZero();
            }

            if (j1 != j2) {
              // modification on a copy
              auto h_mat2_mod = h_mat2;

              auto handle_ds2 = storage::make_handle(data, ids2);
              auto& bc_vel_2 = storage::prop<"bc_velocities_0">(handle_ds2);

              // zero columns in h_mat2_mod / BC DOFs in ds2
              for (auto bc_local_idx : bc_vel_2) {
                h_mat2_mod.col(bc_local_idx).setZero();
              }

              // insertion
              set_value(ct_h_matrix, i_ct, j1, h_mat1_mod);
              set_value(ct_h_matrix, i_ct, j2, h_mat2_mod);
            }
            else {
              // self interaction
              set_value(ct_h_matrix, i_ct, j2, h_mat1_mod);
            }

            i_ct++;
          }
        }

        size_t i_rt = 0;
        for (auto [activation, h_mat1, h_mat2, ids1, ids2, rel] :
             view::zip(rt_ct_activations, rt_ct_h_mat1s, rt_ct_h_mat2s,
                       rt_ct_ids1s, rt_ct_ids2s, rt_ct_rels)) {
          if (activation) {
            // ct / rt activation (i.e disk/fem)
            auto j1 = ct_indices[ids1.value()];
            auto j2 = rt_elem.sum_dofs()[rt_indices[ids2.value()]];

            // BC velocities for ds1
            auto handle_ds1 = storage::make_handle(data, ids1);
            auto& bc_vel_1 = storage::prop<"bc_velocities_0">(handle_ds1);

            // modification on a copy
            auto h_mat1_mod = h_mat1;

            // zero columns in h_mat1_mod / BC DOFs in ds1
            for (auto bc_local_idx : bc_vel_1) {
              h_mat1_mod.col(bc_local_idx).setZero();
            }

            // insertion
            set_value(ct_h_matrix, i_ct, j1, h_mat1_mod);

            auto handle_ds2 = storage::make_handle(data, ids2);
            auto& bc_vel_2 = storage::prop<"bc_velocities_0">(handle_ds2);

            // contact index in original mesh
            auto global_index = variant::visit(
                data, rel,
                [&](match::handle<collision::diskmesh_r> auto rrel) {
                  return rrel.mesh().global_indices()[rrel.contact_index()];
                });

            // modification on a copy
            for (auto i = 0; i < algebra::nrows(h_mat2); ++i) {
              for (auto j = 0; j < algebra::ncols(h_mat2); ++j) {
                set_value(rt_h_matrix, i + i_rt * rt_elem.nslaw_size(),
                          j + j2 + global_index, h_mat2(i, j));
              }

              // bc dofs in ds2
              for (auto bc_local_idx : bc_vel_2) {
                set_value(rt_h_matrix, i + i_rt * rt_elem.nslaw_size(),
                          j2 + bc_local_idx, 0.);
              }
            }

            i_ct++;
            i_rt++;
          }
        }
      }

      void compute_w_matrix(auto step, auto time_step)
      {
        if constexpr (with_k_matrix()) {
          // stiffness matrix is present
          if (algebra::nnz(k_matrix_assembled()) == 0) {
            // fem systems are not involved in contact
            compute_w_matrix_with_diagonal_mass_matrix();
          }
          else {
            // fem systems are involved in contact
            algebra::add(1., mass_matrix_assembled(),
                         time_step * time_step * theta() * theta(),
                         k_matrix_assembled(), iteration_matrix_assembled());

            // H (M+ h^2 \theta^2 K)^-1 H^t
            compute_kkt_matrix(h_matrix_assembled(),
                               iteration_matrix_assembled(),
                               w_matrix_assembled());
          }
        }
        else {
          // fem systems are not present in data
          compute_w_matrix_with_diagonal_mass_matrix();
        }
      }

      void compute_w_matrix_with_diagonal_mass_matrix()
      {
        // mass matrix is assumed to be diagonal
        using env = decltype(self()->env());

        // xxx properties on h_matrix1, if any, are lost here
        auto tmp_matrix = typename traits::config<env>::template convert<
            some::assembled_matrix<some::transposed_matrix<
                typename ct_osi_t::h_matrix1>>>::type{};

        // assumption: compile time element is the first one.
        auto ct_elem = std::get<0>(elements());

        auto&& h_mat = ct_elem.h_matrix_assembled();
        auto&& m_mat = ct_elem.mass_matrix_assembled();
        auto&& w_mat = ct_elem.w_matrix_assembled();

        resize(tmp_matrix, size1(h_mat), size0(h_mat));
        if constexpr (with_k_matrix()) {
          solve_linear_system_with_transpose(cast_to_diag_mat(m_mat), h_mat,
                                             tmp_matrix);
        }
        else {
          solve_linear_system_with_transpose(m_mat, h_mat, tmp_matrix);
        }
        prod(h_mat, tmp_matrix, w_mat);
      }

      // compute H vfree
      void compute_q_nsp_vector_assembled(auto step)
      {
        auto& data = self()->data();
        auto ct_elem = std::get<0>(elements());
        auto rt_elem = std::get<1>(elements());

        auto& ct_ydots_next =
            storage::attr_values<ct_interaction, "ydot">(data, step + 1);
        auto& ct_activations =
            storage::prop_values<ct_interaction, "activation">(data, step);

        auto& rt_ct_ydots_next =
            storage::attr_values<rt_ct_interaction, "ydot">(data, step + 1);
        auto& rt_ct_activations =
            storage::prop_values<rt_ct_interaction, "activation">(data, step);

        auto k = 0;
        for (auto [ydot_next, activation] :
             view::zip(ct_ydots_next, ct_activations)) {
          if (activation) {
            set_value(ct_elem.q_nsp_vector_assembled(), k++, ydot_next);
          }
        }

        for (auto [ydot_next, activation] :
             view::zip(rt_ct_ydots_next, rt_ct_activations)) {
          if (activation) {
            set_value(ct_elem.q_nsp_vector_assembled(), k++, ydot_next);
          }
        }
      }

      void compute_input(auto time_step)
      {
        auto& h_matrix = h_matrix_assembled();
        auto& lambda = lambda_vector_assembled();
        auto& p0 = p0_vector_assembled();
        auto& velo = velocity_vector_assembled();
        auto& mass_matrix = mass_matrix_assembled();

        resize(p0, size1(h_matrix));
        resize(velo, size1(h_matrix));

        transpose(h_matrix);
        prodt1(h_matrix, lambda, p0);

        if constexpr (with_k_matrix()) {
          if (algebra::nnz(k_matrix_assembled()) > 0) {
            solve_linear_system(iteration_matrix_assembled(), p0, velo);
          }
          else {
            solve_linear_system(mass_matrix, p0, velo);
          }
        }
        else {
          solve_linear_system(mass_matrix, p0, velo);
        }
      }

      auto methods()
      {
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;
        using scalar = typename env_t::scalar;

        return collect(
            method("compute_input",
                   &interface<Handle>::compute_input<scalar>),
            method("compute_w_matrix",
                   &interface<Handle>::compute_w_matrix<indice, scalar>));
      }
    };
  };
};
}  // namespace siconos::simul
