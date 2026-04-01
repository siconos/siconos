#pragma once

#include "interaction.hpp"
#include "moreau_jean_element.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/moreau_jean_assembled.hpp"
#include "siconos/simul/moreau_jean_element.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/utils/print.hpp"
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

    using attributes =
        gather<attribute<"elements", elements_t>,
               attribute<"assembled_osi", some::item_ref<assembled_osi_t>>>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      auto assembled_osi()
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
              elem.ds_offset() = ods;
              elem.inter_offset() = ointer;
              ods += elem.total_dofs(step);

              ointer += elem.number_of_interactions() * elem.nslaw_size();

              if constexpr (elem.nslaw_with_friction()) {
                num_fric += elem.number_of_interactions();
              }
            });

        indice_t raw_ds_size = ods;
        indice_t raw_inter_size = ointer;

        if constexpr (with_k_matrix()) {
          algebra::resize(k_matrix_assembled(), raw_ds_size, raw_ds_size);
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

      void compute_w_matrix(auto step, auto time_step)
      {
        if constexpr (with_k_matrix()) {
          // only simulation with fem contains runtime inter

          auto m_matrix =
              algebra::add(1., mass_matrix_assembled(),
                           time_step * time_step * theta() * theta(),
                           k_matrix_assembled());

          compute_kkt_matrix(h_matrix_assembled(), m_matrix,
                             w_matrix_assembled());
        }
        else {
          // fem systems are not present in data
          //compute_kkt_matrix(h_matrix_assembled(), mass_matrix_assembled(),
          //                   w_matrix_assembled());
          compute_w_matrix_with_diagonal_mass_matrix();
        }
      }

      void compute_w_matrix_with_diagonal_mass_matrix()
      {
        // mass matrix is assumed to be diagonal
        auto& data = self()->data();
        using info_t = storage::get_info_t<decltype(data)>;

        using env = typename info_t::env;

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
        solve_linear_system_with_transpose(m_mat, h_mat, tmp_matrix);
        prod(h_mat, tmp_matrix, w_mat);
      }

      void compute_input() { assembled_osi().compute_input(); }

      auto methods()
      {
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;
        using scalar = typename env_t::scalar;

        return collect(
            method("compute_input", &interface<Handle>::compute_input),
            method("compute_w_matrix",
                   &interface<Handle>::compute_w_matrix<indice, scalar>));
      }
    };
  };
};
}  // namespace siconos::simul
