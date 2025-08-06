#pragma once

#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/moreau_jean_assembled.hpp"
#include "siconos/simul/moreau_jean_element.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/utils/print.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

template <typename Topology>
struct one_step_integrator {
  using indice_t = std::size_t;

  using ct_interaction = typename Topology::fixed_dof_interaction;
  using rt_interaction = typename Topology::dynamic_dof_interaction;

  using ct_system = typename Topology::fixed_dof_system;
  using rt_system = typename Topology::dynamic_dof_system;

  struct moreau_jean : item<> {
    using topology = Topology;

    using systems_t = gather<ct_system, rt_system>;
    using interactions_t = gather<ct_interaction, rt_interaction>;

    using ct_osi_t =
        moreau_jean_element<ct_system, ct_interaction, moreau_jean_assembled>;
    using rt_osi_t =
        moreau_jean_element<rt_system, rt_interaction, moreau_jean_assembled>;

    using raw_elements =
        ground::tuple<some::item_ref<ct_osi_t>, some::item_ref<rt_osi_t>>;

    using elements = decltype(ground::unpack(
        ground::filter(
            raw_elements{},
            [](const auto& h) constexpr {
              using t = typename std::decay_t<decltype(h)>::type;
              return ground::bool_c<!std::derived_from<t, empty_item>>;
            }),
        []<typename... Elems>(Elems...) { return some::tuple<Elems...>{}; }));

    using assembled_osi_t = moreau_jean_assembled;

    using items = gather<topology, ct_osi_t, rt_osi_t, moreau_jean_assembled>;

    using attributes =
        gather<attribute<"elements", elements>,
               attribute<"assembled_osi", some::item_ref<assembled_osi_t>>>;

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      auto assembled_osi()
      {
        return storage::handle(self()->data(),
                               attr<"assembled_osi">(*self()));
      }
      void __init__()
      {
        assembled_osi() = storage::add<assembled_osi_t>(self()->data());

        ground::for_each(elements(), [&]<typename Elem>(Elem elem) {
          elem = storage::add<typename Elem::type>(self()->data());
          elem.assembled_osi() = self()->assembled_osi();
        });
      }

      void __del__()
      {
        // not implemented
        throw std::runtime_error("Stable pointers need to be implemented");
        storage::remove(assembled_osi());

        ground::for_each(elements(), [&]<typename Elem>(Elem elem) {
          storage::remove(elem);
        });
      }

      void initialize(auto step)
      {
        ground::for_each(elements(),
                         [&](auto elem) { elem.initialize(step); });
      }

      decltype(auto) theta() { return assembled_osi().theta(); }

      decltype(auto) gamma() { return assembled_osi().gamma(); }

      decltype(auto) constraint_activation_threshold()
      {
        return assembled_osi().constraint_activation_threshold();
      }

      decltype(auto) elements()
      {
        return ground::transform(
            storage::attr<"elements">(*self()),
            [&](auto elem) { return storage::handle(self()->data(), elem); });
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

      void assemble_setup()
      {
        using env_t = decltype(self()->env());
        using indice_t = typename env_t::indice;

        indice_t ods = 0;
        indice_t ointer = 0;

        ground::for_each(elements(), [&ods, &ointer](auto elem) {
          elem.ds_offset() = ods;
          elem.inter_offset() = ointer;

          ods += elem.number_of_involved_ds() * elem.dof();
          ointer += elem.number_of_interactions() * elem.nslaw_size();
        });

        assembled_osi().assemble_setup(ods, ointer);
      }

      void compute_w_matrix(auto step)
      {
        auto& data = self()->data();

        // get the number of interactions where dof is defined at
        // runtime
        auto runtime_inters =
          storage::prop_values<rt_interaction, "nds">(data, step);

        if (std::size(runtime_inters) > 0) {
          // general case: mass matrix is not diagonal
          assembled_osi().compute_w_matrix(step);
        }
        else {
          // mass matrix is assumed to be diagonal
          auto& data = self()->data();
          using info_t = storage::get_info_t<decltype(data)>;

          using env = typename info_t::env;
          auto tmp_matrix = typename traits::config<env>::template convert<
              some::assembled_matrix<some::transposed_matrix<
                  typename ct_osi_t::h_matrix1>>>::type{};

          auto ct_elem = std::get<0>(elements());

          auto&& h_mat = ct_elem.h_matrix_assembled();
          auto&& m_mat = ct_elem.mass_matrix_assembled();
          auto&& w_mat = ct_elem.w_matrix_assembled();

          resize(tmp_matrix, size1(h_mat), size0(h_mat));
          solvet(m_mat, h_mat, tmp_matrix);
          resize(w_mat, size0(h_mat), size1(tmp_matrix));

          prod(h_mat, tmp_matrix, w_mat);

        }
      }

      void compute_input() { assembled_osi().compute_input(); }

      auto methods()
      {
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;

        return collect(
            method("compute_input", &interface<Handle>::compute_input),
            method("compute_w_matrix",
                   &interface<Handle>::compute_w_matrix<indice>));
      }
    };
  };
};
}  // namespace siconos::simul
