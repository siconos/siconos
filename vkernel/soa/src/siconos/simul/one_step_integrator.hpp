#pragma once

#include "interaction.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/moreau_jean_kinds.hpp"
#include "siconos/simul/moreau_jean_schedule.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/utils/range.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

// the item types of a tuple of item_ref that are not empty items
template <typename... Refs>
constexpr auto non_empty(mp::tuple<Refs...> refs) {
  return mp::unpack(
      mp::filter(refs,
                 [](const auto& h) constexpr {
                   using t = typename std::decay_t<decltype(h)>::type;
                   return mp::bool_c<!std::derived_from<t, empty_item>>;
                 }),
      []<typename... Kinds>(Kinds...) { return mp::tuple<typename Kinds::type...>{}; });
}

// a tuple of item_ref to the item types of a tuple
template <typename... Items>
constexpr auto refs(mp::tuple<Items...>) {
  return some::tuple<some::item_ref<Items>...>{};
}

template <typename Topology>
struct one_step_integrator {
  using ct_interaction = typename Topology::fixed_dof_interaction;
  using rt_ct_interaction = typename Topology::dynamic_dof_fixed_dof_interaction;
  using rt_rt_interaction = typename Topology::dynamic_dof_dynamic_dof_interaction;

  using ct_system = typename Topology::fixed_dof_system;
  using rt_system = typename Topology::dynamic_dof_system;

  struct moreau_jean : item {
    using topology = Topology;

    // the kinds of systems and of interactions (see moreau_jean_kinds.hpp)
    using ct_systems_t = moreau_jean_systems<ct_system, Topology>;
    using rt_systems_t = moreau_jean_systems<rt_system, Topology>;

    using ct_interactions_t =
        moreau_jean_interactions<ct_interaction, ct_systems_t, ct_systems_t, Topology>;
    using rt_ct_interactions_t =
        moreau_jean_interactions<rt_ct_interaction, ct_systems_t, rt_systems_t, Topology>;
    using rt_rt_interactions_t =
        moreau_jean_interactions<rt_rt_interaction, rt_systems_t, rt_systems_t, Topology>;

    // the non-empty kinds, in the order of the assembled problem
    using system_kinds_t = decltype(non_empty(
        mp::tuple<some::item_ref<ct_systems_t>, some::item_ref<rt_systems_t>>{}));
    using interaction_kinds_t = decltype(non_empty(
        mp::tuple<some::item_ref<ct_interactions_t>, some::item_ref<rt_ct_interactions_t>,
                  some::item_ref<rt_rt_interactions_t>>{}));

    // the stages of a time step (see moreau_jean_schedule.hpp)
    using time_step_schedule = simul::moreau_jean::time_step_schedule;

    using items = decltype(mp::concat(mp::concat(mp::tuple<topology>{}, system_kinds_t{}),
                                      interaction_kinds_t{}));

    using size_1_t = some::indice_value<1>;
    using assembled_matrix =
        some::assembled_matrix<some::matrix<some::scalar, size_1_t, size_1_t>>;
    using assembled_vector = some::assembled_vector<some::vector<some::scalar, size_1_t>>;

    struct attributes {
      some::scalar theta;
      some::scalar gamma;
      some::scalar constraint_activation_threshold;
      // the assembled problem
      assembled_matrix mass_matrix_assembled;
      assembled_matrix iteration_matrix_assembled;
      assembled_matrix k_matrix_assembled;
      assembled_matrix h_matrix_assembled;
      assembled_matrix w_matrix_assembled;
      assembled_vector lambda_vector_assembled;
      assembled_vector ydot_vector_assembled;
      assembled_vector velocity_vector_assembled;
      assembled_vector p0_vector_assembled;
      assembled_vector q_nsp_vector_assembled;
      assembled_vector mu_vector_assembled;
      // the kinds
      decltype(refs(system_kinds_t{})) systems;
      decltype(refs(interaction_kinds_t{})) interactions;
    };

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      void __init__() {
        auto add = [&]<typename Kind>(Kind kind) {
          kind = storage::add<typename Kind::type>(self()->data());
        };
        mp::for_each(systems(), add);
        mp::for_each(interactions(), add);

        // the gap at the middle of the step for the active set
        gamma() = 0.5;
      }

      void __del__() {
        // not implemented
        throw std::runtime_error("Stable pointers need to be implemented");
      }

      decltype(auto) theta() { return storage::attr<"theta">(*self()); }
      decltype(auto) gamma() { return storage::attr<"gamma">(*self()); }
      decltype(auto) constraint_activation_threshold() {
        return storage::attr<"constraint_activation_threshold">(*self());
      }

      decltype(auto) mass_matrix_assembled() {
        return storage::attr<"mass_matrix_assembled">(*self());
      }
      decltype(auto) iteration_matrix_assembled() {
        return storage::attr<"iteration_matrix_assembled">(*self());
      }
      decltype(auto) k_matrix_assembled() {
        return storage::attr<"k_matrix_assembled">(*self());
      }
      decltype(auto) h_matrix_assembled() {
        return storage::attr<"h_matrix_assembled">(*self());
      }
      decltype(auto) w_matrix_assembled() {
        return storage::attr<"w_matrix_assembled">(*self());
      }
      decltype(auto) lambda_vector_assembled() {
        return storage::attr<"lambda_vector_assembled">(*self());
      }
      decltype(auto) ydot_vector_assembled() {
        return storage::attr<"ydot_vector_assembled">(*self());
      }
      decltype(auto) velocity_vector_assembled() {
        return storage::attr<"velocity_vector_assembled">(*self());
      }
      decltype(auto) p0_vector_assembled() {
        return storage::attr<"p0_vector_assembled">(*self());
      }
      decltype(auto) q_nsp_vector_assembled() {
        return storage::attr<"q_nsp_vector_assembled">(*self());
      }
      decltype(auto) mu_vector_assembled() {
        return storage::attr<"mu_vector_assembled">(*self());
      }

      // the kinds of systems, of interactions, as handles
      decltype(auto) systems() {
        return mp::transform(storage::attr<"systems">(*self()), [&](auto kind) {
          return storage::make_handle(self()->data(), kind);
        });
      }

      decltype(auto) interactions() {
        return mp::transform(storage::attr<"interactions">(*self()), [&](auto kind) {
          return storage::make_handle(self()->data(), kind);
        });
      }

      // the kind of systems Kind. Returned by value: systems() is a fresh
      // temporary tuple
      template <typename Kind>
      auto systems_of() {
        constexpr auto n = mp::unpack(system_kinds_t{}, []<typename... Kinds>(Kinds...) {
          std::size_t i = 0;
          ((std::is_same_v<Kinds, Kind> ? false : (++i, true)) && ...);
          return i;
        });
        return std::get<n>(systems());
      }

      void initialize(auto step) {
        mp::for_each(systems(), [&](auto kind) { kind.initialize(step); });
        mp::for_each(interactions(), [&](auto kind) { kind.initialize(step); });
      }

      static constexpr auto with_k_matrix() {
        return mp::any_of(system_kinds_t{}, []<typename Kind>(Kind) {
          return typename Kind::template interface<Handle>{}.with_k_matrix();
        });
      }

      static constexpr auto with_friction() {
        return mp::any_of(interaction_kinds_t{}, []<typename Kind>(Kind) {
          return typename Kind::template interface<Handle>{}.nslaw_with_friction();
        });
      }

      // the sizes of the assembled problem, the offsets of the kinds in it
      void assemble_setup(auto step) {
        using indice_t = typename decltype(self()->env())::indice;

        indice_t ods = 0;       // offset ds
        indice_t ointer = 0;    // offset inter
        indice_t num_fric = 0;  // total number of friction nslaws for mu vector size.

        mp::for_each(systems(), [&ods, &step](auto kind) {
          kind.compute_total_dofs(step);
          kind.ds_offset() = ods;
          ods += kind.total_dofs();
        });

        mp::for_each(interactions(), [&ointer, &num_fric](auto kind) {
          kind.inter_offset() = ointer;
          ointer += kind.number_of_interactions() * kind.nslaw_size();
          if constexpr (kind.nslaw_with_friction()) {
            num_fric += kind.number_of_interactions();
          }
        });

        indice_t raw_ds_size = ods;
        indice_t raw_inter_size = ointer;

        if constexpr (with_k_matrix()) {
          algebra::resize(k_matrix_assembled(), raw_ds_size, raw_ds_size);
          algebra::resize(iteration_matrix_assembled(), raw_ds_size, raw_ds_size);
        }
        algebra::resize(mass_matrix_assembled(), raw_ds_size, raw_ds_size);
        algebra::resize(h_matrix_assembled(), raw_inter_size, raw_ds_size);
        algebra::resize(w_matrix_assembled(), raw_inter_size, raw_inter_size);
        algebra::resize(lambda_vector_assembled(), raw_inter_size);
        algebra::resize(velocity_vector_assembled(), raw_ds_size);
        algebra::resize(ydot_vector_assembled(), raw_inter_size);
        algebra::resize(q_nsp_vector_assembled(), raw_inter_size);
        algebra::resize(p0_vector_assembled(), raw_ds_size);

        if constexpr (with_friction()) {
          if (num_fric > 0) {
            algebra::resize(mu_vector_assembled(), num_fric);
          }
        }
      }

      // the active set: mark the active interactions and their systems, then
      // the ranks by scans (kernels::offsets): the row of an active
      // interaction among the active interactions of its kind, the index of
      // an involved system among the involved systems of its kind, both in
      // the order of the items. The number of involved systems.
      auto compute_active_interactions(auto step, auto h) {
        mp::for_each(systems(), [&](auto kind) { kind.clear_involved(step); });

        // a system may be involved by interactions of several kinds (a rigid
        // one by a rigid-FEM interaction)
        mp::for_each(interactions(), [&](auto kind) {
          kind.mark_active(step, h, gamma(), constraint_activation_threshold());
        });

        typename decltype(self()->env())::indice involved = 0;
        mp::for_each(systems(), [&](auto kind) { involved += kind.index_involved(step); });
        return involved;
      }

      // the blocks of H of the active interactions of a kind
      void assemble_h_matrix(auto step, auto kind) {
        using kind_t = typename decltype(kind)::type;
        kind.assemble_h_matrix(step, h_matrix_assembled(),
                               systems_of<typename kind_t::systems1>(),
                               systems_of<typename kind_t::systems2>());
      }

      // the involved FEM systems as factored diagonal blocks of the assembled
      // iteration matrix
      std::vector<algebra::factored_block> factored_blocks(auto step, auto time_step) {
        std::vector<algebra::factored_block> blocks;
        mp::for_each(systems(), [&](auto kind) {
          kind.collect_factored_blocks(step, time_step, theta(), blocks);
        });
        return blocks;
      }

      // drop the FEM factorizations, to be called after changing the mass
      // or stiffness matrix of a FEM system (h and theta changes are
      // detected)
      void invalidate_fem_operators() {
        mp::for_each(systems(), [](auto kind) { kind.invalidate_fem_operators(); });
      }

      void compute_w_matrix(auto step, auto time_step) {
        if constexpr (with_k_matrix()) {
          // stiffness matrix is present
          if (algebra::nnz(k_matrix_assembled()) == 0) {
            // fem systems are not involved in contact
            compute_w_matrix_with_diagonal_mass_matrix();
          } else {
            // fem systems are involved in contact
            algebra::add(1., mass_matrix_assembled(),
                         time_step * time_step * theta() * theta(), k_matrix_assembled(),
                         iteration_matrix_assembled());

            // H (M+ h^2 \theta^2 K)^-1 H^t, with the factorizations kept
            // on the FEM systems
            auto blocks = factored_blocks(step, time_step);
            compute_kkt_matrix(h_matrix_assembled(), iteration_matrix_assembled(),
                               w_matrix_assembled(), blocks);
          }
        } else {
          // fem systems are not present in data
          compute_w_matrix_with_diagonal_mass_matrix();
        }
      }

      void compute_w_matrix_with_diagonal_mass_matrix() {
        // mass matrix is assumed to be diagonal
        using env = decltype(self()->env());

        // the blocks of the rigid systems and of the interactions between
        // them, the first kinds
        auto ct_systems = std::get<0>(systems());
        auto ct_interactions = std::get<0>(interactions());

        // xxx properties on h_matrix1, if any, are lost here
        auto tmp_matrix =
            typename traits::config<env>::template convert<some::assembled_matrix<
                some::transposed_matrix<typename ct_interactions_t::h_matrix1>>>::type{};

        auto h_mat =
            ct_interactions.h_matrix_view(h_matrix_assembled(), ct_systems.ds_offset());
        auto m_mat = ct_systems.mass_matrix_view(mass_matrix_assembled());
        auto w_mat = ct_interactions.w_matrix_view(w_matrix_assembled());

        resize(tmp_matrix, size1(h_mat), size0(h_mat));
        if constexpr (with_k_matrix()) {
          solve_linear_system_with_transpose(cast_to_diag_mat(m_mat), h_mat, tmp_matrix);
        } else {
          solve_linear_system_with_transpose(m_mat, h_mat, tmp_matrix);
        }
        prod(h_mat, tmp_matrix, w_mat);
      }

      void compute_input(auto step, auto time_step) {
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
            // block diagonal: rigid systems are divided, FEM blocks use the
            // factorizations kept on the systems
            auto blocks = factored_blocks(step, time_step);
            solve_block_diagonal_system(iteration_matrix_assembled(), p0, velo, blocks);
          } else {
            solve_linear_system(mass_matrix, p0, velo);
          }
        } else {
          solve_linear_system(mass_matrix, p0, velo);
        }
      }

      auto methods() {
        using env_t = decltype(self()->env());
        using indice = typename env_t::indice;
        using scalar = typename env_t::scalar;

        return collect(
            method("compute_input", &interface<Handle>::compute_input<indice, scalar>),
            method("invalidate_fem_operators", &interface<Handle>::invalidate_fem_operators),
            method("compute_w_matrix", &interface<Handle>::compute_w_matrix<indice, scalar>));
      }
    };
  };
};
}  // namespace siconos::simul
