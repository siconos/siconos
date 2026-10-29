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

  using rt_ct_interaction = typename Topology::dynamic_dof_fixed_dof_interaction;

  using rt_rt_interaction = typename Topology::dynamic_dof_dynamic_dof_interaction;

  using ct_system = typename Topology::fixed_dof_system;
  using rt_system = typename Topology::dynamic_dof_system;

  struct moreau_jean : item {
    using topology = Topology;

    using systems_t = gather<ct_system, rt_system>;
    using interactions_t = gather<ct_interaction, rt_ct_interaction, rt_rt_interaction>;

    using ct_osi_t =
        moreau_jean_element<ct_system, ct_interaction, moreau_jean_assembled, Topology>;
    using rt_ct_osi_t =
        moreau_jean_element<rt_system, rt_ct_interaction, moreau_jean_assembled, Topology>;
    using rt_rt_osi_t =
        moreau_jean_element<rt_system, rt_rt_interaction, moreau_jean_assembled, Topology>;

    using all_elements_t = mp::tuple<some::item_ref<ct_osi_t>, some::item_ref<rt_ct_osi_t>,
                                     some::item_ref<rt_rt_osi_t>>;

    using elements_t = decltype(mp::unpack(
        mp::filter(all_elements_t{},
                   [](const auto& h) constexpr {
                     using t = typename std::decay_t<decltype(h)>::type;
                     return mp::bool_c<!std::derived_from<t, empty_item>>;
                   }),
        []<typename... Elems>(Elems...) { return some::tuple<Elems...>{}; }));

    using assembled_osi_t = moreau_jean_assembled;

    // Only include non-empty element types in items
    using filtered_element_types_t = decltype(mp::unpack(
        mp::filter(all_elements_t{},
                   [](const auto& h) constexpr {
                     using t = typename std::decay_t<decltype(h)>::type;
                     return mp::bool_c<!std::derived_from<t, empty_item>>;
                   }),
        []<typename... Elems>(Elems...) { return mp::tuple<typename Elems::type...>{}; }));

    using items = decltype(mp::concat(mp::tuple<topology, ct_osi_t, assembled_osi_t>{},
                                      filtered_element_types_t{}));

    struct attributes {
      elements_t elements;
      some::item_ref<assembled_osi_t> assembled_osi;
    };

    template <typename Handle>
    struct interface : default_interface<Handle> {
      using default_interface<Handle>::self;

      decltype(auto) assembled_osi() {
        return storage::make_ref_handle(self()->data(), attr<"assembled_osi">(*self()));
      }
      void __init__() {
        auto osi = storage::add<assembled_osi_t>(self()->data());
        assembled_osi() = osi;

        mp::for_each(elements(), [&]<typename Elem>(Elem elem) {
          elem = storage::add<typename Elem::type>(self()->data());
          elem.assembled_osi() = osi;
        });
      }

      void __del__() {
        // not implemented
        throw std::runtime_error("Stable pointers need to be implemented");
        storage::remove(assembled_osi());

        mp::for_each(elements(), [&]<typename Elem>(Elem elem) { storage::remove(elem); });
      }

      void initialize(auto step) {
        mp::for_each(elements(), [&](auto elem) { elem.initialize(step); });
      }

      decltype(auto) theta() { return assembled_osi().theta(); }

      decltype(auto) gamma() { return assembled_osi().gamma(); }

      decltype(auto) constraint_activation_threshold() {
        return assembled_osi().constraint_activation_threshold();
      }

      static constexpr auto index_elements() { return elements_t{}; }

      decltype(auto) elements() {
        return mp::transform(storage::attr<"elements">(*self()), [&](auto elem) {
          return storage::make_handle(self()->data(), elem);
        });
      }

      template <size_t N, typename Func>
      decltype(auto) visit_element(Func&& func) {
        return func(std::get<N>(elements()));
      }

      decltype(auto) mass_matrix_assembled() {
        return assembled_osi().mass_matrix_assembled();
      };

      decltype(auto) iteration_matrix_assembled() {
        return assembled_osi().iteration_matrix_assembled();
      }

      decltype(auto) k_matrix_assembled() { return assembled_osi().k_matrix_assembled(); };

      decltype(auto) h_matrix_assembled() { return assembled_osi().h_matrix_assembled(); };

      decltype(auto) w_matrix_assembled() { return assembled_osi().w_matrix_assembled(); };

      decltype(auto) lambda_vector_assembled() {
        return assembled_osi().lambda_vector_assembled();
      };

      decltype(auto) velocity_vector_assembled() {
        return assembled_osi().velocity_vector_assembled();
      };

      decltype(auto) p0_vector_assembled() { return assembled_osi().p0_vector_assembled(); };

      decltype(auto) q_nsp_vector_assembled() {
        return assembled_osi().q_nsp_vector_assembled();
      };

      decltype(auto) y_vector_assembled() { return assembled_osi().y_vector_assembled(); };

      decltype(auto) ydot_vector_assembled() {
        return assembled_osi().ydot_vector_assembled();
      };

      decltype(auto) mu_vector_assembled() { return assembled_osi().mu_vector_assembled(); };

      static constexpr auto with_runtime_dofs() {
        return mp::any_of(typename elements_t::types{}, []<typename Elem>(Elem) {
          return typename Elem::type::template interface<Handle>{}.runtime_dof();
        });
      }

      static constexpr auto with_k_matrix() {
        return mp::any_of(typename elements_t::types{}, []<typename Elem>(Elem) {
          return typename Elem::type::template interface<Handle>{}.system_with_k_matrix();
        });
      }

      static constexpr auto with_friction() {
        return mp::any_of(typename elements_t::types{}, []<typename Elem>(Elem) {
          return typename Elem::type::template interface<Handle>{}.nslaw_with_friction();
        });
      }

      // whether this topology declares a real (non-empty) rt_ct_interaction
      // (fem/disk contact) resp. rt_rt_interaction (fem/fem "self" contact,
      // e.g. against a fixed ground) - distinct from with_k_matrix(), which
      // is true for ANY rt-family element (rt_ct or rt_rt) sharing a
      // k_matrix-bearing system, regardless of which one(s) are present.
      static constexpr auto with_rt_ct_interaction() {
        return !std::derived_from<rt_ct_interaction, empty_item>;
      }

      static constexpr auto with_rt_rt_interaction() {
        return !std::derived_from<rt_rt_interaction, empty_item>;
      }

      // rt_rt_osi_t's position in elements() is 2 when rt_ct_osi_t is also
      // present (declaration order ct, rt_ct, rt_rt), or 1 if rt_ct is
      // empty_item and thus filtered out of elements(). Returned by value
      // (not decltype(auto)): elements() is a fresh temporary tuple each
      // call, so std::get<N>(elements()) is a reference into a temporary
      // that dies at the end of this function - returning it by decltype
      // would hand the caller a dangling reference.
      auto rt_rt_element() {
        if constexpr (with_rt_ct_interaction()) {
          return std::get<2>(elements());
        } else {
          return std::get<1>(elements());
        }
      }

      void assemble_setup(auto step) {
        using env_t = decltype(self()->env());
        using indice_t = typename env_t::indice;

        indice_t ods = 0;       // offset ds
        indice_t ointer = 0;    // offset inter
        indice_t num_fric = 0;  // total number of friction nslaws for mu vector size.

        mp::for_each(elements(), [&ods, &ointer, &num_fric, &step](auto elem) {
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
          algebra::resize(iteration_matrix_assembled(), raw_ds_size, raw_ds_size);
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

      auto compute_active_interactions(auto step, auto h) {
        auto& data = self()->data();
        using env = decltype(self()->env());
        using indice = typename env::indice;

        auto ct_elem = std::get<0>(elements());

        // all ds -> not involved
        ct_elem.template run<simul::moreau_jean::not_involved>(step, nullptr);
        auto& ct_involveds = storage::prop_values<ct_system, "involved", Topology>(data, step);

        // the index sets are rebuilt with the properties
        if constexpr (with_index_sets<Topology>) {
          index_set<ct_interaction, active_set_tag, Topology>(data).clear();
          index_set<ct_system, involved_set_tag, Topology>(data).clear();
          if constexpr (with_rt_ct_interaction()) {
            index_set<rt_ct_interaction, active_set_tag, Topology>(data).clear();
          }
          if constexpr (with_rt_rt_interaction()) {
            index_set<rt_rt_interaction, active_set_tag, Topology>(data).clear();
          }
          if constexpr (!std::derived_from<rt_system, empty_item>) {
            index_set<rt_system, involved_set_tag, Topology>(data).clear();
          }
        }

        indice ct_ds_counter = 0;
        indice ct_inter_counter = 0;

        indice rt_ds_counter = 0;
        indice rt_ct_inter_counter = 0;
        indice rt_rt_inter_counter = 0;

        // the system ds is involved: the next index of its kind if it is not
        // already (a rigid system touched by a rigid-FEM interaction gets a
        // rigid index)
        auto involve = [&data, &ct_ds_counter, &rt_ds_counter](auto ds) {
          using ds_t = typename decltype(ds)::type;
          indice& counter = std::is_same_v<ds_t, ct_system> ? ct_ds_counter : rt_ds_counter;
          if (!prop<"involved", Topology>(ds)) {
            prop<"involved", Topology>(ds) = true;
            prop<"index", Topology>(ds) = counter;
            if constexpr (with_index_sets<Topology>) {
              index_set<ds_t, involved_set_tag, Topology>(data).emplace(ds.index().value(),
                                                                        counter);
            }
            counter++;
          }
        };

        // the active set of the interactions of an element, rows counted by
        // counter: the activation of each interaction, then in order the rows
        // of the active ones and the indices of their systems
        auto build_active_set = [&](auto elem, indice& counter) {
          using inter_t = typename decltype(elem)::type::interaction;
          auto activate = [&data, &counter](auto i) {
            const indice row = counter++;
            if constexpr (with_index_sets<Topology>) {
              index_set<inter_t, active_set_tag, Topology>(data).emplace(i, row);
            }
          };
          elem.template run<
              kernels::fused<simul::moreau_jean::activation, simul::moreau_jean::activate>>(
              step, simul::moreau_jean::active_set{
                        h, 0.5, self()->constraint_activation_threshold(), activate, involve});
        };

        build_active_set(ct_elem, ct_inter_counter);

        std::print(
            "  [compute_active_interactions] total number of ct ds: "
            "{}, "
            "total "
            "number of "
            "ct interactions: {}\n",
            std::size(ct_involveds),
            std::size(
                storage::prop_values<ct_interaction, "activation", Topology>(data, step)));

        if constexpr (with_k_matrix()) {
          // reset rt_system involvement once per step, regardless of which
          // rt-family interaction type(s) (rt_ct and/or rt_rt) are actually
          // present - both key off this same shared per-rt_system property.
          std::get<1>(elements())
              .template run<simul::moreau_jean::not_involved>(step, nullptr);
        }

        if constexpr (with_rt_ct_interaction()) {
          auto rt_elem = std::get<1>(elements());
          // /!\ a rt system can only be in interaction with a ct system (ds1)
          build_active_set(rt_elem, rt_ct_inter_counter);

          std::print(
              "  [compute_active_interactions] total number of rt ds: "
              "{}, "
              "total "
              "number of "
              "rt ct interactions: {}\n",
              std::size(storage::prop_values<rt_system, "involved", Topology>(data, step)),
              std::size(storage::prop_values<rt_ct_interaction, "activation", Topology>(
                  data, step)));

          rt_elem.number_of_interactions() = rt_ct_inter_counter;
          rt_elem.number_of_involved_ds() = rt_ds_counter;
        }

        if constexpr (with_rt_rt_interaction()) {
          // fem "self" contact (e.g. against a fixed ground): ds1 == ds2,
          // both indexed into the same rt_system as rt_ct above - continue
          // (not reset) rt_ds_counter so a rt_system ds already marked
          // involved via rt_ct keeps its index.
          auto rt_rt_elem = rt_rt_element();
          build_active_set(rt_rt_elem, rt_rt_inter_counter);

          std::print(
              "  [compute_active_interactions] number of activated rt rt "
              "interactions: {}\n",
              rt_rt_inter_counter);

          rt_rt_elem.number_of_interactions() = rt_rt_inter_counter;
          rt_rt_elem.number_of_involved_ds() = rt_ds_counter;
        }

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

        ct_elem.number_of_interactions() = ct_inter_counter;
        ct_elem.number_of_involved_ds() = ct_ds_counter;

        return ct_ds_counter + rt_ds_counter;
      }

      auto assemble_h_matrix_for_involved_ds(auto step) {
        using env_t = decltype(self()->env());
        using scalar = typename env_t::scalar;
        using matrix_1x1_t = typename env_t::template matrix<scalar, 1, 1>;

        auto& data = self()->data();

        // assumption: compile time element is the first one.
        auto ct_elem = std::get<0>(elements());
        auto&& ct_h_matrix = ct_elem.h_matrix_assembled();

        // ct / ct activation (i.e disk/disk)
        ct_elem.template run<simul::moreau_jean::h_blocks>(
            step, simul::moreau_jean::block_matrix{ct_h_matrix});

        // the rt ct rows follow the ct rows
        [[maybe_unused]] const std::size_t n_ct = ct_elem.number_of_interactions();

        if constexpr (with_rt_ct_interaction()) {
          // ct / rt activation (i.e disk/fem)
          auto rt_elem = std::get<1>(elements());
          auto&& rt_h_matrix =
              algebra::mat_view<matrix_1x1_t>(assembled_osi().h_matrix_assembled(),
                                              rt_elem.inter_offset(), rt_elem.ds_offset());
          rt_elem.template run<simul::moreau_jean::fem_h_blocks>(
              step, simul::moreau_jean::fem_h_matrix{
                        data, rt_h_matrix, ct_h_matrix, rt_elem.sum_dofs(), n_ct,
                        static_cast<std::size_t>(rt_elem.nslaw_size())});
        }

        if constexpr (with_rt_rt_interaction()) {
          auto rt_rt_elem = rt_rt_element();
          auto&& rt_rt_h_matrix = algebra::mat_view<matrix_1x1_t>(
              assembled_osi().h_matrix_assembled(), rt_rt_elem.inter_offset(),
              rt_rt_elem.ds_offset());
          rt_rt_elem.template run<simul::moreau_jean::fem_self_h_block>(
              step,
              simul::moreau_jean::fem_block_matrix{rt_rt_h_matrix, rt_rt_elem.sum_dofs()});
        }
      }

      // the involved FEM systems as factored diagonal blocks of the
      // assembled iteration matrix
      std::vector<algebra::factored_block> factored_blocks(auto step, auto time_step) {
        std::vector<algebra::factored_block> blocks;
        mp::for_each(elements(), [&](auto elem) {
          elem.collect_factored_blocks(step, time_step, blocks);
        });
        return blocks;
      }

      // drop the FEM factorizations, to be called after changing the mass
      // or stiffness matrix of a FEM system (h and theta changes are
      // detected)
      void invalidate_fem_operators() {
        if constexpr (with_k_matrix()) {
          auto& ops =
              storage::prop_values<rt_system, "fem_operator", Topology>(self()->data(), 0);
          for (auto& op : ops) op.reset();
        }
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

        // xxx properties on h_matrix1, if any, are lost here
        auto tmp_matrix =
            typename traits::config<env>::template convert<some::assembled_matrix<
                some::transposed_matrix<typename ct_osi_t::h_matrix1>>>::type{};

        // assumption: compile time element is the first one.
        auto ct_elem = std::get<0>(elements());

        auto&& h_mat = ct_elem.h_matrix_assembled();
        auto&& m_mat = ct_elem.mass_matrix_assembled();
        auto&& w_mat = ct_elem.w_matrix_assembled();

        resize(tmp_matrix, size1(h_mat), size0(h_mat));
        if constexpr (with_k_matrix()) {
          solve_linear_system_with_transpose(cast_to_diag_mat(m_mat), h_mat, tmp_matrix);
        } else {
          solve_linear_system_with_transpose(m_mat, h_mat, tmp_matrix);
        }
        prod(h_mat, tmp_matrix, w_mat);
      }

      // compute H vfree: the free output derivatives of the active
      // interactions, the rows of each interaction type following the rows of
      // the previous ones
      void compute_q_nsp_vector_assembled(auto step) {
        auto ct_elem = std::get<0>(elements());

        // /!\ with ct_elem offsets (0), should be directly in global
        // q_nsp_vector_assembled()
        auto q = ct_elem.q_nsp_vector_assembled();
        std::size_t k0 = 0;

        ct_elem.template run<simul::moreau_jean::free_output>(
            step, simul::moreau_jean::assembled_at{q, k0});
        k0 += ct_elem.number_of_interactions();

        if constexpr (with_rt_ct_interaction()) {
          auto rt_elem = std::get<1>(elements());
          rt_elem.template run<simul::moreau_jean::free_output>(
              step, simul::moreau_jean::assembled_at{q, k0});
          k0 += rt_elem.number_of_interactions();
        }

        if constexpr (with_rt_rt_interaction()) {
          rt_rt_element().template run<simul::moreau_jean::free_output>(
              step, simul::moreau_jean::assembled_at{q, k0});
        }
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
