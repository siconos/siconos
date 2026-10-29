#pragma once

#include <fstream>
#include <print>

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/index_sets.hpp"
#include "siconos/simul/moreau_jean_kernels.hpp"
#include "siconos/simul/relation_sets.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

// Topology is the item that declares the attached storages of the systems and
// interactions (involved, index, ds1, ...): their lookups name it as declarer.
template <typename System, typename Inter, typename OsiAssembled, typename Topology>
struct moreau_jean_element : item {
  using system = System;
  using interaction = Inter;
  using osi_assembled_t = OsiAssembled;
  using topology = Topology;
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

  using properties =
      gather<minv_storage_attachment, storage::keep<minv_storage_attachment, 1>,
             storage::keep<attr_t<system, "fext">, 1>, storage::keep<attr_t<system, "q">, 2>,
             storage::keep<attr_t<system, "velocity">, 2>, storage::keep<y, 2>,
             storage::keep<ydot, 2>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    static constexpr bool runtime_dof() {
      using data_t = typename Handle::data_t;
      return model::runtime_dof<decltype(storage::make_handle<system>(
          std::declval<data_t>(), storage::index<system, int>{0}))>();
    }

    decltype(auto) sum_dofs() { return storage::attr<"sum_dofs">(*self()); }

    void compute_total_dofs(auto step) {
      if constexpr (runtime_dof()) {
        // runtime degrees of freedom: sum_dofs()[k] is the first dof of the
        // involved system of rank k (the rank is the order in which the
        // systems got involved, not the order of the systems): the dofs of
        // each system at k + 1, then a scan
        using involved_systems = kernels::involved<kernels::systems>;
        auto n = kernels::count<involved_systems, kernel_context>(self()->data(), step);
        sum_dofs().assign(n + 1, 0);
        run<moreau_jean::dofs>(step, moreau_jean::dof_offsets{sum_dofs()});
        kernels::inclusive_scan(sum_dofs());
      } else {
        // compile time systems => sum dofs useless
      }
    }

    decltype(auto) total_dofs() {
      if constexpr (runtime_dof()) {
        return sum_dofs()[sum_dofs().size() - 1];
      } else {
        // compile time degree of freedom: same dof for all systems of this
        // element
        using env_t = decltype(self()->env());
        return (typename traits::config<env_t>::template convert<
                    typename interaction::dof>::type{}
                    .value) *
               number_of_involved_ds();
      }
    }

    decltype(auto) assembled_osi() {
      return storage::make_ref_handle(self()->data(), storage::attr<"assembled_osi">(*self()));
    }

    decltype(auto) ds_offset() { return storage::attr<"ds_offset">(*self()); }
    decltype(auto) inter_offset() { return storage::attr<"inter_offset">(*self()); }

    decltype(auto) nslaw_size() { return nslaw_size_t{}.value; }

    static constexpr auto nslaw_with_friction() {
      return std::derived_from<nslaw, model::newton_impact_friction>;
    }

    static constexpr auto system_with_k_matrix() { return model::has_k_matrix(system{}); }

    decltype(auto) number_of_involved_ds() { return attr<"number_of_involved_ds">(*self()); };
    decltype(auto) number_of_interactions() {
      return attr<"number_of_interactions">(*self());
    };
    decltype(auto) theta() { return assembled_osi().theta(); }
    decltype(auto) gamma() { return assembled_osi().gamma(); }
    decltype(auto) constraint_activation_threshold() {
      return assembled_osi().constraint_activation_threshold();
    }
    decltype(auto) mass_matrix_assembled() {
      auto mass_matrix_storage =
          get_storage_type(self()->data(), attr_t<system, "mass_matrix">{});

      return algebra::mat_view(mass_matrix_storage, assembled_osi().mass_matrix_assembled(),
                               ds_offset(), ds_offset());
    }

    decltype(auto) k_matrix_assembled() {
      auto k_matrix_storage = get_storage_type(self()->data(), attr_t<system, "k_matrix">{});

      return algebra::mat_view(k_matrix_storage, assembled_osi().k_matrix_assembled(),
                               ds_offset(), ds_offset());
    }

    decltype(auto) h_matrix_assembled() {
      auto h_matrix2_storage = get_storage_type(self()->data(), h_matrix2{});

      return algebra::mat_view(h_matrix2_storage, assembled_osi().h_matrix_assembled(),
                               inter_offset(), ds_offset());
    }

    decltype(auto) w_matrix_assembled() {
      auto w_matrix_storage =
          convert_storage_type(system{}, self()->data(),
                               some::matrix<some::scalar, nth_t<0, typename h_matrix1::sizes>,
                                            nth_t<0, typename h_matrix1::sizes>>{});

      return algebra::mat_view(w_matrix_storage, assembled_osi().w_matrix_assembled(),
                               inter_offset(), inter_offset());
    }

    decltype(auto) q_nsp_vector_assembled() {
      auto lambda_storage = get_storage_type(self()->data(), lambda{});

      return algebra::vec_view(lambda_storage, assembled_osi().q_nsp_vector_assembled(),
                               inter_offset());
    }
    decltype(auto) velocity_vector_assembled() {
      auto velocity_storage = get_storage_type(self()->data(), velocity{});

      return algebra::vec_view(velocity_storage, assembled_osi().velocity_vector_assembled(),
                               ds_offset());
    }
    decltype(auto) lambda_vector_assembled() {
      auto lambda_storage = get_storage_type(self()->data(), lambda{});

      return algebra::vec_view(lambda_storage, assembled_osi().lambda_vector_assembled(),
                               inter_offset());
    }
    decltype(auto) p0_vector_assembled() {
      using env_t = decltype(self()->env());
      if constexpr (!runtime_dof()) {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, typename interaction::dof>>::type;

        return algebra::vec_view<vec_t>(assembled_osi().p0_vector_assembled(), ds_offset());
      } else {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<1>>>::type;

        return algebra::vec_view<vec_t>(assembled_osi().p0_vector_assembled(), ds_offset());
      }
    }
    decltype(auto) y_vector_assembled() {
      auto y_storage = get_storage_type(self()->data(), y{});

      return algebra::vec_view(y_storage, assembled_osi().y_vector_assembled(),
                               inter_offset());
    }

    decltype(auto) ydot_vector_assembled() {
      auto ydot_storage = get_storage_type(self()->data(), ydot{});

      return algebra::vec_view(ydot_storage, assembled_osi().ydot_vector_assembled(),
                               inter_offset());
    }

    decltype(auto) mu_vector_assembled() {
      if constexpr (nslaw_with_friction()) {
        using env_t = decltype(self()->env());
        using vec_mu_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<1>>>::type;
        return algebra::vec_view<vec_mu_t>(assembled_osi().mu_vector_assembled(),
                                           inter_offset() / 2);
      } else {
        // cf
        // https://stackoverflow.com/questions/38304847/constexpr-if-and-static-assert
        []<bool flag = false>() { static_assert(flag, "no mu vector with this nslaw"); }();
      }
    }

    void initialize(auto step) {
      run<moreau_jean::zero_velocity>(step, nullptr);

      // useless at initialisation ?
      run<moreau_jean::zero_multipliers>(step, nullptr);

      compute_iteration_matrix(step);
      compute_h_matrices(step);
    }

    void compute_iteration_matrix(auto step) {
      if constexpr (system_with_k_matrix()) {
        // M^-1 (f - K q) with the factorization of M kept on the system
        using scalar = typename decltype(self()->env())::scalar;
        auto factor = [&](auto i) -> auto& { return fem_m_factor(step, i); };
        run<moreau_jean::fem_minv_f<algebra::unbounded_vector<scalar>>>(
            step, moreau_jean::fem_factors{factor});
      } else {
        run<moreau_jean::minv_f>(step, nullptr);
      }
    }

    // factorizations of the FEM system i, kept on the system for the whole
    // simulation (see algebra::fem_operator)
    algebra::fem_operator& fem_operator(auto step, std::size_t i) {
      auto& ops = storage::prop_values<system, "fem_operator", topology>(self()->data(), step);
      auto& op = ops[i];
      if (!op) op = std::make_shared<algebra::fem_operator>();
      return *op;
    }

    // M + h^2 theta^2 K, refactorized only if h or theta change
    algebra::sparse_factor& fem_w_factor(auto step, std::size_t i, auto h) {
      auto& op = fem_operator(step, i);
      double theta = self()->theta();
      if (!op.w || op.h != h || op.theta != theta) {
        auto& data = self()->data();
        auto& m = storage::attr_values<system, "mass_matrix">(data, step)[i];
        auto& k = storage::attr_values<system, "k_matrix">(data, step)[i];
        algebra::sparse_factor::matrix_t w = m;
        w += (h * h * theta * theta) * k;
        w.makeCompressed();
        op.w = std::make_shared<algebra::sparse_factor>(w);
        op.h = h;
        op.theta = theta;
      }
      return *op.w;
    }

    algebra::sparse_factor& fem_m_factor(auto step, std::size_t i) {
      auto& op = fem_operator(step, i);
      if (!op.m) {
        algebra::sparse_factor::matrix_t m =
            storage::attr_values<system, "mass_matrix">(self()->data(), step)[i];
        m.makeCompressed();
        op.m = std::make_shared<algebra::sparse_factor>(m);
      }
      return *op.m;
    }

    // the involved FEM systems as factored diagonal blocks of the assembled
    // iteration matrix, at the offsets used by the assembly
    void collect_factored_blocks(auto step, auto h,
                                 std::vector<algebra::factored_block>& blocks) {
      if constexpr (system_with_k_matrix() && runtime_dof()) {
        auto& data = self()->data();
        auto& mass_matrices = storage::attr_values<system, "mass_matrix">(data, step);
        for_each_involved(step, [&](auto i, auto index) {
          blocks.push_back({static_cast<std::ptrdiff_t>(ds_offset() + sum_dofs()[index]),
                            static_cast<std::ptrdiff_t>(algebra::ncols(mass_matrices[i])),
                            &fem_w_factor(step, i, h)});
        });
      }
    }

    void compute_h_matrices(auto step) {
      run<moreau_jean::jacobians>(step, moreau_jean::at_step{step});
    }

    static constexpr bool relation_sets() {
      return requires { requires topology::relation_sets; };
    }

    // the relation sets hold every interaction (see relation_sets.hpp): once
    // per step, before the kernels on the relations
    void sync_relation_sets(auto step) {
      if constexpr (relation_sets()) {
        simul::sync_relation_sets<interaction, topology>(self()->data(), step);
      }
    }

    // f(i, index) for every involved system i, index being its index among
    // the involved systems
    void for_each_involved(auto step, auto&& f) {
      simul::for_each_involved<system, topology>(self()->data(), step, f);
    }

    void update_iteration_matrix(auto step) {
      using data_t = const std::decay_t<decltype(self()->data())>;
      if constexpr (!storage::has_property_t<attr_t<system, "fext">, property::time_invariant,
                                             data_t>() ||
                    system_with_k_matrix()) {
        // constant fext => constant iteration matrix
        compute_iteration_matrix(step);
      }
    }

    static constexpr bool time_invariant_h_matrices() {
      using data_t = const std::decay_t<typename Handle::data_t>;
      return storage::has_property_t<interaction, property::time_invariant, data_t>();
    }

    void update_h_matrices(auto step) {
      if constexpr (!time_invariant_h_matrices()) {
        compute_h_matrices(step);
      }
    };

    // jacobians for the assembly: already up to date at this step if
    // update_h_matrices computed them (q(step) has not changed since)
    void ensure_h_matrices(auto step) {
      if constexpr (time_invariant_h_matrices()) {
        compute_h_matrices(step);
      }
    }

    // the context of the kernels of this element (see simul/kernels.hpp)
    struct kernel_context {
      using system = System;
      using interaction = Inter;
      using element = moreau_jean_element;
      using topology = Topology;
    };

    template <typename Kernel>
    void run(auto step, auto&& p) {
      kernels::run<Kernel, kernel_context>(self()->data(), step, p);
    }

    // the parameters of the kernels (see moreau_jean_kernels.hpp)
    auto parameters(auto h) {
      using scalar = typename decltype(self()->env())::scalar;
      return moreau_jean::time_step<decltype(h), scalar>{h, self()->theta()};
    }

    // the part of an assembled vector that belongs to the item of a rank
    static auto part() {
      return [](auto& assembled, auto index, const auto& like) {
        if constexpr (match::fixed_size_vector<velocity>) {
          return get_vector(assembled, index);
        } else {
          return get_vector(assembled, index, like.size());
        }
      };
    }

    // free state: v(step + 1) and q(step + 1) without the contact impulses
    void compute_free_state(auto step, auto h) {
      auto p = parameters(h);

      if constexpr (system_with_k_matrix()) {
        // W dv = h theta F(step + 1) + h (1 - theta) F(step), F = f - K q,
        // W = M + h^2 theta^2 K factorized once for the simulation
        using env_t = decltype(self()->env());
        using vector = typename env_t::template unbounded_vector<typename env_t::scalar>;
        auto factor = [&](auto i) -> auto& { return fem_w_factor(step, i, h); };
        run<kernels::fused<moreau_jean::fem_free_velocity<vector>, moreau_jean::position>>(
            step, moreau_jean::fem_time_step{p.h, p.theta, factor});
      } else {
        run<kernels::fused<moreau_jean::free_velocity, moreau_jean::position>>(step, p);
      }

      run<moreau_jean::boundary_conditions>(step, p);
    }

    void compute_output(auto step) {
      run<moreau_jean::output>(step, moreau_jean::at_step{step});
    }

    void keep_lambdas(auto step) {
      run<moreau_jean::keep_multipliers>(
          step, moreau_jean::multipliers{part(), lambda_vector_assembled()});
    }

    // add the contact impulses to v(step + 1) and update q(step + 1), for the
    // involved systems only: the others keep their free state
    void update_velocities_and_positions(auto step, auto h) {
      auto base = parameters(h);
      auto p = moreau_jean::time_step{base.h, base.theta, part(), velocity_vector_assembled()};

      run<kernels::fused<moreau_jean::contact_velocity, moreau_jean::involved_position>>(step,
                                                                                         p);

      // the involved systems may have broken them
      run<moreau_jean::boundary_conditions>(step, p);
    }

    void update_positions(auto step, auto h) {
      run<moreau_jean::position>(step, parameters(h));
    }

    auto assemble_vectors(auto step) {
      if constexpr (nslaw_with_friction()) {
        run<kernels::fused<moreau_jean::assemble_multipliers, moreau_jean::assemble_friction>>(
            step, moreau_jean::multipliers{part(), lambda_vector_assembled(),
                                           ydot_vector_assembled(), mu_vector_assembled()});
      } else {
        run<moreau_jean::assemble_multipliers>(
            step, moreau_jean::multipliers{part(), lambda_vector_assembled(),
                                           ydot_vector_assembled()});
      }
    }

    void assemble_mass_matrix_for_involved_ds(auto step) {
      // size may be 0
      run<moreau_jean::diagonal_block<"mass_matrix">>(
          step, moreau_jean::block_matrix{mass_matrix_assembled()});
    }

    void assemble_k_matrix_for_involved_ds(auto step) {
      if constexpr (system_with_k_matrix()) {
        // size may be 0
        run<moreau_jean::diagonal_block<"k_matrix">>(
            step, moreau_jean::block_matrix{k_matrix_assembled()});
      }
      // else the system is rigid and the sparse k_matrix is not filled.
    }

    // nonsmooth law effect
    void nsl_effect_on_free_output(auto step) { run<moreau_jean::impact_law>(step, nullptr); }
  };
};

template <typename System, typename OsiAssembled, typename Topology>
struct moreau_jean_element<System, empty_item, OsiAssembled, Topology> : empty_item {
  using attributes = gather<>;
};

template <typename T, typename Topology>
struct moreau_jean_element<empty_item, empty_item, T, Topology> : empty_item {
  using attributes = gather<>;
};

}  // namespace siconos::simul
