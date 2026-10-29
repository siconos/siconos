#pragma once

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

/* The kinds of the Moreau-Jean integrator: a kind of systems (rigid with a
   fixed number of dofs, FEM with a dynamic one) or a kind of interactions
   (between rigid systems, rigid-FEM, FEM-FEM). A kind runs the kernels on its
   items (see moreau_jean_kernels.hpp) and keeps where they are in the
   assembled problem; the assembled matrices and vectors are those of the
   integrator, given to the kind.

   Topology is the item that declares the attached storages of the systems
   and interactions (involved, index, ds1, ...): their lookups name it as
   declarer. */

// a kind of systems
template <typename System, typename Topology>
struct moreau_jean_systems : item {
  using system = System;
  using topology = Topology;
  using items = gather<System>;

  using q = attr_t<system, "q">;
  using velocity = attr_t<system, "velocity">;
  using fext = attr_t<system, "fext">;

  struct attributes {
    some::unbounded_collection<some::indice> sum_dofs;
    some::indice ds_offset;
    some::indice number_of_involved_ds;
  };

  using minv_storage_attachment = storage::attached<system, symbol<"minv_f">, fext>;

  using properties =
      gather<minv_storage_attachment, storage::keep<minv_storage_attachment, 1>,
             storage::keep<fext, 1>, storage::keep<q, 2>, storage::keep<velocity, 2>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    static constexpr bool runtime_dof() {
      using data_t = typename Handle::data_t;
      return model::runtime_dof<decltype(storage::make_handle<system>(
          std::declval<data_t>(), storage::index<system, int>{0}))>();
    }

    static constexpr auto with_k_matrix() { return model::has_k_matrix(system{}); }

    decltype(auto) sum_dofs() { return storage::attr<"sum_dofs">(*self()); }
    decltype(auto) ds_offset() { return storage::attr<"ds_offset">(*self()); }
    decltype(auto) number_of_involved_ds() {
      return storage::attr<"number_of_involved_ds">(*self());
    }

    // the context of the kernels of this kind (see simul/kernels.hpp)
    struct kernel_context {
      using system = System;
      using kind = moreau_jean_systems;
      using topology = Topology;
    };

    template <typename Kernel>
    void run(auto step, auto&& p) {
      kernels::run<Kernel, kernel_context>(self()->data(), step, p);
    }

    // the parameters of the kernels (see moreau_jean_kernels.hpp)
    auto parameters(auto h, auto theta) {
      using scalar = typename decltype(self()->env())::scalar;
      return moreau_jean::time_step<decltype(h), scalar>{h, theta};
    }

    // the part of an assembled vector that belongs to the system of a rank
    static auto part() {
      return [](auto& assembled, auto index, const auto& like) {
        if constexpr (match::fixed_size_vector<velocity>) {
          return get_vector(assembled, index);
        } else {
          return get_vector(assembled, index, like.size());
        }
      };
    }

    // the blocks of the systems in the assembled matrices and vectors
    decltype(auto) mass_matrix_view(auto& assembled) {
      auto storage = get_storage_type(self()->data(), attr_t<system, "mass_matrix">{});
      return algebra::mat_view(storage, assembled, ds_offset(), ds_offset());
    }

    decltype(auto) k_matrix_view(auto& assembled) {
      auto storage = get_storage_type(self()->data(), attr_t<system, "k_matrix">{});
      return algebra::mat_view(storage, assembled, ds_offset(), ds_offset());
    }

    decltype(auto) velocity_view(auto& assembled) {
      auto storage = get_storage_type(self()->data(), velocity{});
      return algebra::vec_view(storage, assembled, ds_offset());
    }

    decltype(auto) p0_view(auto& assembled) {
      using env_t = decltype(self()->env());
      if constexpr (!runtime_dof()) {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, typename topology::dof>>::type;
        return algebra::vec_view<vec_t>(assembled, ds_offset());
      } else {
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<1>>>::type;
        return algebra::vec_view<vec_t>(assembled, ds_offset());
      }
    }

    void initialize(auto step) {
      run<moreau_jean::zero_velocity>(step, nullptr);
      compute_iteration_matrix(step);
    }

    // the dofs of the involved systems
    void compute_total_dofs(auto step) {
      if constexpr (runtime_dof()) {
        // runtime degrees of freedom: sum_dofs()[k] is the first dof of the
        // involved system of rank k: the dofs of each system at k + 1, then
        // a scan
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
        // compile time degree of freedom: the same for all the systems
        using env_t = decltype(self()->env());
        return (typename traits::config<env_t>::template convert<
                    typename topology::dof>::type{}
                    .value) *
               number_of_involved_ds();
      }
    }

    // the active set: no system is involved before the marks of the
    // interactions
    void clear_involved(auto step) {
      run<moreau_jean::not_involved>(step, nullptr);
      if constexpr (with_index_sets<Topology>) {
        index_set<system, involved_set_tag, Topology>(self()->data()).clear();
      }
    }

    // after the marks: the index of the involved systems, their rank among
    // the involved systems of this kind in the order of the systems; the
    // number of involved systems
    auto index_involved(auto step) {
      auto& data = self()->data();
      using indice = typename decltype(self()->env())::indice;
      auto offsets = kernels::offsets<indice>(
          storage::prop_values<system, "involved", Topology>(data, step));
      run<moreau_jean::involved_index>(step, moreau_jean::ranks{offsets});
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto index) {
          index_set<system, involved_set_tag, Topology>(data).emplace(i, index);
        };
        run<moreau_jean::involved_entry>(step, moreau_jean::ranks{offsets, emplace});
      }
      number_of_involved_ds() = offsets.back();
      return offsets.back();
    }

    // f(i, index) for every involved system i, index being its index among
    // the involved systems
    void for_each_involved(auto step, auto&& f) {
      simul::for_each_involved<system, topology>(self()->data(), step, f);
    }

    void compute_iteration_matrix(auto step) {
      if constexpr (with_k_matrix()) {
        // M^-1 (f - K q) with the factorization of M kept on the system
        using scalar = typename decltype(self()->env())::scalar;
        auto factor = [&](auto i) -> auto& { return fem_m_factor(step, i); };
        run<moreau_jean::fem_minv_f<algebra::unbounded_vector<scalar>>>(
            step, moreau_jean::fem_factors{factor});
      } else {
        run<moreau_jean::minv_f>(step, nullptr);
      }
    }

    void update_iteration_matrix(auto step) {
      using data_t = const std::decay_t<decltype(self()->data())>;
      if constexpr (!storage::has_property_t<fext, property::time_invariant, data_t>() ||
                    with_k_matrix()) {
        // constant fext => constant iteration matrix
        compute_iteration_matrix(step);
      }
    }

    // free state: v(step + 1) and q(step + 1) without the contact impulses
    void compute_free_state(auto step, auto h, auto theta) {
      auto p = parameters(h, theta);

      if constexpr (with_k_matrix()) {
        // W dv = h theta F(step + 1) + h (1 - theta) F(step), F = f - K q,
        // W = M + h^2 theta^2 K factorized once for the simulation
        using env_t = decltype(self()->env());
        using vector = typename env_t::template unbounded_vector<typename env_t::scalar>;
        auto factor = [&](auto i) -> auto& { return fem_w_factor(step, i, h, theta); };
        run<kernels::fused<moreau_jean::fem_free_velocity<vector>, moreau_jean::position>>(
            step, moreau_jean::fem_time_step{p.h, p.theta, factor});
      } else {
        run<kernels::fused<moreau_jean::free_velocity, moreau_jean::position>>(step, p);
      }

      run<moreau_jean::boundary_conditions>(step, p);
    }

    // add the contact impulses to v(step + 1) and update q(step + 1), for the
    // involved systems only: the others keep their free state
    void update_velocities_and_positions(auto step, auto h, auto theta, auto& velocities) {
      auto base = parameters(h, theta);
      auto p = moreau_jean::time_step{base.h, base.theta, part(), velocity_view(velocities)};

      run<kernels::fused<moreau_jean::contact_velocity, moreau_jean::involved_position>>(step,
                                                                                         p);

      // the involved systems may have broken them
      run<moreau_jean::boundary_conditions>(step, p);
    }

    void assemble_mass_matrix(auto step, auto& assembled) {
      // size may be 0
      run<moreau_jean::diagonal_block<"mass_matrix">>(
          step, moreau_jean::block_matrix{mass_matrix_view(assembled)});
    }

    void assemble_k_matrix(auto step, auto& assembled) {
      if constexpr (with_k_matrix()) {
        // size may be 0
        run<moreau_jean::diagonal_block<"k_matrix">>(
            step, moreau_jean::block_matrix{k_matrix_view(assembled)});
      }
      // else the systems are rigid and the sparse k_matrix is not filled.
    }

    // the factorizations of the FEM system i, kept on the system for the
    // whole simulation (see algebra::fem_operator)
    algebra::fem_operator& fem_operator(auto step, std::size_t i) {
      auto& ops = storage::prop_values<system, "fem_operator", topology>(self()->data(), step);
      auto& op = ops[i];
      if (!op) op = std::make_shared<algebra::fem_operator>();
      return *op;
    }

    // M + h^2 theta^2 K, refactorized only if h or theta change
    algebra::sparse_factor& fem_w_factor(auto step, std::size_t i, auto h, double theta) {
      auto& op = fem_operator(step, i);
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

    // drop the FEM factorizations, to be called after changing the mass or
    // stiffness matrix of a FEM system (h and theta changes are detected)
    void invalidate_fem_operators() {
      if constexpr (with_k_matrix()) {
        auto& ops = storage::prop_values<system, "fem_operator", topology>(self()->data(), 0);
        for (auto& op : ops) op.reset();
      }
    }

    // the involved FEM systems as factored diagonal blocks of the assembled
    // iteration matrix, at the offsets used by the assembly
    void collect_factored_blocks(auto step, auto h, auto theta,
                                 std::vector<algebra::factored_block>& blocks) {
      if constexpr (with_k_matrix() && runtime_dof()) {
        auto& data = self()->data();
        auto& mass_matrices = storage::attr_values<system, "mass_matrix">(data, step);
        for_each_involved(step, [&](auto i, auto index) {
          blocks.push_back({static_cast<std::ptrdiff_t>(ds_offset() + sum_dofs()[index]),
                            static_cast<std::ptrdiff_t>(algebra::ncols(mass_matrices[i])),
                            &fem_w_factor(step, i, h, theta)});
        });
      }
    }
  };
};

template <typename Topology>
struct moreau_jean_systems<empty_item, Topology> : empty_item {
  using attributes = gather<>;
};

// a kind of interactions, between the systems of the kinds Systems1 (ds1) and
// Systems2 (ds2)
template <typename Interaction, typename Systems1, typename Systems2, typename Topology>
struct moreau_jean_interactions : item {
  using interaction = Interaction;
  using systems1 = Systems1;
  using systems2 = Systems2;
  using topology = Topology;
  using items = gather<Interaction>;

  using nslaw = typename interaction::nslaw;
  using nslaw_size_t = typename interaction::nslaw_size;
  using y = attr_t<interaction, "y">;
  using ydot = attr_t<interaction, "ydot">;
  using lambda = attr_t<interaction, "lambda">;
  using h_matrix1 = attr_t<interaction, "h_matrix1">;
  using h_matrix2 = attr_t<interaction, "h_matrix2">;

  struct attributes {
    some::indice inter_offset;
    some::indice number_of_interactions;
  };

  using properties = gather<storage::keep<y, 2>, storage::keep<ydot, 2>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) inter_offset() { return storage::attr<"inter_offset">(*self()); }
    decltype(auto) number_of_interactions() {
      return storage::attr<"number_of_interactions">(*self());
    }

    static constexpr auto nslaw_size() { return nslaw_size_t{}.value; }

    static constexpr auto nslaw_with_friction() {
      return std::derived_from<nslaw, model::newton_impact_friction>;
    }

    // the context of the kernels of this kind (see simul/kernels.hpp)
    struct kernel_context {
      using interaction = Interaction;
      using kind = moreau_jean_interactions;
      using topology = Topology;
    };

    template <typename Kernel>
    void run(auto step, auto&& p) {
      kernels::run<Kernel, kernel_context>(self()->data(), step, p);
    }

    // the part of an assembled vector that belongs to the interaction of a
    // rank
    static auto part() {
      return [](auto& assembled, auto index, const auto&) {
        return get_vector(assembled, index);
      };
    }

    // the rows of the interactions in the assembled matrices and vectors
    decltype(auto) h_matrix_view(auto& assembled, auto ds_offset) {
      auto storage = get_storage_type(self()->data(), h_matrix2{});
      return algebra::mat_view(storage, assembled, inter_offset(), ds_offset);
    }

    decltype(auto) w_matrix_view(auto& assembled) {
      auto storage =
          convert_storage_type(typename systems1::system{}, self()->data(),
                               some::matrix<some::scalar, nth_t<0, typename h_matrix1::sizes>,
                                            nth_t<0, typename h_matrix1::sizes>>{});
      return algebra::mat_view(storage, assembled, inter_offset(), inter_offset());
    }

    decltype(auto) lambda_view(auto& assembled) {
      auto storage = get_storage_type(self()->data(), lambda{});
      return algebra::vec_view(storage, assembled, inter_offset());
    }

    decltype(auto) ydot_view(auto& assembled) {
      auto storage = get_storage_type(self()->data(), ydot{});
      return algebra::vec_view(storage, assembled, inter_offset());
    }

    decltype(auto) mu_view(auto& assembled) {
      using env_t = decltype(self()->env());
      using vec_mu_t = traits::config<env_t>::template convert<
          some::vector<some::scalar, some::indice_value<1>>>::type;
      return algebra::vec_view<vec_mu_t>(assembled, inter_offset() / 2);
    }

    void initialize(auto step) {
      // useless at initialisation ?
      run<moreau_jean::zero_multipliers>(step, nullptr);

      compute_h_matrices(step);
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

    void compute_h_matrices(auto step) {
      run<moreau_jean::jacobians>(step, moreau_jean::at_step{step});
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

    void compute_output(auto step) {
      run<moreau_jean::output>(step, moreau_jean::at_step{step});
    }

    // the active set: the active interactions, their rank among the active
    // interactions of this kind, the systems they involve; the number of
    // active interactions
    auto mark_active(auto step, auto h, auto gamma, auto threshold) {
      auto& data = self()->data();
      using indice = typename decltype(self()->env())::indice;
      if constexpr (with_index_sets<Topology>) {
        index_set<interaction, active_set_tag, Topology>(data).clear();
      }
      run<kernels::fused<moreau_jean::activation, moreau_jean::mark_involved>>(
          step, moreau_jean::active_set{h, static_cast<decltype(h)>(gamma),
                                        static_cast<decltype(h)>(threshold)});
      auto rows = kernels::offsets<indice>(
          storage::prop_values<interaction, "activation", Topology>(data, step));
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto row) {
          index_set<interaction, active_set_tag, Topology>(data).emplace(i, row);
        };
        run<moreau_jean::active_entry>(step, moreau_jean::ranks{rows, emplace});
      }
      number_of_interactions() = rows.back();
      return rows.back();
    }

    // the blocks of H of the active interactions, the systems being at the
    // columns of their kinds
    void assemble_h_matrix(auto step, auto& assembled, auto ds1_systems, auto ds2_systems) {
      using env_t = decltype(self()->env());
      using matrix_1x1_t = typename env_t::template matrix<typename env_t::scalar, 1, 1>;

      if constexpr (std::same_as<typename systems2::system,
                                 typename topology::fixed_dof_system>) {
        // between rigid systems (i.e disk/disk)
        run<moreau_jean::h_blocks>(step, moreau_jean::block_matrix{h_matrix_view(
                                             assembled, ds2_systems.ds_offset())});
      } else if constexpr (std::same_as<systems1, systems2>) {
        // self contact of a FEM system
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        run<moreau_jean::fem_self_h_block>(
            step, moreau_jean::fem_block_matrix{fem_matrix, ds2_systems.sum_dofs()});
      } else {
        // a rigid system (ds1) and a FEM one (ds2) (i.e disk/fem)
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        auto rigid_storage = get_storage_type(self()->data(), h_matrix1{});
        auto rigid_matrix = algebra::mat_view(rigid_storage, assembled, inter_offset(),
                                              ds1_systems.ds_offset());
        run<moreau_jean::fem_h_blocks>(
            step, moreau_jean::fem_h_matrix{self()->data(), fem_matrix, rigid_matrix,
                                            ds2_systems.sum_dofs(),
                                            static_cast<std::size_t>(nslaw_size())});
      }
    }

    void assemble_vectors(auto step, auto& lambdas, auto& ydots, auto& mus) {
      if constexpr (nslaw_with_friction()) {
        run<kernels::fused<moreau_jean::assemble_multipliers, moreau_jean::assemble_friction>>(
            step, moreau_jean::multipliers{part(), lambda_view(lambdas), ydot_view(ydots),
                                           mu_view(mus)});
      } else {
        run<moreau_jean::assemble_multipliers>(
            step, moreau_jean::multipliers{part(), lambda_view(lambdas), ydot_view(ydots)});
      }
    }

    // nonsmooth law effect
    void nsl_effect_on_free_output(auto step) { run<moreau_jean::impact_law>(step, nullptr); }

    // H v_free: the free output derivatives in q of the problem
    void free_output(auto step, auto& assembled) {
      run<moreau_jean::free_output>(step, moreau_jean::block_vector{lambda_view(assembled)});
    }

    void keep_lambdas(auto step, auto& lambdas) {
      run<moreau_jean::keep_multipliers>(
          step, moreau_jean::multipliers{part(), lambda_view(lambdas)});
    }
  };
};

template <typename Systems1, typename Systems2, typename Topology>
struct moreau_jean_interactions<empty_item, Systems1, Systems2, Topology> : empty_item {
  using attributes = gather<>;
};

}  // namespace siconos::simul
