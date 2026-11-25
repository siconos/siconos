// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#pragma once

#include "siconos/vkernel/algebra/numerics_backend.hpp"
#include "siconos/vkernel/collision/diskmesh_r.hpp"
#include "siconos/vkernel/model/lagrangian_ds.hpp"
#include "siconos/vkernel/model/lagrangian_r.hpp"
#include "siconos/vkernel/model/nslaws.hpp"
#include "siconos/vkernel/simul/contacts.hpp"
#include "siconos/vkernel/simul/index_sets.hpp"
#include "siconos/vkernel/simul/relation_sets.hpp"
#include "siconos/vkernel/storage/get.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"
#include "siconos/vkernel/storage/schedule.hpp"
#include "siconos/vkernel/storage/sparse_set.hpp"
#include "siconos/vkernel/utils/variant.hpp"

namespace siconos::vkernel::simul {

/* The kinds of an integrator: a kind of systems (rigid with a fixed number
   of dofs, FEM with a dynamic one) or a kind of interactions (between rigid
   systems, rigid-FEM, FEM-FEM). A kind runs kernels on its items and keeps
   where they are in the assembled problem; the assembled matrices and
   vectors are those of the integrator, given to the kind. The stages of a
   time step run for each kind (for_systems, for_interactions, at the end).

   Topology is the item that declares the attached storages of the systems
   and interactions (involved, index, ds1, ...): their lookups name it as
   declarer. Storage is what the scheme of the integrator needs on the items
   of the kind (time levels, attached storages: system_storage and
   interaction_storage of a scheme, see moreau_jean.hpp): the properties of
   the kind, which declares them. */

// a kind of systems
template <typename System, typename Topology, typename Storage>
struct system_kind : item {
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

  using properties = Storage;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    static constexpr bool runtime_dof() {
      using data_t = typename Handle::data_t;
      return model::runtime_dof<decltype(storage::make_handle<system>(
          std::declval<data_t>(), storage::index<system, int>{0}))>();
    }

    static constexpr auto with_k_matrix() { return model::has_k_matrix(system{}); }

    // f constant over the simulation
    static constexpr bool time_invariant_fext() {
      using data_t = const std::decay_t<typename Handle::data_t>;
      return storage::has_property_t<fext, property::time_invariant, data_t>();
    }

    decltype(auto) sum_dofs() { return storage::attr<"sum_dofs">(*self()); }
    decltype(auto) ds_offset() { return storage::attr<"ds_offset">(*self()); }
    decltype(auto) number_of_involved_ds() {
      return storage::attr<"number_of_involved_ds">(*self());
    }

    // the context of the kernels of this kind (see simul/kernels.hpp)
    struct kernel_context {
      using system = System;
      using kind = system_kind;
      using topology = Topology;
    };

    template <typename Kernel>
    void run(auto step, auto&& p) {
      kernels::run<Kernel, kernel_context>(self()->data(), step, p);
    }

    // where the involved system of a rank is in the assembled vectors and
    // matrices of this kind, in their blocks: its rank when the systems have
    // the same dofs, its first dof otherwise (the blocks are scalars)
    auto at() {
      if constexpr (runtime_dof()) {
        return [&sum_dofs = sum_dofs()](auto index) { return sum_dofs[index]; };
      } else {
        return [](auto index) { return index; };
      }
    }

    // the part of an assembled vector that belongs to the system of a rank
    auto part() {
      return [at = at()](auto& assembled, auto index, const auto& like) {
        if constexpr (match::fixed_size_vector<velocity>) {
          return get_vector(assembled, at(index));
        } else {
          return get_vector(assembled, at(index), like.size());
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

    // the dofs of the involved systems
    void compute_total_dofs(auto step) {
      if constexpr (runtime_dof()) {
        // runtime degrees of freedom: sum_dofs()[k] is the first dof of the
        // involved system of rank k: the dofs of each system at k + 1, then
        // a scan
        using involved_systems = kernels::involved<kernels::systems>;
        auto n = kernels::count<involved_systems, kernel_context>(self()->data(), step);
        sum_dofs().assign(n + 1, 0);
        run<contacts::compute_dofs>(step, contacts::dof_offsets{sum_dofs()});
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
      run<contacts::clear_involved>(step, nullptr);
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
      run<contacts::set_involved_index>(step, contacts::ranks{offsets});
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto index) {
          index_set<system, involved_set_tag, Topology>(data).emplace(i, index);
        };
        run<contacts::collect_involved_systems>(step, contacts::ranks{offsets, emplace});
      }
      number_of_involved_ds() = offsets.back();
      return offsets.back();
    }

    // f(i, index) for every involved system i, index being its index among
    // the involved systems
    void for_each_involved(auto step, auto&& f) {
      simul::for_each_involved<system, topology>(self()->data(), step, f);
    }

    // the FEM blocks with their Dirichlet conditions (the rigid ones: a
    // diagonal M and zero columns of H suffice for the blocked components)
    void assemble_mass_matrix(auto step, auto& assembled) {
      // size may be 0
      if constexpr (with_k_matrix()) {
        run<contacts::assemble_constrained_diagonal_block<"mass_matrix", 1>>(
            step, contacts::diagonal_blocks{mass_matrix_view(assembled), at()});
      } else {
        run<contacts::assemble_diagonal_block<"mass_matrix">>(
            step, contacts::diagonal_blocks{mass_matrix_view(assembled), at()});
      }
    }

    void assemble_k_matrix(auto step, auto& assembled) {
      if constexpr (with_k_matrix()) {
        // size may be 0
        run<contacts::assemble_constrained_diagonal_block<"k_matrix", 0>>(
            step, contacts::diagonal_blocks{k_matrix_view(assembled), at()});
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

    // the blocked velocity components of the system i (its Dirichlet
    // conditions), f(blocked); nothing if it has none
    void with_blocked_dofs(auto step, std::size_t i, auto&& f) {
      auto& values =
          storage::prop_values<system, "bc_velocities_0", topology>(self()->data(), step);
      if constexpr (storage::is_sparse_set_v<std::decay_t<decltype(values)>>) {
        if (const auto* blocked = values.try_get(i)) f(*blocked);
      } else {
        f(values[i]);
      }
    }

    // M + h^2 theta^2 K, refactorized only if h or theta change; with the
    // Dirichlet conditions of the system (algebra::constrain_dofs)
    algebra::sparse_factor& fem_w_factor(auto step, std::size_t i, auto h, double theta) {
      auto& op = fem_operator(step, i);
      if (!op.w || op.h != h || op.theta != theta) {
        auto& data = self()->data();
        auto& m = storage::attr_values<system, "mass_matrix">(data, step)[i];
        auto& k = storage::attr_values<system, "k_matrix">(data, step)[i];
        algebra::sparse_factor::matrix_t w = m;
        w += (h * h * theta * theta) * k;
        with_blocked_dofs(
            step, i, [&](const auto& blocked) { algebra::constrain_dofs(w, blocked, 1.); });
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
        with_blocked_dofs(
            step, i, [&](const auto& blocked) { algebra::constrain_dofs(m, blocked, 1.); });
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

template <typename Topology, typename Storage>
struct system_kind<empty_item, Topology, Storage> : empty_item {
  using attributes = gather<>;
};

// a kind of interactions, between the systems of the kinds Systems1 (ds1) and
// Systems2 (ds2)
template <typename Interaction, typename Systems1, typename Systems2, typename Topology,
          typename Storage>
struct interaction_kind : item {
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

  using properties = Storage;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) inter_offset() { return storage::attr<"inter_offset">(*self()); }
    decltype(auto) number_of_interactions() {
      return storage::attr<"number_of_interactions">(*self());
    }

    static constexpr auto nslaw_size() { return nslaw_size_t{}.value; }

    // between rigid systems (i.e disk/disk), not FEM ones
    static constexpr bool between_rigid_systems() {
      return std::same_as<typename systems2::system, typename topology::fixed_dof_system>;
    }

    static constexpr auto nslaw_with_friction() {
      return std::derived_from<nslaw, model::newton_impact_friction>;
    }

    // the rolling rows, rolling resistance mu_r (see
    // model::newton_impact_rolling_friction)
    static constexpr auto nslaw_with_rolling_friction() {
      return std::derived_from<nslaw, model::newton_impact_rolling_friction>;
    }

    static constexpr int rolling_rows() {
      if constexpr (nslaw_with_rolling_friction()) {
        return nslaw::rolling_size;
      } else {
        return 0;
      }
    }

    // f(rows), rows (a std::integral_constant) the rows of a contact of this
    // kind in the assembled problem: those of the law, or those of the
    // frictional contact problem when no active contact resists rolling
    // (rolling false), the rolling rows left out
    static void with_assembled_rows(bool rolling, auto&& f) {
      constexpr int size = nslaw_size_t::value;
      if constexpr (rolling_rows() > 0) {
        if (!rolling) {
          f(std::integral_constant<int, size - rolling_rows()>{});
          return;
        }
      }
      f(std::integral_constant<int, size>{});
    }

    // Coulomb friction on u(k + theta), see model::fremond_impact_friction
    static constexpr auto nslaw_is_fremond() {
      return std::derived_from<nslaw, model::fremond_impact_friction>;
    }

    // the context of the kernels of this kind (see simul/kernels.hpp)
    struct kernel_context {
      using interaction = Interaction;
      using kind = interaction_kind;
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

    // the rows of the interactions in the assembled matrices and vectors,
    // Rows rows per interaction (see with_assembled_rows)

    // the jacobians of Rows rows (H for a lazy instantiation: the jacobian of
    // a FEM self contact has no fixed number of columns)
    template <int Rows, typename H = h_matrix1>
    using rows_of_h_matrix =
        some::matrix<some::scalar, some::indice_value<Rows>, nth_t<1, typename H::sizes>>;

    template <int Rows = nslaw_size_t::value>
    decltype(auto) h_matrix_view(auto& assembled, auto ds_offset,
                                 std::integral_constant<int, Rows> = {}) {
      if constexpr (Rows == nslaw_size_t::value) {
        auto storage = get_storage_type(self()->data(), h_matrix1{});
        return algebra::mat_view(storage, assembled, inter_offset(), ds_offset);
      } else {
        auto storage = convert_storage_type(typename systems1::system{}, self()->data(),
                                            rows_of_h_matrix<Rows>{});
        return algebra::mat_view(storage, assembled, inter_offset(), ds_offset);
      }
    }

    template <int Rows = nslaw_size_t::value>
    decltype(auto) w_matrix_view(auto& assembled, std::integral_constant<int, Rows> = {}) {
      auto storage = convert_storage_type(
          typename systems1::system{}, self()->data(),
          some::matrix<some::scalar, some::indice_value<Rows>, some::indice_value<Rows>>{});
      return algebra::mat_view(storage, assembled, inter_offset(), inter_offset());
    }

    // an assembled matrix for the transposed jacobians, xxx properties on
    // h_matrix1, if any, are lost here
    template <int Rows = nslaw_size_t::value>
    auto transposed_h_matrix(std::integral_constant<int, Rows> = {}) {
      using env = decltype(self()->env());
      if constexpr (Rows == nslaw_size_t::value) {
        return typename traits::config<env>::template convert<
            some::assembled_matrix<some::transposed_matrix<h_matrix1>>>::type{};
      } else {
        return typename traits::config<env>::template convert<
            some::assembled_matrix<some::transposed_matrix<rows_of_h_matrix<Rows>>>>::type{};
      }
    }

    // the multipliers, the output derivatives (and the free ones, the impact
    // terms) of Rows rows
    template <int Rows = nslaw_size_t::value>
    decltype(auto) lambda_view(auto& assembled, std::integral_constant<int, Rows> = {}) {
      if constexpr (Rows == nslaw_size_t::value) {
        auto storage = get_storage_type(self()->data(), lambda{});
        return algebra::vec_view(storage, assembled, inter_offset());
      } else {
        using env_t = decltype(self()->env());
        using vec_t = traits::config<env_t>::template convert<
            some::vector<some::scalar, some::indice_value<Rows>>>::type;
        return algebra::vec_view<vec_t>(assembled, inter_offset());
      }
    }

    template <int Rows = nslaw_size_t::value>
    decltype(auto) ydot_view(auto& assembled, std::integral_constant<int, Rows> rows = {}) {
      if constexpr (Rows == nslaw_size_t::value) {
        auto storage = get_storage_type(self()->data(), ydot{});
        return algebra::vec_view(storage, assembled, inter_offset());
      } else {
        return lambda_view(assembled, rows);
      }
    }

    // a coefficient per interaction: mu, mu_r
    template <int Rows = nslaw_size_t::value>
    decltype(auto) mu_view(auto& assembled, std::integral_constant<int, Rows> = {}) {
      using env_t = decltype(self()->env());
      using vec_mu_t = traits::config<env_t>::template convert<
          some::vector<some::scalar, some::indice_value<1>>>::type;
      return algebra::vec_view<vec_mu_t>(assembled, inter_offset() / Rows);
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
      run<contacts::compute_jacobians>(
          step, contacts::at_step<decltype(step), rolling_rows()>{step});
    }

    static constexpr bool time_invariant_h_matrices() {
      using data_t = const std::decay_t<typename Handle::data_t>;
      return storage::has_property_t<interaction, property::time_invariant, data_t>();
    }

    void compute_output(auto step) {
      run<contacts::compute_output>(step, contacts::at_step{step});
    }

    void compute_output_velocity(auto step) {
      run<contacts::compute_output_velocity>(step, contacts::at_step{step});
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
      run<kernels::fused<contacts::compute_activation, contacts::mark_involved>>(
          step, contacts::active_set{h, static_cast<decltype(h)>(gamma),
                                     static_cast<decltype(h)>(threshold)});
      auto rows = kernels::offsets<indice>(
          storage::prop_values<interaction, "activation", Topology>(data, step));
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto row) {
          index_set<interaction, active_set_tag, Topology>(data).emplace(i, row);
        };
        run<contacts::collect_active_interactions>(step, contacts::ranks{rows, emplace});
      }
      number_of_interactions() = rows.back();
      return rows.back();
    }

    // does an active interaction resist rolling (see
    // contacts::mark_rolling_resistance)
    bool resists_rolling(auto step) {
      if constexpr (rolling_rows() > 0) {
        contacts::rolling_resistances<rolling_rows()> p;
        run<contacts::mark_rolling_resistance>(step, p);
        return p.found;
      } else {
        return false;
      }
    }

    // the blocks of H of the active interactions, the systems being at the
    // columns of their kinds
    template <int Rows>
    void assemble_h_matrix(auto step, auto& assembled, auto ds1_systems, auto ds2_systems,
                           std::integral_constant<int, Rows> rows) {
      using env_t = decltype(self()->env());
      using matrix_1x1_t = typename env_t::template matrix<typename env_t::scalar, 1, 1>;

      if constexpr (between_rigid_systems()) {
        // between rigid systems (i.e disk/disk)
        run<contacts::assemble_h_blocks>(step, contacts::block_matrix{h_matrix_view(
                                                   assembled, ds2_systems.ds_offset(), rows)});
      } else if constexpr (std::same_as<systems1, systems2>) {
        // self contact of a FEM system
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        run<contacts::assemble_fem_self_h_block>(
            step, contacts::fem_block_matrix{fem_matrix, ds2_systems.sum_dofs(), rows});
      } else {
        // a rigid system (ds1) and a FEM one (ds2) (i.e disk/fem)
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        auto rigid_matrix = h_matrix_view(assembled, ds1_systems.ds_offset(), rows);
        run<contacts::assemble_fem_h_blocks>(
            step, contacts::fem_h_matrix{self()->data(), fem_matrix, rigid_matrix,
                                         ds2_systems.sum_dofs()});
      }
    }

    // lambda, ydot, mu, and mu_r if the rolling rows are assembled
    template <int Rows>
    void assemble_vectors(auto step, auto& lambdas, auto& ydots, auto& mus, auto& mu_rs,
                          std::integral_constant<int, Rows> rows) {
      if constexpr (nslaw_with_friction()) {
        run<kernels::fused<contacts::assemble_multipliers, contacts::assemble_friction>>(
            step, contacts::multipliers{part(), lambda_view(lambdas, rows),
                                        ydot_view(ydots, rows), mu_view(mus, rows)});
      } else {
        run<contacts::assemble_multipliers>(
            step,
            contacts::multipliers{part(), lambda_view(lambdas, rows), ydot_view(ydots, rows)});
      }
      if constexpr (rolling_rows() > 0 && Rows == nslaw_size_t::value) {
        run<contacts::assemble_rolling_resistance>(
            step,
            contacts::rolling_resistances<rolling_rows(), decltype(mu_view(mu_rs, rows))>{
                mu_view(mu_rs, rows)});
      }
    }

    template <int Rows>
    void keep_lambdas(auto step, auto& lambdas, std::integral_constant<int, Rows> rows) {
      run<contacts::keep_multipliers>(
          step, contacts::multipliers{part(), lambda_view(lambdas, rows)});
    }
  };
};

template <typename Systems1, typename Systems2, typename Topology, typename Storage>
struct interaction_kind<empty_item, Systems1, Systems2, Topology, Storage> : empty_item {
  using attributes = gather<>;
};

/* The time step of a simulation as a schedule (see storage/schedule.hpp).

   The context of a step: the time stepping, its integrator, the step, the
   time step h, the number of involved systems once the active set is
   known, and for a solver written as a schedule (repeat) the number of its
   iterations and the change of the last one. for_systems (for_interactions) runs its stages
   for each kind of systems (of interactions) of the integrator, with ctx.kind the kind. */

template <typename TimeStepping, typename Integrator, typename Step, typename H>
struct step_context {
  TimeStepping time_stepping;
  Integrator integrator;
  Step step;
  H h;
  Step involved = 0;
  Step iterations = 0;
  H change = 0;
};

template <typename Ctx, typename Kind>
struct kind_context : Ctx {
  Kind kind;
};

// stages run for each kind, each kind after the other
void for_each_kind(auto& ctx, auto&& kinds, auto&& run_stages) {
  storage::mp::for_each(kinds, [&](auto kind) {
    kind_context<std::decay_t<decltype(ctx)>, decltype(kind)> c{ctx, kind};
    run_stages(c);
    ctx.involved = c.involved;
    ctx.iterations = c.iterations;
    ctx.change = c.change;
  });
}

template <typename... Stages>
struct for_systems {
  static constexpr std::string_view name{"for each kind of systems"};
  using stages = std::tuple<Stages...>;
  static void run(auto& ctx, auto&& run_stages) {
    for_each_kind(ctx, ctx.integrator.systems(), run_stages);
  }
};

template <typename... Stages>
struct for_interactions {
  static constexpr std::string_view name{"for each kind of interactions"};
  using stages = std::tuple<Stages...>;
  static void run(auto& ctx, auto&& run_stages) {
    for_each_kind(ctx, ctx.integrator.interactions(), run_stages);
  }
};

}  // namespace siconos::vkernel::simul
