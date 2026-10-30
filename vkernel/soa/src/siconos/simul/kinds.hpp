#pragma once

#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/contacts.hpp"
#include "siconos/simul/index_sets.hpp"
#include "siconos/simul/relation_sets.hpp"
#include "siconos/storage/get.hpp"
#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/schedule.hpp"
#include "siconos/utils/variant.hpp"

namespace siconos::simul {

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
        run<contacts::dofs>(step, contacts::dof_offsets{sum_dofs()});
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
      run<contacts::not_involved>(step, nullptr);
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
      run<contacts::involved_index>(step, contacts::ranks{offsets});
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto index) {
          index_set<system, involved_set_tag, Topology>(data).emplace(i, index);
        };
        run<contacts::involved_entry>(step, contacts::ranks{offsets, emplace});
      }
      number_of_involved_ds() = offsets.back();
      return offsets.back();
    }

    // f(i, index) for every involved system i, index being its index among
    // the involved systems
    void for_each_involved(auto step, auto&& f) {
      simul::for_each_involved<system, topology>(self()->data(), step, f);
    }

    void assemble_mass_matrix(auto step, auto& assembled) {
      // size may be 0
      run<contacts::diagonal_block<"mass_matrix">>(
          step, contacts::diagonal_blocks{mass_matrix_view(assembled), at()});
    }

    void assemble_k_matrix(auto step, auto& assembled) {
      if constexpr (with_k_matrix()) {
        // size may be 0
        run<contacts::diagonal_block<"k_matrix">>(
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

    // an assembled matrix for the transposed jacobians, xxx properties on
    // h_matrix1, if any, are lost here
    auto transposed_h_matrix() {
      using env = decltype(self()->env());
      return typename traits::config<env>::template convert<
          some::assembled_matrix<some::transposed_matrix<h_matrix1>>>::type{};
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
      run<contacts::jacobians>(step, contacts::at_step{step});
    }

    static constexpr bool time_invariant_h_matrices() {
      using data_t = const std::decay_t<typename Handle::data_t>;
      return storage::has_property_t<interaction, property::time_invariant, data_t>();
    }

    void compute_output(auto step) { run<contacts::output>(step, contacts::at_step{step}); }

    // the active set: the active interactions, their rank among the active
    // interactions of this kind, the systems they involve; the number of
    // active interactions
    auto mark_active(auto step, auto h, auto gamma, auto threshold) {
      auto& data = self()->data();
      using indice = typename decltype(self()->env())::indice;
      if constexpr (with_index_sets<Topology>) {
        index_set<interaction, active_set_tag, Topology>(data).clear();
      }
      run<kernels::fused<contacts::activation, contacts::mark_involved>>(
          step, contacts::active_set{h, static_cast<decltype(h)>(gamma),
                                     static_cast<decltype(h)>(threshold)});
      auto rows = kernels::offsets<indice>(
          storage::prop_values<interaction, "activation", Topology>(data, step));
      if constexpr (with_index_sets<Topology>) {
        auto emplace = [&data](auto i, auto row) {
          index_set<interaction, active_set_tag, Topology>(data).emplace(i, row);
        };
        run<contacts::active_entry>(step, contacts::ranks{rows, emplace});
      }
      number_of_interactions() = rows.back();
      return rows.back();
    }

    // the blocks of H of the active interactions, the systems being at the
    // columns of their kinds
    void assemble_h_matrix(auto step, auto& assembled, auto ds1_systems, auto ds2_systems) {
      using env_t = decltype(self()->env());
      using matrix_1x1_t = typename env_t::template matrix<typename env_t::scalar, 1, 1>;

      if constexpr (between_rigid_systems()) {
        // between rigid systems (i.e disk/disk)
        run<contacts::h_blocks>(
            step, contacts::block_matrix{h_matrix_view(assembled, ds2_systems.ds_offset())});
      } else if constexpr (std::same_as<systems1, systems2>) {
        // self contact of a FEM system
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        run<contacts::fem_self_h_block>(
            step, contacts::fem_block_matrix{fem_matrix, ds2_systems.sum_dofs()});
      } else {
        // a rigid system (ds1) and a FEM one (ds2) (i.e disk/fem)
        auto fem_matrix = algebra::mat_view<matrix_1x1_t>(assembled, inter_offset(),
                                                          ds2_systems.ds_offset());
        auto rigid_storage = get_storage_type(self()->data(), h_matrix1{});
        auto rigid_matrix = algebra::mat_view(rigid_storage, assembled, inter_offset(),
                                              ds1_systems.ds_offset());
        run<contacts::fem_h_blocks>(
            step, contacts::fem_h_matrix{self()->data(), fem_matrix, rigid_matrix,
                                         ds2_systems.sum_dofs(),
                                         static_cast<std::size_t>(nslaw_size())});
      }
    }

    void assemble_vectors(auto step, auto& lambdas, auto& ydots, auto& mus) {
      if constexpr (nslaw_with_friction()) {
        run<kernels::fused<contacts::assemble_multipliers, contacts::assemble_friction>>(
            step, contacts::multipliers{part(), lambda_view(lambdas), ydot_view(ydots),
                                        mu_view(mus)});
      } else {
        run<contacts::assemble_multipliers>(
            step, contacts::multipliers{part(), lambda_view(lambdas), ydot_view(ydots)});
      }
    }

    void keep_lambdas(auto step, auto& lambdas) {
      run<contacts::keep_multipliers>(step,
                                      contacts::multipliers{part(), lambda_view(lambdas)});
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

}  // namespace siconos::simul
