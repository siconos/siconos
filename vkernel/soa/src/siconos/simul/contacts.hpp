#pragma once

#include <cassert>
#include <cstddef>
#include <type_traits>

#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
// after the models
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/simul/kernels.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::contacts {

/* The contacts of a time step, whatever the scheme: the relations (the
   jacobians H, the outputs y and ydot), the active set (the interactions
   whose gap may close, the systems they involve) and the assembly of the
   problem of the step (H, M, K, lambda, ydot, mu of the active
   interactions). The kernels, then the stages on the kinds (see kinds.hpp)
   and on the integrator. */

using namespace kernels;

// the storage of the contacts on the interactions of a kind (the properties
// of the kind, see kinds.hpp): y and ydot at two time levels, the outputs at
// step and the free outputs at step + 1 (the active set reads y, ydot at
// step after the free outputs, the impact law ydot at step and step + 1)
template <typename Interaction>
using interaction_storage = gather<storage::time_levels<attr_t<Interaction, "y">, 2>,
                                   storage::time_levels<attr_t<Interaction, "ydot">, 2>>;

// the parameters of the kernels

// an assembled matrix of FEM systems and the dof offsets of the involved ones
template <typename Matrix, typename Offsets>
struct fem_block_matrix {
  Matrix matrix;
  Offsets& sum_dofs;
};

// the active set: gamma for the gap at the middle of the step, the
// activation threshold
template <typename H, typename S>
struct active_set {
  H h;
  S gamma;
  S threshold;
};

// the ranks of the flagged items of a column, offsets[i] (the scan of the
// flags, see kernels::offsets), and the index set of the flagged items,
// emplace(i, rank)
template <typename Offsets, typename Emplace = std::nullptr_t>
struct ranks {
  Offsets& offsets;
  Emplace emplace = nullptr;
};

// the H matrix of the FEM systems in contact with rigid ones: its columns of
// the FEM systems (matrix), of the rigid systems (rigid_matrix), the dof
// offsets of the involved FEM systems
template <typename Data, typename Matrix, typename RigidMatrix, typename Offsets,
          typename Indice>
struct fem_h_matrix {
  Data& data;
  Matrix matrix;
  RigidMatrix rigid_matrix;
  Offsets& sum_dofs;
  Indice nslaw_size;
};

// the step of the relations
template <typename Step>
struct at_step {
  Step step;
};

// the assembled vectors of the multipliers, of the output derivatives and of
// the friction coefficients, part as in time_step
template <typename Part, typename Lambdas, typename Ydots = std::nullptr_t,
          typename Mus = std::nullptr_t>
struct multipliers {
  Part part;
  Lambdas lambdas;
  Ydots ydots = nullptr;
  Mus mus = nullptr;
};

// an assembled block matrix
template <typename Matrix>
struct block_matrix {
  Matrix matrix;
};

// an assembled block diagonal matrix, the block of the system of a rank at
// at(rank)
template <typename Matrix, typename At>
struct diagonal_blocks {
  Matrix matrix;
  At at;
};

// sum_dofs[k] is the first dof of the involved system of rank k
template <typename Offsets>
struct dof_offsets {
  Offsets& sum_dofs;
};

// the kernels

// a block of the assembled matrix: the attribute Name of an involved system
template <storage::pattern::string_literal Name>
struct assemble_diagonal_block : kernel<involved<systems>, rank, in<Name>>, assembles {
  static void apply(auto& p, auto index, const auto& block) {
    set_value(p.matrix, p.at(index), p.at(index), block);
  }
};

// the number of dofs of an involved system, before the scan that gives the
// offsets (see system_kind::compute_total_dofs)
struct compute_dofs : kernel<involved<systems>, rank, in<"mass_matrix">>, assembles {
  static void apply(auto& p, auto index, const auto& mass_matrix) {
    p.sum_dofs[index + 1] = algebra::ncols(mass_matrix);
  }
};

// the multipliers start from zero
struct zero_multipliers
    : kernel<all<interactions>, out<"lambda">, out<"ydot">, out<"ydot", next>> {
  static void apply(const auto&, auto& lambda, auto& ydot, auto& ydot_next) {
    algebra::set_zero(lambda);
    algebra::set_zero(ydot);
    algebra::set_zero(ydot_next);
  }
};

// an active interaction in the assembled vectors
struct assemble_multipliers : kernel<active<interactions>, rank, in<"lambda">, in<"ydot">>,
                              assembles {
  static void apply(auto& p, auto index, const auto& lambda, const auto& ydot) {
    set_value(p.lambdas, index, lambda);
    set_value(p.ydots, index, ydot);
  }
};

struct assemble_friction : kernel<active<interactions>, rank, through<"nslaw", "mu">>,
                           assembles {
  static void apply(auto& p, auto index, const auto& mu) { set_value(p.mus, index, mu); }
};

// the multipliers of an active interaction from the solution
struct keep_multipliers : kernel<active<interactions>, rank, out<"lambda">> {
  static void apply(auto& p, auto index, auto& lambda) {
    lambda = p.part(p.lambdas, index, lambda);
  }
};

// the relations: their computations go through handles (host only), and
// those of the nonlinear relations through narrowphase caches shared by all
// the relations of a step (reads_world)

// the jacobians of an interaction
struct compute_jacobians
    : kernel<relations<interactions>, relation, referenced<"ds1", by_topology>,
             referenced<"ds2", by_topology>, inout<"h_matrix1">, inout<"h_matrix2">>,
      reads_world,
      host_only {
  static void apply(const auto& p, auto& rel, auto ds1, auto ds2, auto& h1, auto& h2) {
    using rel_t = std::decay_t<decltype(rel)>;
    // the first of relation1, relation2 wins: a linear relation that is a
    // relation1 (lagrangian_r) gets the one ds jacobian whatever its nds
    if constexpr (match::relation1<rel_t>) {
      rel.compute_jachq(p.step, ds1, h1);
    } else if constexpr (match::relation2<rel_t>) {
      rel.compute_jachq(p.step, ds1, ds2, h1, h2);
    } else {
      assert(false);
    }
  }
};

// y and ydot of an interaction from the state of its systems
struct compute_output : kernel<relations<interactions>, relation, in<"nds", by_topology>,
                               referenced<"ds1", by_topology>, referenced<"ds2", by_topology>,
                               in<"h_matrix1">, in<"h_matrix2">, out<"y">, out<"ydot">>,
                        reads_world,
                        host_only {
  static void apply(const auto& p, auto& rel, auto nds, auto ds1, auto ds2, const auto& h1,
                    const auto& h2, auto& y, auto& ydot) {
    using rel_t = std::decay_t<decltype(rel)>;
    const auto step = p.step;

    // y[1] (the tangential gap row) is only written by the linear relations,
    // the nonlinear ones set y[0] from compute_h(): io::contact_points must
    // not read an uninitialized y[1] (HDF5 column 15)
    algebra::set_zero(y);

    // the first that matches wins (as mp::overload)
    if constexpr (match::linear_relation<rel_t>) {
      y = h1 * storage::attr<"q">(ds1, step) + rel.b();
      ydot = h1 * storage::attr<"velocity">(ds1, step);
      if (nds == 2) {
        y += h2 * storage::attr<"q">(ds2, step);
        ydot += h2 * storage::attr<"velocity">(ds2, step);
      }
    } else if constexpr (match::handle<rel_t, collision::mesh_relation>) {
      using env = decltype(ds2.env());
      using indice = typename env::indice;
      using vector_4d = typename env::template vector<typename env::scalar, 4>;

      y[0] = rel.compute_h(step, ds1, ds2);

      // the nodes of the contact segment
      indice k = 4 * rel.contact_index();
      const auto& nodes = rel.mesh().global_indices();
      const auto& v2 = storage::attr<"velocity">(ds2, step);
      vector_4d local_velocity = {v2[nodes[k]], v2[nodes[k + 1]], v2[nodes[k + 2]],
                                  v2[nodes[k + 3]]};

      ydot = h1 * storage::attr<"velocity">(ds1, step);
      ydot += h2 * local_velocity;
    } else if constexpr (match::relation1<rel_t>) {
      y[0] = rel.compute_h(step, ds1);
      ydot = h1 * storage::attr<"velocity">(ds1, step);
    } else if constexpr (match::relation2<rel_t>) {
      y[0] = rel.compute_h(step, ds1, ds2);
      ydot = h1 * storage::attr<"velocity">(ds1, step);
      ydot += h2 * storage::attr<"velocity">(ds2, step);
    }
  }
};

// the active set

// before the active set: no system is involved
struct clear_involved : kernel<all<systems>, out<"involved", by_topology>> {
  static void apply(const auto&, auto& involved) { involved = false; }
};

// an interaction is active when its gap at the middle of the step may be
// below the threshold
struct compute_activation
    : kernel<all<interactions>, in<"y">, in<"ydot">, out<"activation", by_topology>> {
  static void apply(const auto& p, const auto& y, const auto& ydot, auto& activation) {
    activation = ((y + p.gamma * p.h * ydot)(0) <= p.threshold);
  }
};

// the systems of an active interaction are involved: several interactions
// may write the same system, they all write true
struct mark_involved
    : kernel<all<interactions>, in<"activation", by_topology>, in<"nds", by_topology>,
             out_through<"ds1", "involved", by_topology>,
             out_through<"ds2", "involved", by_topology>> {
  static void apply(const auto&, bool activation, auto nds, auto& involved1, auto& involved2) {
    if (activation) {
      involved2 = true;
      if (nds == 2) {
        involved1 = true;
      }
    }
  }
};

// the index of an involved system: its rank among the involved systems of
// its type, in the order of the systems
struct set_involved_index : kernel<all<systems>, item_index, in<"involved", by_topology>,
                                   out<"index", by_topology>> {
  static void apply(const auto& p, auto i, bool involved, auto& index) {
    if (involved) {
      index = p.offsets[i];
    }
  }
};

// the active interactions and the involved systems in their index sets (host
// containers, filled in the order of the items)
struct collect_active_interactions
    : kernel<all<interactions>, item_index, in<"activation", by_topology>>,
      host_only,
      ordered {
  static void apply(auto& p, auto i, bool activation) {
    if (activation) {
      p.emplace(i, p.offsets[i]);
    }
  }
};

struct collect_involved_systems
    : kernel<all<systems>, item_index, in<"involved", by_topology>, in<"index", by_topology>>,
      host_only,
      ordered {
  static void apply(auto& p, auto i, bool involved, auto index) {
    if (involved) {
      p.emplace(i, index);
    }
  }
};

// the blocks of H of an active interaction between rigid systems: its
// jacobians without the columns of the blocked velocity components, at the
// rows of the interaction and the columns of its systems
struct assemble_h_blocks
    : kernel<active<interactions>, rank, in<"h_matrix1">, in<"h_matrix2">,
             through<"ds1", "index", by_topology>, through<"ds2", "index", by_topology>,
             through<"ds1", "bc_velocities_0", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>>,
      assembles {
  static void apply(auto& p, auto row, const auto& h1, const auto& h2, auto j1, auto j2,
                    const auto& blocked1, const auto& blocked2) {
    auto h1_mod = h1;
    auto h2_mod = h2;
    for (auto k : blocked1) {
      h1_mod.col(k).setZero();
    }
    if (j1 != j2) {
      for (auto k : blocked2) {
        h2_mod.col(k).setZero();
      }
      set_value(p.matrix, row, j1, h1_mod);
      set_value(p.matrix, row, j2, h2_mod);
    } else {
      // self interaction
      set_value(p.matrix, row, j2, h1_mod);
    }
  }
};

// the blocks of H of an active interaction between a rigid system (ds1) and
// a FEM one (ds2): the rigid block, then the FEM columns of the nodes of the
// contact segment of its mesh relation
struct assemble_fem_h_blocks
    : kernel<active<interactions>, rank, in<"h_matrix1">, in<"h_matrix2">, in<"relation">,
             through<"ds1", "index", by_topology>, through<"ds2", "index", by_topology>,
             through<"ds1", "bc_velocities_0", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>>,
      host_only,
      assembles {
  static void apply(auto& p, auto row, const auto& h1, const auto& h2, const auto& rel,
                    auto j1, auto i2, const auto& blocked1, const auto& blocked2) {
    const auto j2 = p.sum_dofs[i2];
    const auto rows = row * p.nslaw_size;

    auto h1_mod = h1;
    for (auto k : blocked1) {
      h1_mod.col(k).setZero();
    }
    set_value(p.rigid_matrix, row, j1, h1_mod);

    variant::visit(p.data, rel,
                   mp::overload(
                       [&](match::handle<collision::mesh_relation> auto rrel) {
                         const auto& nodes = rrel.mesh().global_indices();
                         const auto base = 4 * rrel.contact_index();
                         for (auto i = 0; i < algebra::nrows(h2); ++i) {
                           for (auto c = 0; c < 4; ++c) {
                             set_value(p.matrix, i + rows, nodes[base + c] + j2, h2(i, c));
                           }
                         }
                       },
                       [](auto) {}));

    // the blocked velocity components of ds2
    for (auto i = 0; i < algebra::nrows(h2); ++i) {
      for (auto k : blocked2) {
        set_value(p.matrix, i + rows, j2 + k, 0.);
      }
    }
  }
};

// the block of H of an active self contact of a FEM system (ds1 == ds2)
struct assemble_fem_self_h_block
    : kernel<active<interactions>, rank, in<"h_matrix1">, in<"nds", by_topology>,
             through<"ds2", "index", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>>,
      assembles {
  static void apply(auto& p, auto row, const auto& h1, [[maybe_unused]] auto nds, auto i2,
                    const auto& blocked2) {
    // only ds1 == ds2 (nds == 1) is supported (rt_lagrangian_r is a relation1)
    assert(nds == 1);
    const auto j2 = p.sum_dofs[i2];
    auto h1_mod = h1;
    for (auto k : blocked2) {
      h1_mod.col(k).setZero();
    }
    set_value(p.matrix, row, j2, h1_mod);
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::condition;
using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct set_initial_state : stage<"Set lambda, ydot to zero and compute the jacobians H"> {
  static void run(auto& c) {
    c.kind.template run<zero_multipliers>(c.step, nullptr);
    c.kind.compute_h_matrices(c.step);
  }
};

struct sync_relation_sets
    : stage<"Synchronize the relation sets: they hold every interaction">,
      host_code {
  static void run(auto& c) { c.kind.sync_relation_sets(c.step); }
};

struct compute_jacobians : stage<"Compute the jacobians H, if the relations change"> {
  static void run(auto& c) {
    if constexpr (!std::decay_t<decltype(c.kind)>::time_invariant_h_matrices()) {
      c.kind.compute_h_matrices(c.step);
    }
  }
};

struct compute_outputs : stage<"Compute the outputs y, ydot at step"> {
  static void run(auto& c) { c.kind.compute_output(c.step); }
};

struct compute_free_outputs : stage<"Compute the free outputs y, ydot at step + 1"> {
  static void run(auto& c) { c.kind.compute_output(c.step + 1); }
};

// jacobians for the assembly: already up to date at this step if the
// jacobians stage computed them (q(step) has not changed since)
struct compute_time_invariant_jacobians
    : stage<"Compute the jacobians H, if the relations do not change"> {
  static void run(auto& c) {
    if constexpr (std::decay_t<decltype(c.kind)>::time_invariant_h_matrices()) {
      c.kind.compute_h_matrices(c.step);
    }
  }
};

struct assemble_h_matrix : stage<"Assemble H of the active interactions"> {
  static void run(auto& c) { c.integrator.assemble_h_matrix(c.step, c.kind); }
};

struct assemble_multipliers : stage<"Assemble lambda, ydot, mu of the active interactions"> {
  static void run(auto& c) {
    c.kind.assemble_vectors(c.step, c.integrator.lambda_vector_assembled(),
                            c.integrator.ydot_vector_assembled(),
                            c.integrator.mu_vector_assembled());
  }
};

struct keep_multipliers : stage<"Keep lambda of the active interactions from the solution"> {
  static void run(auto& c) {
    c.kind.keep_lambdas(c.step, c.integrator.lambda_vector_assembled());
  }
};

// on a kind of systems

struct assemble_mass_matrix : stage<"Assemble M of the involved systems"> {
  static void run(auto& c) {
    c.kind.assemble_mass_matrix(c.step, c.integrator.mass_matrix_assembled());
  }
};

struct assemble_stiffness_matrix : stage<"Assemble K of the involved FEM systems"> {
  static void run(auto& c) {
    c.kind.assemble_k_matrix(c.step, c.integrator.k_matrix_assembled());
  }
};

// on the integrator

struct compute_active_set : stage<
                                "Compute the active set: the interactions whose gap may "
                                "close, the systems they involve">,
                            host_code {
  static void run(auto& c) {
    c.involved = c.integrator.compute_active_interactions(c.step, c.h);
  }
};

struct with_contacts : condition<"there are contacts"> {
  static bool holds(const auto& c) { return c.involved > 0; }
};

struct setup_assembly
    : stage<"Set up the assembly: the sizes and offsets of the assembled problem">,
      host_code {
  static void run(auto& c) { c.integrator.assemble_setup(c.step); }
};

}  // namespace stages

}  // namespace siconos::simul::contacts
