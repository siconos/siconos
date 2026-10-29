#pragma once

#include <cassert>
#include <cstddef>
#include <type_traits>

#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/collision/diskmesh_r.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/simul/kernels.hpp"

namespace siconos::simul::moreau_jean {

/* The Moreau-Jean time integration as kernels (see storage/kernels.hpp and
   simul/kernels.hpp). The parameters p of a kernel are one of the structs
   below, given by the kinds (moreau_jean_kinds.hpp) when they run the kernel. */

using namespace kernels;

// the parameters of the kernels

// the time step, the theta of the theta-method, and for the contact impulse:
// the assembled velocity correction and the part of it that belongs to a
// system, part(assembled, rank, like)
template <typename H, typename S, typename Part = std::nullptr_t,
          typename Velocities = std::nullptr_t>
struct time_step {
  H h;
  S theta;
  Part part = nullptr;
  Velocities velocities = nullptr;
};

// the time step, the theta of the theta-method and the factorizations of
// the FEM systems kept on them: factor(i) of the system i
template <typename H, typename S, typename Factor>
struct fem_time_step {
  H h;
  S theta;
  Factor factor;
};

// the factorizations of the FEM systems: factor(i) of the system i
template <typename Factor>
struct fem_factors {
  Factor factor;
};

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

// an assembled vector
template <typename Vector>
struct block_vector {
  Vector vector;
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

// an assembled block diagonal matrix
template <typename Matrix>
struct block_matrix {
  Matrix matrix;
};

// sum_dofs[k] is the first dof of the involved system of rank k
template <typename Offsets>
struct dof_offsets {
  Offsets& sum_dofs;
};

// the systems

// theta-method on the positions: q(step + 1) from the velocities at step and
// step + 1, on the systems of Space
template <typename Space>
struct position_on
    : kernel<Space, in<"q">, in<"velocity">, in<"velocity", next>, out<"q", next>> {
  static void apply(const auto& p, const auto& q, const auto& v, const auto& v_next,
                    auto& q_next) {
    q_next = q + p.h * (p.theta * v_next + (1.0 - p.theta) * v);
  }
};

using position = position_on<all<systems>>;
using involved_position = position_on<involved<systems>>;

// M^-1 f of a rigid system, f being constant over the step
struct minv_f : kernel<all<systems>, in<"mass_matrix">, in<"fext">, out<"minv_f", by_kind>,
                       out<"minv_f", by_kind, next>> {
  static void apply(const auto&, auto& mass_matrix, const auto& f, auto& minv_f,
                    auto& minv_f_next) {
    minv_f = f;
    algebra::solve_in_place(mass_matrix, minv_f);
    minv_f_next = minv_f;
  }
};

// free flight of a rigid system (K = 0)
struct free_velocity : kernel<all<systems>, in<"velocity">, in<"minv_f", by_kind>,
                              in<"minv_f", by_kind, next>, out<"velocity", next>> {
  static void apply(const auto& p, const auto& v, const auto& minv_f, const auto& minv_f_next,
                    auto& v_next) {
    v_next = v + p.h * p.theta * minv_f_next + p.h * (1 - p.theta) * minv_f;
  }
};

// M^-1 (f - K q) of a FEM system, with the factorization of M kept on the
// system, f being constant over the step
template <typename Vector>
struct fem_minv_f : kernel<all<systems>, item_index, in<"k_matrix">, in<"q">, in<"fext">,
                           out<"minv_f", by_kind>, out<"minv_f", by_kind, next>>,
                    host_only {
  static void apply(auto& p, auto i, const auto& k, const auto& q, const auto& f, auto& minv_f,
                    auto& minv_f_next) {
    Vector fmkq = f - k * q;
    minv_f = p.factor(i).solve(fmkq);
    minv_f_next = minv_f;
  }
};

// free flight of a FEM system: W dv = h theta F(step + 1) + h (1 - theta)
// F(step), F = f - K q, W = M + h^2 theta^2 K factorized once on the system
template <typename Vector>
struct fem_free_velocity : kernel<all<systems>, item_index, in<"k_matrix">, in<"q">,
                                  in<"velocity">, in<"fext">, out<"velocity", next>>,
                           host_only {
  static void apply(auto& p, auto i, const auto& k, const auto& q, const auto& v,
                    const auto& f, auto& v_next) {
    Vector tf = f - k * q;
    Vector tf_next = f - k * (q + p.h * v);
    Vector rhs = p.h * p.theta * tf_next + p.h * (1 - p.theta) * tf;
    Vector dv = p.factor(i).solve(rhs);
    v_next = v + dv;
  }
};

// the contact impulse of an involved system
struct contact_velocity : kernel<involved<systems>, rank, inout<"velocity", next>> {
  static void apply(auto& p, auto index, auto& v_next) {
    v_next += p.part(p.velocities, index, v_next);
  }
};

// the boundary conditions of a system: its blocked velocity components are
// zero, and its position follows. Only the systems with boundary conditions
// are visited when bc_velocities_0 is a sparse set.
struct boundary_conditions
    : kernel<entries<systems, "bc_velocities_0", by_topology>, entry, in<"q">, in<"velocity">,
             inout<"velocity", next>, out<"q", next>> {
  static void apply(const auto& p, const auto& blocked, const auto& q, const auto& v,
                    auto& v_next, auto& q_next) {
    if (blocked.size() == 0) return;
    for (auto k : blocked) {
      v_next(k) = 0.0;
    }
    position::apply(p, q, v, v_next, q_next);
  }
};

// a block of the assembled matrix: the attribute Name of an involved system
template <storage::pattern::string_literal Name>
struct diagonal_block : kernel<involved<systems>, rank, in<Name>>, assembles {
  static void apply(auto& p, auto index, const auto& block) {
    set_value(p.matrix, index, index, block);
  }
};

// the number of dofs of an involved system, before the scan that gives the
// offsets (see moreau_jean_systems::compute_total_dofs)
struct dofs : kernel<involved<systems>, rank, in<"mass_matrix">>, assembles {
  static void apply(auto& p, auto index, const auto& mass_matrix) {
    p.sum_dofs[index + 1] = algebra::ncols(mass_matrix);
  }
};

// the interactions

// the velocities of the next step and the multipliers start from zero
struct zero_velocity : kernel<all<systems>, out<"velocity", next>> {
  static void apply(const auto&, auto& v_next) { algebra::set_zero(v_next); }
};

struct zero_multipliers
    : kernel<all<interactions>, out<"lambda">, out<"ydot">, out<"ydot", next>> {
  static void apply(const auto&, auto& lambda, auto& ydot, auto& ydot_next) {
    algebra::set_zero(lambda);
    algebra::set_zero(ydot);
    algebra::set_zero(ydot_next);
  }
};

// the Newton impact law on the free output derivative
struct impact_law
    : kernel<all<interactions>, in<"ydot">, through<"nslaw", "e">, inout<"ydot", next>> {
  static void apply(const auto&, const auto& ydot, const auto& e, auto& ydot_next) {
    ydot_next += e * ydot;
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
struct jacobians
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
struct output : kernel<relations<interactions>, relation, in<"nds", by_topology>,
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
struct not_involved : kernel<all<systems>, out<"involved", by_topology>> {
  static void apply(const auto&, auto& involved) { involved = false; }
};

// an interaction is active when its gap at the middle of the step may be
// below the threshold
struct activation
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
struct involved_index : kernel<all<systems>, item_index, in<"involved", by_topology>,
                               out<"index", by_topology>> {
  static void apply(const auto& p, auto i, bool involved, auto& index) {
    if (involved) {
      index = p.offsets[i];
    }
  }
};

// the active interactions and the involved systems in their index sets (host
// containers, filled in the order of the items)
struct active_entry : kernel<all<interactions>, item_index, in<"activation", by_topology>>,
                      host_only,
                      ordered {
  static void apply(auto& p, auto i, bool activation) {
    if (activation) {
      p.emplace(i, p.offsets[i]);
    }
  }
};

struct involved_entry
    : kernel<all<systems>, item_index, in<"involved", by_topology>, in<"index", by_topology>>,
      host_only,
      ordered {
  static void apply(auto& p, auto i, bool involved, auto index) {
    if (involved) {
      p.emplace(i, index);
    }
  }
};

// the assembled nonsmooth problem

// the free output derivative of an active interaction in q of the problem
struct free_output : kernel<active<interactions>, rank, in<"ydot", next>>, assembles {
  static void apply(auto& p, auto row, const auto& ydot_next) {
    set_value(p.vector, row, ydot_next);
  }
};

// the blocks of H of an active interaction between rigid systems: its
// jacobians without the columns of the blocked velocity components, at the
// rows of the interaction and the columns of its systems
struct h_blocks
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
struct fem_h_blocks
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
struct fem_self_h_block : kernel<active<interactions>, rank, in<"h_matrix1">,
                                 in<"nds", by_topology>, through<"ds2", "index", by_topology>,
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

}  // namespace siconos::simul::moreau_jean
