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
   below, given by moreau_jean_element when it runs the kernel. */

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
struct minv_f : kernel<all<systems>, in<"mass_matrix">, in<"fext">, out<"minv_f", by_element>,
                       out<"minv_f", by_element, next>> {
  static void apply(const auto&, auto& mass_matrix, const auto& f, auto& minv_f,
                    auto& minv_f_next) {
    minv_f = f;
    algebra::solve_in_place(mass_matrix, minv_f);
    minv_f_next = minv_f;
  }
};

// free flight of a rigid system (K = 0)
struct free_velocity : kernel<all<systems>, in<"velocity">, in<"minv_f", by_element>,
                              in<"minv_f", by_element, next>, out<"velocity", next>> {
  static void apply(const auto& p, const auto& v, const auto& minv_f, const auto& minv_f_next,
                    auto& v_next) {
    v_next = v + p.h * p.theta * minv_f_next + p.h * (1 - p.theta) * minv_f;
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
struct diagonal_block : kernel<involved<systems>, rank, in<Name>> {
  static void apply(auto& p, auto index, const auto& block) {
    set_value(p.matrix, index, index, block);
  }
};

// the number of dofs of an involved system, before the scan that gives the
// offsets (see moreau_jean_element::compute_total_dofs)
struct dofs : kernel<involved<systems>, rank, in<"mass_matrix">> {
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
struct assemble_multipliers : kernel<active<interactions>, rank, in<"lambda">, in<"ydot">> {
  static void apply(auto& p, auto index, const auto& lambda, const auto& ydot) {
    set_value(p.lambdas, index, lambda);
    set_value(p.ydots, index, ydot);
  }
};

struct assemble_friction : kernel<active<interactions>, rank, through<"nslaw", "mu">> {
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
      reads_world {
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
                reads_world {
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

}  // namespace siconos::simul::moreau_jean
