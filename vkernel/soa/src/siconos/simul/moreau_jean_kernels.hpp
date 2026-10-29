#pragma once

#include "siconos/simul/kernels.hpp"

namespace siconos::simul::moreau_jean {

/* The Moreau-Jean time integration as kernel bricks (see storage/kernels.hpp and
   simul/kernels.hpp).

   p is the parameters of the step:
     p.h, p.theta        the time step and the theta of the theta-method,
     p.contact_velocity  (index, v): the velocity correction of the involved
                         system of rank index, from the solution of the
                         nonsmooth problem. */

using namespace kernels;

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

using position = position_on<all<system>>;
using involved_position = position_on<involved<system>>;

// free flight of a rigid system (K = 0), M^-1 f being kept by the element
struct free_velocity : kernel<all<system>, in<"velocity">, in<"minv_f", element>,
                              in<"minv_f", element, next>, out<"velocity", next>> {
  static void apply(const auto& p, const auto& v, const auto& minv_f, const auto& minv_f_next,
                    auto& v_next) {
    v_next = v + p.h * p.theta * minv_f_next + p.h * (1 - p.theta) * minv_f;
  }
};

// the contact impulse of an involved system
struct contact_velocity : kernel<involved<system>, rank, inout<"velocity", next>> {
  static void apply(const auto& p, auto index, auto& v_next) {
    v_next += p.contact_velocity(index, v_next);
  }
};

// the boundary conditions of a system: its blocked velocity components are
// zero, and its position follows. Only the systems with boundary conditions
// are visited when bc_velocities_0 is a sparse set.
struct boundary_conditions
    : kernel<entries<system, "bc_velocities_0", kernels::topology>, entry, in<"q">,
             in<"velocity">, inout<"velocity", next>, out<"q", next>> {
  static void apply(const auto& p, const auto& blocked, const auto& q, const auto& v,
                    auto& v_next, auto& q_next) {
    if (blocked.size() == 0) return;
    for (auto k : blocked) {
      v_next(k) = 0.0;
    }
    position::apply(p, q, v, v_next, q_next);
  }
};

}  // namespace siconos::simul::moreau_jean
