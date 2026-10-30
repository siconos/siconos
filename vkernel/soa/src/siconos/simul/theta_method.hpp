#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/algebra/numerics.hpp"
#include "siconos/simul/kernels.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::theta_method {

/* The theta-method on the smooth dynamics M dv = (f - K q) dt, the part of
   the Moreau-Jean schemes (local and global) without contacts. With A = M +
   h^2 theta^2 K the iteration matrix and F = f - K q:

     v_free(k + 1) = v(k) + A^-1 h (theta F(k + 1) + (1 - theta) F(k))
     q(k + 1)      = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   For a rigid system (K = 0, A = M) with f constant over the step, M^-1 f is
   kept on the system (minv_f): v_free = v + h theta M^-1 f(k + 1) + h (1 -
   theta) M^-1 f(k). The blocked velocity components (boundary conditions)
   are zero. */

using namespace kernels;

// the storage of the theta-method on the systems of a kind (the properties of
// the kind, see kinds.hpp): M^-1 f, declared by the kind, at one time level
// (f is constant over the step: minv_f at step and at step + 1 are the same
// value), q and the velocity at two levels, step and step + 1
template <typename System>
using system_storage =
    gather<storage::attached<System, symbol<"minv_f">, attr_t<System, "fext">>,
           storage::time_levels<attr_t<System, "q">, 2>,
           storage::time_levels<attr_t<System, "velocity">, 2>>;

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

// the velocities of the next step start from zero
struct zero_velocity : kernel<all<systems>, out<"velocity", next>> {
  static void apply(const auto&, auto& v_next) { algebra::set_zero(v_next); }
};

// the stages (see storage/schedule.hpp), on a kind of systems (see kinds.hpp)

namespace stages {

using storage::schedule::stage;

// the parameters of the kernels: the time step h, the theta of the method
auto parameters(auto kind, auto h, auto theta) {
  using scalar = typename decltype(kind.env())::scalar;
  return time_step<decltype(h), scalar>{h, theta};
}

// M^-1 f, f being constant over the step: M^-1 (f - K q) for the FEM
// systems, with the factorization of M kept on them
void compute_minv_f(auto kind, auto step) {
  if constexpr (kind.with_k_matrix()) {
    using scalar = typename decltype(kind.env())::scalar;
    auto factor = [&](auto i) -> auto& { return kind.fem_m_factor(step, i); };
    kind.template run<fem_minv_f<algebra::unbounded_vector<scalar>>>(step,
                                                                     fem_factors{factor});
  } else {
    kind.template run<minv_f>(step, nullptr);
  }
}

struct initial_state : stage<"v at step + 1 is zero, M^-1 f"> {
  static void run(auto& c) {
    c.kind.template run<zero_velocity>(c.step, nullptr);
    compute_minv_f(c.kind, c.step);
  }
};

struct iteration_matrix : stage<"M^-1 f, if f or K change"> {
  static void run(auto& c) {
    using kind_t = std::decay_t<decltype(c.kind)>;
    if constexpr (kind_t::with_k_matrix() || !kind_t::time_invariant_fext()) {
      compute_minv_f(c.kind, c.step);
    }
  }
};

struct free_state : stage<"the free state: v, q at step + 1 without the contact impulses"> {
  static void run(auto& c) {
    auto p = parameters(c.kind, c.h, c.integrator.theta());

    if constexpr (std::decay_t<decltype(c.kind)>::with_k_matrix()) {
      // A dv = h theta F(step + 1) + h (1 - theta) F(step), A factorized
      // once for the simulation
      using env_t = decltype(c.kind.env());
      using vector = typename env_t::template unbounded_vector<typename env_t::scalar>;
      auto factor = [&](auto i) -> auto& {
        return c.kind.fem_w_factor(c.step, i, c.h, c.integrator.theta());
      };
      c.kind.template run<fused<fem_free_velocity<vector>, position>>(
          c.step, fem_time_step{p.h, p.theta, factor});
    } else {
      c.kind.template run<fused<free_velocity, position>>(c.step, p);
    }

    c.kind.template run<boundary_conditions>(c.step, p);
  }
};

}  // namespace stages

}  // namespace siconos::simul::theta_method
