#pragma once

#include <SolverOptions.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <stdexcept>

#include "siconos/simul/contacts.hpp"
#include "siconos/simul/kinds.hpp"
#include "siconos/simul/moreau_jean.hpp"
#include "siconos/simul/theta_method.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::anitescu_tasora {

/* The Anitescu-Tasora scheme (Anitescu & Tasora 2010, Tasora & Anitescu
   2011), a reference for the comparisons with Moreau-Jean, not a scheme for
   production: the Coulomb law is replaced by its convexification, the cone
   complementarity problem (CCP)

     lambda in K,  u in K*,  u . lambda = 0
     K = {|lambda_T| <= mu lambda_N},  K* = {u_N >= mu |u_T|}

   whose solutions slide with a normal velocity u_N = mu |u_T| > 0 instead
   of u_N = 0: the artificial dilatancy, the spurious separation of a
   sliding body. With mu = 0 both laws are the same.

   The CCP is solved without assembly, interaction by interaction (projected
   Gauss-Seidel, matrix free): for each active interaction, with v the
   velocities of its systems so far,

     u      = H v + e ydot(k)
     lambda <- P_K(lambda - eta u),  eta = 2 / trace(H M^-1 H^T)
     v      <- v + M^-1 H^T (the change of lambda)

   until the largest change of lambda in a sweep is below the tolerance of
   the solver options (dparam TOL, relative to the largest lambda), at most
   iparam MAX_ITER sweeps. Warm start: lambda starts from its value at the
   previous step (zero for an interaction that was not active), v from
   v_free + M^-1 H^T lambda. Nothing is assembled but the contact impulses H^T lambda of the
   involved systems, for the outputs.

   Everything else is the Moreau-Jean step (moreau_jean.hpp): the theta
   method, the active set, the Newton impact law, q(k + 1) from v(k + 1).
   The stabilization of Anitescu-Tasora (the gap / h term) is not used, so
   that a comparison isolates the convexified cone. Rigid systems only (a
   diagonal M). */

using namespace kernels;

// the closed form projection on the 2D cone K of friction coefficient mu
void project_on_cone(auto& lambda, auto mu) {
  const auto n = lambda(0);
  const auto t = std::abs(lambda(1));
  if (t <= mu * n) {
    return;  // in K
  }
  if (mu * t <= -n) {
    lambda(0) = 0;  // in the polar cone: the vertex
    lambda(1) = 0;
    return;
  }
  // on the boundary ray of the side of lambda_T
  const auto n_proj = (n + mu * t) / (1 + mu * mu);
  lambda(1) = std::copysign(mu * n_proj, lambda(1));
  lambda(0) = n_proj;
}

// the parameters of a sweep: its largest change of lambda and largest lambda
template <typename S>
struct sweep_state {
  S change = 0;
  S largest = 0;
};

// the multiplier of an interaction that is not active is zero: it starts
// from zero when it is active again
struct inactive_multiplier
    : kernel<all<interactions>, in<"activation", by_topology>, inout<"lambda">> {
  static void apply(const auto&, bool activation, auto& lambda) {
    if (!activation) {
      algebra::set_zero(lambda);
    }
  }
};

// M^-1 H^T without the blocked velocity components: the velocity change of a
// system for the multiplier of one of its interactions
auto velocity_change(auto& m, auto h, const auto& blocked, const auto& lambda, const auto& v) {
  for (auto k : blocked) h.col(k).setZero();
  std::decay_t<decltype(v)> dv = h.transpose() * lambda;
  algebra::solve_in_place(m, dv);
  return dv;
}

// the warm start: the systems of an active interaction get the impulse of
// its multiplier at the previous step
struct warm_start
    : kernel<active<interactions>, in<"nds", by_topology>, in<"h_matrix1">, in<"h_matrix2">,
             in<"lambda">, through<"ds1", "mass_matrix", by_topology>,
             through<"ds2", "mass_matrix", by_topology>,
             through<"ds1", "bc_velocities_0", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>,
             add_through<"ds1", "velocity", by_topology, next>,
             add_through<"ds2", "velocity", by_topology, next>> {
  static void apply(const auto&, auto nds, const auto& h1, const auto& h2, const auto& lambda,
                    auto& m1, auto& m2, const auto& blocked1, const auto& blocked2, auto& v1,
                    auto& v2) {
    v1 += velocity_change(m1, h1, blocked1, lambda, v1);
    if (nds == 2) {
      v2 += velocity_change(m2, h2, blocked2, lambda, v2);
    }
  }
};

// one projected Gauss-Seidel step on an active interaction between rigid
// systems: its multiplier, and the velocities of its systems by accumulation
// (add_through). A one system interaction (nds == 1) has ds1 == ds2 and its
// jacobian h_matrix1 (as contacts::h_blocks).
struct projected_gauss_seidel
    : kernel<active<interactions>, in<"nds", by_topology>, in<"h_matrix1">, in<"h_matrix2">,
             in<"ydot">, through<"nslaw", "e">, through<"nslaw", "mu">,
             through<"ds1", "mass_matrix", by_topology>,
             through<"ds2", "mass_matrix", by_topology>,
             through<"ds1", "bc_velocities_0", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>, inout<"lambda">,
             add_through<"ds1", "velocity", by_topology, next>,
             add_through<"ds2", "velocity", by_topology, next>>,
      ordered {
  static void apply(auto& p, auto nds, const auto& h1, const auto& h2, const auto& ydot,
                    const auto& e, const auto& mu, auto& m1, auto& m2, const auto& blocked1,
                    const auto& blocked2, auto& lambda, auto& v1, auto& v2) {
    using vector = std::decay_t<decltype(v1)>;
    const bool two = nds == 2;

    // the jacobians without the columns of the blocked velocity components
    auto g1 = h1;
    auto g2 = h2;
    for (auto k : blocked1) g1.col(k).setZero();
    for (auto k : blocked2) g2.col(k).setZero();

    // M^-1 G^T, row by row of G
    auto minv_gt = [](auto& m, const auto& g, auto r) {
      vector w = g.row(r).transpose();
      algebra::solve_in_place(m, w);
      return w;
    };
    const vector a1n = minv_gt(m1, g1, 0), a1t = minv_gt(m1, g1, 1);
    const vector a2n = two ? minv_gt(m2, g2, 0) : vector::Zero();
    const vector a2t = two ? minv_gt(m2, g2, 1) : vector::Zero();

    // trace(H M^-1 H^T)
    auto trace = g1.row(0).dot(a1n) + g1.row(1).dot(a1t);
    if (two) trace += g2.row(0).dot(a2n) + g2.row(1).dot(a2t);
    if (trace == 0) {
      return;  // every velocity component it acts on is blocked
    }

    // u = H v + e ydot(k)
    auto u = (g1 * v1).eval();
    if (two) u += g2 * v2;
    u += e * ydot;

    auto previous = lambda;
    lambda -= (2 / trace) * u;
    project_on_cone(lambda, mu);
    const auto d = (lambda - previous).eval();

    v1 += d(0) * a1n + d(1) * a1t;
    if (two) v2 += d(0) * a2n + d(1) * a2t;

    p.change = std::max(p.change, d.template lpNorm<Eigen::Infinity>());
    p.largest = std::max(p.largest, lambda.template lpNorm<Eigen::Infinity>());
  }
};

// the contact impulse H^T lambda of an active interaction between rigid
// systems, added to the assembled impulses of its systems (two interactions
// may add to the same system: in order). A one system interaction: ds1.
struct contact_impulse
    : kernel<active<interactions>, in<"nds", by_topology>, in<"h_matrix1">, in<"h_matrix2">,
             in<"lambda">, through<"ds1", "index", by_topology>,
             through<"ds2", "index", by_topology>,
             through<"ds1", "bc_velocities_0", by_topology>,
             through<"ds2", "bc_velocities_0", by_topology>>,
      assembles,
      ordered {
  static void apply(auto& p, auto nds, const auto& h1, const auto& h2, const auto& lambda,
                    auto j1, auto j2, const auto& blocked1, const auto& blocked2) {
    auto g1 = h1;
    auto g2 = h2;
    for (auto k : blocked1) g1.col(k).setZero();
    for (auto k : blocked2) g2.col(k).setZero();
    auto add = [&p](auto j, const auto& impulse) {
      set_value(p.vector, j, get_vector(p.vector, j) + impulse);
    };
    if (nds == 2) {
      add(j1, (g1.transpose() * lambda).eval());
      add(j2, (g2.transpose() * lambda).eval());
    } else {
      add(j1, (g1.transpose() * lambda).eval());
    }
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::condition;
using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct warm_start : stage<
                        "lambda at the previous step, zero if not active; "
                        "v <- v + M^-1 H^T lambda, the warm start"> {
  static void run(auto& c) {
    c.kind.template run<inactive_multiplier>(c.step, nullptr);
    if constexpr (std::decay_t<decltype(c.kind)>::between_rigid_systems()) {
      c.kind.template run<anitescu_tasora::warm_start>(c.step, nullptr);
    }
  }
};

struct sweep : stage<
                   "lambda <- P_K(lambda - eta (H v + e ydot)), v <- v + M^-1 H^T dlambda, "
                   "interaction by interaction"> {
  static void run(auto& c) {
    if constexpr (std::decay_t<decltype(c.kind)>::between_rigid_systems()) {
      using scalar = typename decltype(c.kind.env())::scalar;
      sweep_state<scalar> p;
      c.kind.template run<projected_gauss_seidel>(c.step, p);
      // the largest over the kinds, relative to the largest lambda
      c.change = std::max(c.change, p.change / std::max(p.largest, scalar(1e-300)));
    } else if (c.kind.number_of_interactions() > 0) {
      throw std::runtime_error("anitescu_tasora: contacts with FEM systems are not supported");
    }
  }
};

struct contact_impulses : stage<"p0 = H^T lambda, the contact impulses (outputs)"> {
  static void run(auto& c) {
    if constexpr (std::decay_t<decltype(c.kind)>::between_rigid_systems()) {
      auto rigid_systems = std::get<0>(c.integrator.systems());
      c.kind.template run<contact_impulse>(
          c.step, moreau_jean::block_vector{
                      rigid_systems.p0_view(c.integrator.p0_vector_assembled())});
    }
  }
};

// on a kind of systems

struct position_update : stage<"q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    auto p = theta_method::stages::parameters(c.kind, c.h, c.integrator.theta());
    // their velocities are already updated, by the sweeps
    c.kind.template run<theta_method::involved_position>(c.step, p);
    c.kind.template run<theta_method::boundary_conditions>(c.step, p);
  }
};

// on the integrator

// the offsets of the involved systems in the assembled vectors (no matrix)
struct impulses_setup : stage<"the offsets of the involved systems, p0 is zero">, host_code {
  static void run(auto& c) {
    auto osi = c.integrator;
    std::size_t offset = 0;
    storage::mp::for_each(osi.systems(), [&](auto kind) {
      kind.compute_total_dofs(c.step);
      kind.ds_offset() = offset;
      offset += kind.total_dofs();
    });
    algebra::resize(osi.p0_vector_assembled(), offset);  // zeros
  }
};

struct start : stage<"no sweep yet">, host_code {
  static void run(auto& c) { c.iterations = 0; }
};

struct next_sweep : stage<"a new sweep">, host_code {
  static void run(auto& c) {
    ++c.iterations;
    c.change = 0;
  }
};

// the tolerance and the maximum number of sweeps of the solver options
auto tolerance(auto& c) {
  auto options = c.time_stepping.one_step_nonsmooth_problem().options();
  return options.dparam(SICONOS_DPARAM_TOL);
}
auto max_iterations(auto& c) {
  auto options = c.time_stepping.one_step_nonsmooth_problem().options();
  return static_cast<std::size_t>(options.iparam(SICONOS_IPARAM_MAX_ITER));
}

struct converged : condition<"lambda changes less than the tolerance, or too many sweeps"> {
  static bool holds(auto& c) {
    return c.change <= tolerance(c) || c.iterations >= max_iterations(c);
  }
};

// as one_step_nonsmooth_problem::check: the result of the solver
struct check : stage<"the number of sweeps, a warning if lambda did not converge">, host_code {
  static void run(auto& c) {
    auto osnsp = c.time_stepping.one_step_nonsmooth_problem();
    const bool ok = c.change <= tolerance(c);
    storage::attr<"number_of_iterations">(osnsp) = c.iterations;
    storage::attr<"converged">(osnsp) = ok;
    storage::attr<"solver_error">(osnsp) = c.change;
    if (!ok) {
      std::fputs(std::format("warning: the Anitescu-Tasora sweeps did not converge: {} "
                             "sweeps, change {:.3e}, tolerance {:.1e}\n",
                             c.iterations, c.change, tolerance(c))
                     .c_str(),
                 stderr);
    }
  }
};

}  // namespace stages

// the schedules

using storage::schedule::repeat;
using storage::schedule::sequence;
using storage::schedule::when;

namespace th = theta_method::stages;
namespace co = contacts::stages;
namespace at = stages;

// one stage per line: the schedules read as the algorithm
// clang-format off

// the Anitescu-Tasora time step
using time_step_schedule = sequence<
    "Anitescu-Tasora time step",
    for_systems<th::iteration_matrix, th::free_state>,
    for_interactions<co::relation_sets, co::jacobians, co::outputs, co::free_outputs>,
    co::active_set,
    when<co::with_contacts,
         for_interactions<co::time_invariant_jacobians, at::warm_start>,
         at::start,
         repeat<at::converged,
                at::next_sweep,
                for_interactions<at::sweep>>,
         at::check,
         for_systems<at::position_update>,
         at::impulses_setup,
         for_interactions<at::contact_impulses>>>;

// the state before the first step: as moreau_jean
using initial_schedule = moreau_jean::initial_schedule;

// clang-format on

struct scheme {
  using time_step_schedule = anitescu_tasora::time_step_schedule;
  using initial_schedule = anitescu_tasora::initial_schedule;
};

}  // namespace siconos::simul::anitescu_tasora
