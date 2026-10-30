#pragma once

#include <SolverOptions.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <format>
#include <stdexcept>
#include <vector>

#include "siconos/simul/contacts.hpp"
#include "siconos/simul/kinds.hpp"
#include "siconos/simul/moreau_jean.hpp"
#include "siconos/simul/theta_method.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::anitescu_tasora {

/* The Anitescu-Tasora scheme (Anitescu & Tasora 2010, Tasora & Anitescu
   2011), here for the comparisons with Moreau-Jean: the Coulomb law is
   replaced by its convexification, the cone complementarity problem (CCP)

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
   diagonal M).

   The second variant solves the same CCP with APGD, the accelerated
   projected gradient descent of Mazhar, Heyn, Negrut & Tasora (ACM TOG
   2015): the CCP is the optimality condition of the convex problem

     min 1/2 lambda^T W lambda + q^T lambda,  lambda in K

   on the assembled problem of Moreau-Jean (W = H A^-1 H^T, q = H v_free +
   e ydot(k)): the same assembly, the same velocity correction, only the
   solve differs (rigid and FEM systems). Nesterov's extrapolation, the step
   1 / L by backtracking on the Lipschitz estimate L, the adaptive restart;
   warm started from lambda at the previous step; it stops when the natural
   map residual ||lambda - P_K(lambda - (W lambda + q))|| / (1 + ||q||) is
   below the tolerance of the solver options. */

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

// the solution of the CCP (W, q, mu) by APGD: lambda, from its initial
// value; n rows, d = 2 (normal, tangent) per contact, or 1 without friction
struct apgd_result {
  std::size_t iterations = 0;
  double residual = 0;
};

inline apgd_result apgd(NumericsMatrix* w, const double* q, const double* mu, double* lambda,
                        std::size_t n, std::size_t d, double tolerance,
                        std::size_t max_iterations) {
  using vector = std::vector<double>;
  auto dot = [n](const vector& a, const vector& b) {
    double s = 0;
    for (std::size_t k = 0; k < n; ++k) s += a[k] * b[k];
    return s;
  };
  auto project = [&](vector& x) {
    for (std::size_t c = 0; c < n / d; ++c) {
      if (d == 2) {
        Eigen::Map<Eigen::Vector2d> contact(&x[2 * c]);
        project_on_cone(contact, mu[c]);
      } else {
        x[c] = std::max(x[c], 0.);
      }
    }
  };
  auto product = [&](const vector& x, vector& wx) {  // wx = W x
    NM_gemv(1., w, x.data(), 0., wx.data());
  };
  const vector qv(q, q + n);
  const double norm_q = std::sqrt(dot(qv, qv));
  // ||x - P_K(x - (W x + q))|| / (1 + ||q||)
  auto residual = [&](const vector& x, const vector& wx) {
    vector z(n);
    for (std::size_t k = 0; k < n; ++k) z[k] = x[k] - (wx[k] + q[k]);
    project(z);
    for (std::size_t k = 0; k < n; ++k) z[k] -= x[k];
    return std::sqrt(dot(z, z)) / (1 + norm_q);
  };

  vector gamma(lambda, lambda + n), y = gamma, best = gamma, gamma_new(n), wy(n), wg(n), g(n);

  // the first Lipschitz estimate, between gamma and gamma + 1
  vector e(n, 1.), we(n);
  product(e, we);
  double lipschitz = std::max(std::sqrt(dot(we, we) / n), 1e-300);
  double theta = 1;
  apgd_result result{0, INFINITY};

  product(gamma, wg);
  result.residual = residual(gamma, wg);
  while (result.residual > tolerance && result.iterations < max_iterations) {
    ++result.iterations;
    product(y, wy);
    for (std::size_t k = 0; k < n; ++k) g[k] = wy[k] + q[k];

    // the projected gradient step, the step 1 / L shortened until the
    // quadratic upper bound f(x) <= f(y) + g.(x - y) + L/2 |x - y|^2 holds:
    // for a quadratic f, d^T W d <= L d^T d with d = x - y, computed as
    // d.(W x - W y) without the cancellation of the differences of f
    for (;;) {
      for (std::size_t k = 0; k < n; ++k) gamma_new[k] = y[k] - g[k] / lipschitz;
      project(gamma_new);
      product(gamma_new, wg);
      vector dy(n), wdy(n);
      for (std::size_t k = 0; k < n; ++k) {
        dy[k] = gamma_new[k] - y[k];
        wdy[k] = wg[k] - wy[k];
      }
      if (dot(dy, wdy) <= lipschitz * dot(dy, dy)) break;
      lipschitz *= 2;
    }

    const double r = residual(gamma_new, wg);
    if (r < result.residual) {
      result.residual = r;
      best = gamma_new;
    }

    // Nesterov's extrapolation, restarted when the step goes uphill
    const double theta_new = (-theta * theta + theta * std::sqrt(theta * theta + 4)) / 2;
    const double beta = theta * (1 - theta) / (theta * theta + theta_new);
    vector step(n);
    for (std::size_t k = 0; k < n; ++k) step[k] = gamma_new[k] - gamma[k];
    if (dot(g, step) > 0) {
      y = gamma_new;
      theta = 1;
    } else {
      for (std::size_t k = 0; k < n; ++k) y[k] = gamma_new[k] + beta * step[k];
      theta = theta_new;
    }
    gamma = gamma_new;
    lipschitz *= 0.9;
  }
  std::copy(best.begin(), best.end(), lambda);
  return result;
}

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

struct inactive_multipliers : stage<"lambda of the interactions that are not active is zero"> {
  static void run(auto& c) { c.kind.template run<inactive_multiplier>(c.step, nullptr); }
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

// the tolerance and the maximum number of sweeps (iterations) of the solver
// options
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
struct check : stage<"the number of iterations, a warning if lambda did not converge">,
               host_code {
  static void run(auto& c) {
    auto osnsp = c.time_stepping.one_step_nonsmooth_problem();
    const bool ok = c.change <= tolerance(c);
    storage::attr<"number_of_iterations">(osnsp) = c.iterations;
    storage::attr<"converged">(osnsp) = ok;
    storage::attr<"solver_error">(osnsp) = c.change;
    if (!ok) {
      std::fputs(std::format("warning: the Anitescu-Tasora solver did not converge: {} "
                             "iterations, error {:.3e}, tolerance {:.1e}\n",
                             c.iterations, c.change, tolerance(c))
                     .c_str(),
                 stderr);
    }
  }
};

// the CCP on the assembled problem, by APGD
struct apgd_solve : stage<"the cone complementarity problem (W, q, mu) by APGD: lambda">,
                    host_code {
  static void run(auto& c) {
    auto osi = c.integrator;
    auto& lambda = osi.lambda_vector_assembled();
    const std::size_t n = lambda._v->size0;
    const double* mu = nullptr;
    std::size_t d = 1;
    if constexpr (osi.with_friction()) {
      mu = osi.mu_vector_assembled()._v->dense;
      d = 2;
    }
    const auto result =
        apgd(osi.w_matrix_assembled()._m, osi.q_nsp_vector_assembled()._v->dense, mu,
             lambda._v->dense, n, d, tolerance(c), max_iterations(c));
    c.iterations = result.iterations;
    c.change = result.residual;
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

namespace mj = moreau_jean::stages;

// the same CCP by APGD, on the assembled problem of Moreau-Jean
using apgd_time_step_schedule = sequence<
    "Anitescu-Tasora time step, APGD",
    for_systems<th::iteration_matrix, th::free_state>,
    for_interactions<co::relation_sets, co::jacobians, co::outputs, co::free_outputs>,
    co::active_set,
    when<co::with_contacts,
         co::assembly_setup,
         for_interactions<co::time_invariant_jacobians, at::inactive_multipliers, co::h_matrix,
                          co::multipliers, mj::impact_law, mj::free_velocity_output>,
         for_systems<co::mass_matrix, co::stiffness_matrix>,
         mj::delassus_operator,
         at::apgd_solve,
         at::check,
         mj::velocity_correction,
         for_interactions<co::keep_multipliers>,
         for_systems<mj::contact_update>>>;

// the state before the first step: as moreau_jean
using initial_schedule = moreau_jean::initial_schedule;

// clang-format on

struct scheme {
  using time_step_schedule = anitescu_tasora::time_step_schedule;
  using initial_schedule = anitescu_tasora::initial_schedule;
  template <typename System>
  using system_storage = theta_method::system_storage<System>;
  template <typename Interaction>
  using interaction_storage = contacts::interaction_storage<Interaction>;
};

struct apgd_scheme : scheme {
  using time_step_schedule = anitescu_tasora::apgd_time_step_schedule;
};

}  // namespace siconos::simul::anitescu_tasora
