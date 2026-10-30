#pragma once

#include "siconos/algebra/numerics.hpp"
#include "siconos/simul/contacts.hpp"
#include "siconos/simul/kinds.hpp"
#include "siconos/simul/moreau_jean.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/theta_method.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::moreau_jean_gosi {

/* The Moreau-Jean scheme on the global problem (MoreauJeanGOSI in the
   Siconos kernel): the same theta-method and Newton impact law as
   moreau_jean.hpp, but the frictional contact problem is solved on the
   velocities of the systems, without the Delassus operator. A = M + h^2
   theta^2 K the iteration matrix, H the jacobians of the active
   interactions:

     v_free         the free state (theta_method.hpp)
     A v(k + 1) = A v_free + H^T lambda
     u = H v(k + 1) + e ydot(k)
     (u, lambda) in the Coulomb complementarity, friction coefficients mu
     q(k + 1) = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   numerics solves it as (A, H^T, q = A v_free, b = e ydot(k), mu), and
   gives v(k + 1) and lambda: no A^-1 H^T, the FEM systems in contact stay
   sparse. With a linear A it is the same problem as the local one, the
   trajectories are those of moreau_jean.hpp up to the solver tolerance. */

using namespace kernels;
using moreau_jean::block_vector;

// the global problem of the step, kept by the integrator (numerics)
struct global_problem : nonsmooth_problem<GlobalFrictionContactProblem> {};
using theta_method::involved_position;
using theta_method::time_step;

// b of an active interaction: its impact term e ydot(step)
struct impact_term : kernel<active<interactions>, rank, in<"ydot">, through<"nslaw", "e">>,
                     assembles {
  static void apply(auto& p, auto row, const auto& ydot, const auto& e) {
    set_value(p.vector, row, e * ydot);
  }
};

// an assembled vector of the systems, the system of a rank at at(rank)
template <typename Vector, typename At>
struct system_vector {
  Vector vector;
  At at;
};

// the free velocity of an involved system in the assembled velocities
struct free_velocity : kernel<involved<systems>, rank, in<"velocity", next>>, assembles {
  static void apply(auto& p, auto index, const auto& v_free) {
    set_value(p.vector, p.at(index), v_free);
  }
};

// the velocity of an involved system from the solution
struct solution_velocity : kernel<involved<systems>, rank, out<"velocity", next>> {
  static void apply(auto& p, auto index, auto& v_next) {
    v_next = p.part(p.velocities, index, v_next);
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct impact_term : stage<"b = e ydot at step, the impact terms"> {
  static void run(auto& c) {
    c.kind.template run<moreau_jean_gosi::impact_term>(
        c.step, block_vector{c.kind.lambda_view(c.integrator.b_vector_assembled())});
  }
};

// on a kind of systems

struct free_velocities : stage<"v_free of the involved systems"> {
  static void run(auto& c) {
    c.kind.template run<moreau_jean_gosi::free_velocity>(
        c.step, system_vector{c.kind.velocity_view(c.integrator.velocity_vector_assembled()),
                              c.kind.at()});
  }
};

struct velocity_update : stage<"v, q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    auto base = theta_method::stages::parameters(c.kind, c.h, c.integrator.theta());
    auto p = time_step{base.h, base.theta, c.kind.part(),
                       c.kind.velocity_view(c.integrator.velocity_vector_assembled())};

    // only the involved systems: the others keep their free state
    c.kind.template run<fused<solution_velocity, involved_position>>(c.step, p);

    // the involved systems may have broken them
    c.kind.template run<theta_method::boundary_conditions>(c.step, p);
  }
};

// on the integrator

// the global problem with the iteration matrix A: q = A v_free, then v(k + 1)
// in place of v_free, lambda
void solve(auto osi, auto osnsp, auto& a) {
  auto& v = osi.velocity_vector_assembled();
  auto& q = osi.q_nsp_vector_assembled();
  resize(q, size0(v));
  NM_gemv(1., a._m, v._v->dense, 0., q._v->dense);

  // the problem of the integrator, the solver options of the time stepping
  osnsp.solve_global(*osi.global_problem().instance(), a, osi.h_matrix_assembled(), q,
                     osi.b_vector_assembled(), osi.mu_vector_assembled(),
                     osi.lambda_vector_assembled(), osi.ydot_vector_assembled(), v);
}

struct global_solve
    : stage<"the global frictional contact problem (A, H^T, A v_free, b, mu): v, lambda">,
      host_code {
  static void run(auto& c) {
    auto osi = c.integrator;
    auto osnsp = c.time_stepping.one_step_nonsmooth_problem();
    if constexpr (osi.with_k_matrix()) {
      if (algebra::nnz(osi.k_matrix_assembled()) > 0) {
        // A = M + h^2 theta^2 K
        algebra::add(1., osi.mass_matrix_assembled(), c.h * c.h * osi.theta() * osi.theta(),
                     osi.k_matrix_assembled(), osi.iteration_matrix_assembled());
        solve(osi, osnsp, osi.iteration_matrix_assembled());
        return;
      }
    }
    // A = M: no FEM system in contact
    solve(osi, osnsp, osi.mass_matrix_assembled());
  }
};

}  // namespace stages

// the schedules

using storage::schedule::sequence;
using storage::schedule::when;

namespace th = theta_method::stages;
namespace co = contacts::stages;
namespace gosi = stages;

// one stage per line: the schedules read as the algorithm
// clang-format off

// the Moreau-Jean global time step
using time_step_schedule = sequence<
    "Moreau-Jean global time step",
    for_systems<th::iteration_matrix, th::free_state>,
    for_interactions<co::relation_sets, co::jacobians, co::outputs, co::free_outputs>,
    co::active_set,
    when<co::with_contacts,
         co::assembly_setup,
         for_interactions<co::time_invariant_jacobians, co::h_matrix, co::multipliers,
                          gosi::impact_term>,
         for_systems<co::mass_matrix, co::stiffness_matrix, gosi::free_velocities>,
         gosi::global_solve,
         for_interactions<co::keep_multipliers>,
         for_systems<gosi::velocity_update>>>;

// the state before the first step: as moreau_jean
using initial_schedule = moreau_jean::initial_schedule;

// clang-format on

struct scheme {
  using time_step_schedule = moreau_jean_gosi::time_step_schedule;
  using initial_schedule = moreau_jean_gosi::initial_schedule;
  template <typename System>
  using system_storage = theta_method::system_storage<System>;
  template <typename Interaction>
  using interaction_storage = contacts::interaction_storage<Interaction>;
};

}  // namespace siconos::simul::moreau_jean_gosi
