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
     u = H v(k + 1) + e ydot_N(k)
     (u, lambda) in the Coulomb complementarity, friction coefficients mu
     q(k + 1) = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   numerics solves it as (A, H^T, q = A v_free, b = e ydot_N(k), mu), and
   gives v(k + 1) and lambda: no A^-1 H^T, the FEM systems in contact stay
   sparse. With a linear A it is the same problem as the local one, the
   trajectories are those of moreau_jean.hpp up to the solver tolerance. */

using namespace kernels;
using moreau_jean::block_vector;

// the global problem of the step, kept by the integrator (numerics)
struct global_problem : nonsmooth_problem<GlobalFrictionContactProblem> {};
using theta_method::time_step;
using theta_method::update_position_involved;

// b of an active interaction: its impact term e ydot_N(step), no
// tangential restitution
struct assemble_impact_term
    : kernel<active<interactions>, rank, in<"ydot">, through<"nslaw", "e">>,
      assembles {
  static void apply(auto& p, auto row, const auto& ydot, const auto& e) {
    std::decay_t<decltype(ydot)> b = std::decay_t<decltype(ydot)>::Zero();
    b(0) = e * ydot(0);
    set_value(p.vector, row, b);
  }
};

// an assembled vector of the systems, the system of a rank at at(rank)
template <typename Vector, typename At>
struct system_vector {
  Vector vector;
  At at;
};

// the free velocity of an involved system in the assembled velocities
struct assemble_free_velocity : kernel<involved<systems>, rank, in<"velocity", next>>,
                                assembles {
  static void apply(auto& p, auto index, const auto& v_free) {
    set_value(p.vector, p.at(index), v_free);
  }
};

// the velocity of an involved system from the solution
struct set_solution_velocity : kernel<involved<systems>, rank, out<"velocity", next>> {
  static void apply(auto& p, auto index, auto& v_next) {
    v_next = p.part(p.velocities, index, v_next);
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct assemble_impact_terms : stage<"Assemble b = e ydot at step, the impact terms"> {
  static void run(auto& c) {
    c.kind.template run<moreau_jean_gosi::assemble_impact_term>(
        c.step, block_vector{c.kind.lambda_view(c.integrator.b_vector_assembled())});
  }
};

// on a kind of systems

struct assemble_free_velocities : stage<"Assemble v_free of the involved systems"> {
  static void run(auto& c) {
    c.kind.template run<moreau_jean_gosi::assemble_free_velocity>(
        c.step, system_vector{c.kind.velocity_view(c.integrator.velocity_vector_assembled()),
                              c.kind.at()});
  }
};

struct update_solved_state : stage<"Update v, q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    auto base = theta_method::stages::parameters(c.kind, c.h, c.integrator.theta());
    auto p = time_step{base.h, base.theta, c.kind.part(),
                       c.kind.velocity_view(c.integrator.velocity_vector_assembled())};

    // only the involved systems: the others keep their free state
    c.kind.template run<fused<set_solution_velocity, update_position_involved>>(c.step, p);

    // the involved systems may have broken them
    c.kind.template run<theta_method::apply_boundary_conditions>(c.step, p);
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

struct solve_global_problem
    : stage<
          "Solve the global frictional contact problem (A, H^T, A v_free, b, mu): v, lambda">,
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
    for_systems<th::recompute_minv_f, th::compute_free_state>,
    for_interactions<co::sync_relation_sets, co::compute_jacobians, co::compute_outputs, co::compute_free_outputs>,
    co::compute_active_set,
    when<co::with_contacts,
         co::setup_assembly,
         for_interactions<co::compute_time_invariant_jacobians, co::assemble_h_matrix, co::assemble_multipliers,
                          gosi::assemble_impact_terms>,
         for_systems<co::assemble_mass_matrix, co::assemble_stiffness_matrix, gosi::assemble_free_velocities>,
         gosi::solve_global_problem,
         for_interactions<co::keep_multipliers>,
         for_systems<gosi::update_solved_state>>>;

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
