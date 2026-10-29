#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul {

/* The time step of a simulation as a schedule (see storage/schedule.hpp).

   The context of a step: the time stepping, its integrator, the step, the
   time step h, and the number of involved systems once the active set is
   known. for_systems (for_interactions) runs its stages for each kind of
   systems (of interactions) of the integrator, with ctx.kind the kind (see
   moreau_jean_kinds.hpp). */

template <typename TimeStepping, typename Integrator, typename Step, typename H>
struct step_context {
  TimeStepping time_stepping;
  Integrator integrator;
  Step step;
  H h;
  Step involved = 0;
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

namespace moreau_jean::stages {

using storage::schedule::condition;
using storage::schedule::host_code;
using storage::schedule::stage;

// the stages of a kind of systems

struct iteration_matrix : stage<"M^-1 f, if f or K change"> {
  static void run(auto& c) { c.kind.update_iteration_matrix(c.step); }
};

struct free_state : stage<"the free state: v, q at step + 1 without the contact impulses"> {
  static void run(auto& c) { c.kind.compute_free_state(c.step, c.h, c.integrator.theta()); }
};

struct mass_matrix : stage<"M of the involved systems"> {
  static void run(auto& c) {
    c.kind.assemble_mass_matrix(c.step, c.integrator.mass_matrix_assembled());
  }
};

struct stiffness_matrix : stage<"K of the involved FEM systems"> {
  static void run(auto& c) {
    c.kind.assemble_k_matrix(c.step, c.integrator.k_matrix_assembled());
  }
};

struct contact_update : stage<"v, q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    c.kind.update_velocities_and_positions(c.step, c.h, c.integrator.theta(),
                                           c.integrator.velocity_vector_assembled());
  }
};

// the stages of a kind of interactions

struct relation_sets : stage<"the relation sets hold every interaction">, host_code {
  static void run(auto& c) { c.kind.sync_relation_sets(c.step); }
};

struct jacobians : stage<"the jacobians H, if the relations change"> {
  static void run(auto& c) { c.kind.update_h_matrices(c.step); }
};

struct outputs : stage<"the outputs y, ydot at step"> {
  static void run(auto& c) { c.kind.compute_output(c.step); }
};

struct free_outputs : stage<"the free outputs y, ydot at step + 1"> {
  static void run(auto& c) { c.kind.compute_output(c.step + 1); }
};

struct time_invariant_jacobians : stage<"the jacobians H, if the relations do not change"> {
  static void run(auto& c) { c.kind.ensure_h_matrices(c.step); }
};

struct h_matrix : stage<"H of the active interactions"> {
  static void run(auto& c) { c.integrator.assemble_h_matrix(c.step, c.kind); }
};

struct multipliers : stage<"lambda, ydot, mu of the active interactions"> {
  static void run(auto& c) {
    c.kind.assemble_vectors(c.step, c.integrator.lambda_vector_assembled(),
                            c.integrator.ydot_vector_assembled(),
                            c.integrator.mu_vector_assembled());
  }
};

struct impact_law : stage<"the Newton impact law on the free ydot at step + 1"> {
  static void run(auto& c) { c.kind.nsl_effect_on_free_output(c.step); }
};

struct free_velocity_output : stage<"q = H v_free, the free output derivatives"> {
  static void run(auto& c) {
    c.kind.free_output(c.step, c.integrator.q_nsp_vector_assembled());
  }
};

struct keep_multipliers : stage<"lambda of the active interactions from the solution"> {
  static void run(auto& c) {
    c.kind.keep_lambdas(c.step, c.integrator.lambda_vector_assembled());
  }
};

// the stages of the integrator

struct active_set
    : stage<"the active set: the interactions whose gap may close, the systems they involve">,
      host_code {
  static void run(auto& c) {
    c.involved = c.integrator.compute_active_interactions(c.step, c.h);
  }
};

struct contacts : condition<"there are contacts"> {
  static bool holds(const auto& c) { return c.involved > 0; }
};

struct assembly_setup : stage<"the sizes and offsets of the assembled problem">, host_code {
  static void run(auto& c) { c.integrator.assemble_setup(c.step); }
};

struct delassus : stage<"W = H M^-1 H^T, or H (M + h^2 theta^2 K)^-1 H^T">, host_code {
  static void run(auto& c) { c.integrator.compute_w_matrix(c.step, c.h); }
};

struct solve : stage<"the frictional contact problem (W, q, mu): lambda">, host_code {
  static void run(auto& c) {
    using time_stepping_t = typename std::decay_t<decltype(c.time_stepping)>::type;
    c.time_stepping.template solve_nonsmooth_problem<typename time_stepping_t::formulation_t>(
        c.step);
  }
};

struct impulse : stage<"the velocity correction M^-1 H^T lambda">, host_code {
  static void run(auto& c) { c.integrator.compute_input(c.step, c.h); }
};

}  // namespace moreau_jean::stages

namespace moreau_jean {

using storage::schedule::sequence;
using storage::schedule::when;

// the Moreau-Jean time step
using time_step_schedule =
    sequence<"Moreau-Jean time step",
             for_systems<stages::iteration_matrix, stages::free_state>,
             for_interactions<stages::relation_sets, stages::jacobians, stages::outputs,
                              stages::free_outputs>,
             stages::active_set,
             when<stages::contacts, stages::assembly_setup,
                  for_interactions<stages::time_invariant_jacobians, stages::h_matrix,
                                   stages::multipliers, stages::impact_law,
                                   stages::free_velocity_output>,
                  for_systems<stages::mass_matrix, stages::stiffness_matrix>, stages::delassus,
                  stages::solve, stages::impulse, for_interactions<stages::keep_multipliers>,
                  for_systems<stages::contact_update>>>;

}  // namespace moreau_jean

}  // namespace siconos::simul
