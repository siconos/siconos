#pragma once

#include "siconos/storage/mp/mp.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul {

/* The time step of a simulation as a schedule (see storage/schedule.hpp).

   The context of a step: the time stepping, its integrator, the step, the
   time step h, and the number of involved systems once the active set is
   known. per_element runs its stages for each element of the integrator (one
   element per kind of systems and interactions), with ctx.elem the
   element. */

template <typename TimeStepping, typename Integrator, typename Step, typename H>
struct step_context {
  TimeStepping time_stepping;
  Integrator integrator;
  Step step;
  H h;
  Step involved = 0;
};

template <typename Ctx, typename Element>
struct element_context : Ctx {
  Element elem;
};

// stages run for each element, each element after the other
template <storage::pattern::string_literal Name, typename... Stages>
struct per_element {
  static constexpr std::string_view name{Name.value};
  using stages = std::tuple<Stages...>;
  static void run(auto& ctx, auto&& run_stages) {
    storage::mp::for_each(ctx.integrator.elements(), [&](auto elem) {
      element_context<std::decay_t<decltype(ctx)>, decltype(elem)> c{ctx, elem};
      run_stages(c);
      ctx.involved = c.involved;
    });
  }
};

namespace moreau_jean::stages {

using storage::schedule::condition;
using storage::schedule::stage;

// the stages of an element

struct relation_sets : stage<"the relation sets hold every interaction"> {
  static void run(auto& c) { c.elem.sync_relation_sets(c.step); }
};

struct iteration_matrix : stage<"M^-1 f, if f or K change"> {
  static void run(auto& c) { c.elem.update_iteration_matrix(c.step); }
};

struct jacobians : stage<"the jacobians H, if the relations change"> {
  static void run(auto& c) { c.elem.update_h_matrices(c.step); }
};

struct free_state : stage<"the free state: v, q at step + 1 without the contact impulses"> {
  static void run(auto& c) { c.elem.compute_free_state(c.step, c.h); }
};

struct outputs : stage<"the outputs y, ydot at step"> {
  static void run(auto& c) { c.elem.compute_output(c.step); }
};

struct free_outputs : stage<"the free outputs y, ydot at step + 1"> {
  static void run(auto& c) { c.elem.compute_output(c.step + 1); }
};

struct time_invariant_jacobians : stage<"the jacobians H, if the relations do not change"> {
  static void run(auto& c) { c.elem.ensure_h_matrices(c.step); }
};

struct mass_matrix : stage<"M of the involved systems"> {
  static void run(auto& c) { c.elem.assemble_mass_matrix_for_involved_ds(c.step); }
};

struct stiffness_matrix : stage<"K of the involved FEM systems"> {
  static void run(auto& c) { c.elem.assemble_k_matrix_for_involved_ds(c.step); }
};

struct multipliers : stage<"lambda, ydot, mu of the active interactions"> {
  static void run(auto& c) { c.elem.assemble_vectors(c.step); }
};

struct impact_law : stage<"the Newton impact law on the free ydot at step + 1"> {
  static void run(auto& c) { c.elem.nsl_effect_on_free_output(c.step); }
};

struct keep_multipliers : stage<"lambda of the active interactions from the solution"> {
  static void run(auto& c) { c.elem.keep_lambdas(c.step); }
};

struct contact_update : stage<"v, q at step + 1 of the involved systems"> {
  static void run(auto& c) { c.elem.update_velocities_and_positions(c.step, c.h); }
};

// the stages of the integrator

struct active_set
    : stage<"the active set: the interactions whose gap may close, the systems they involve"> {
  static void run(auto& c) {
    c.involved = c.integrator.compute_active_interactions(c.step, c.h);
  }
};

struct contacts : condition<"there are contacts"> {
  static bool holds(const auto& c) { return c.involved > 0; }
};

struct assembly_setup : stage<"the sizes and offsets of the assembled problem"> {
  static void run(auto& c) { c.integrator.assemble_setup(c.step); }
};

struct h_matrix : stage<"H of the active interactions"> {
  static void run(auto& c) { c.integrator.assemble_h_matrix_for_involved_ds(c.step); }
};

struct free_velocity_output : stage<"q = H v_free, the free output derivatives"> {
  static void run(auto& c) { c.integrator.compute_q_nsp_vector_assembled(c.step); }
};

struct delassus : stage<"W = H M^-1 H^T, or H (M + h^2 theta^2 K)^-1 H^T"> {
  static void run(auto& c) { c.integrator.compute_w_matrix(c.step, c.h); }
};

struct solve : stage<"the frictional contact problem (W, q, mu): lambda"> {
  static void run(auto& c) {
    using time_stepping_t = typename std::decay_t<decltype(c.time_stepping)>::type;
    c.time_stepping.template solve_nonsmooth_problem<typename time_stepping_t::formulation_t>(
        c.step);
  }
};

struct impulse : stage<"the velocity correction M^-1 H^T lambda"> {
  static void run(auto& c) { c.integrator.compute_input(c.step, c.h); }
};

}  // namespace moreau_jean::stages

namespace moreau_jean {

using storage::schedule::sequence;
using storage::schedule::when;

// the Moreau-Jean time step
using time_step_schedule = sequence<
    "Moreau-Jean time step",
    per_element<"for each element", stages::relation_sets, stages::iteration_matrix,
                stages::jacobians, stages::free_state, stages::outputs, stages::free_outputs>,
    stages::active_set,
    when<stages::contacts, stages::assembly_setup,
         per_element<"for each element", stages::time_invariant_jacobians>, stages::h_matrix,
         per_element<"for each element", stages::mass_matrix, stages::stiffness_matrix,
                     stages::multipliers, stages::impact_law>,
         stages::free_velocity_output, stages::delassus, stages::solve, stages::impulse,
         per_element<"for each element", stages::keep_multipliers, stages::contact_update>>>;

}  // namespace moreau_jean

}  // namespace siconos::simul
