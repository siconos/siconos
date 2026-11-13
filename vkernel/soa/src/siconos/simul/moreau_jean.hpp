#pragma once

#include <stdexcept>

#include "siconos/algebra/numerics.hpp"
#include "siconos/simul/contacts.hpp"
#include "siconos/simul/kinds.hpp"
#include "siconos/simul/theta_method.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::moreau_jean {

/* The Moreau-Jean scheme (Jean 1999, Acary & Brogliato 2008): the
   theta-method on the smooth dynamics (see theta_method.hpp), the Newton
   impact law and Coulomb friction at the velocity level, on the local
   problem of the contacts. A = M + h^2 theta^2 K the iteration matrix, H the
   jacobians of the active interactions:

     v_free         the free state (theta_method.hpp)
     W  = H A^-1 H^T                   the Delassus operator
     q  = H v_free + e ydot_N(k)       the free output derivatives
     lambda solution of the frictional contact problem (W, q, mu)
     v(k + 1) = v_free + A^-1 H^T lambda
     q(k + 1) = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   With the Frémond impact law (model::fremond_impact_friction, Acary &
   Collins-Craft, JTCAM 2025) Coulomb friction is written on u(k + theta)
   = theta u(k + 1) + (1 - theta) u(k) instead of u(k + 1): the problem of a
   contact is w = theta W lambda + q with q = theta H v_free + (1 - theta)
   ydot(k) + (theta (1 + e) - 1) ydot_N(k). The cone complementarity is
   invariant under a positive scaling of the (w, q) of a contact, so it is
   divided by theta: W is unchanged, and so is the normal row of q, only

     q_T = H_T v_free + (1 - theta) / theta ydot_T(k)

   lambda is the same impulse, the velocity update is unchanged, and the
   Newton and the Frémond contacts can share one problem. With theta = 1 it
   is the Newton impact law.

   The kernels of the scheme, its stages, then the schedules (the time step,
   the initial state). The rest of a step is in contacts.hpp. */

using namespace kernels;
using theta_method::time_step;
using theta_method::update_position_involved;

// an assembled vector
template <typename Vector>
struct block_vector {
  Vector vector;
};

// the contact impulse of an involved system
struct update_contact_next_velocity
    : kernel<involved<systems>, rank, inout<"velocity", next>> {
  static void apply(auto& p, auto index, auto& v_next) {
    v_next += p.part(p.velocities, index, v_next);
  }
};

// the Newton impact law on the free output derivative: e on the normal
// component only, no tangential restitution (et = 0 in the kernel)
struct apply_newton_impact_law
    : kernel<all<interactions>, in<"ydot">, through<"nslaw", "e">, inout<"ydot", next>> {
  static void apply(const auto&, const auto& ydot, const auto& e, auto& ydot_next) {
    ydot_next(0) += e * ydot(0);
  }
};

// the Frémond impact law on the free output derivative, divided by theta
// (see above): the Newton impact law on the normal component, Coulomb
// friction on u(k + theta) through the tangential ones
struct apply_fremond_impact_law
    : kernel<all<interactions>, in<"ydot">, through<"nslaw", "e">, inout<"ydot", next>> {
  static void apply(const auto& p, const auto& ydot, const auto& e, auto& ydot_next) {
    ydot_next(0) += e * ydot(0);
    for (decltype(ydot_next.size()) i = 1; i < ydot_next.size(); ++i) {
      ydot_next(i) += (1 - p.theta) / p.theta * ydot(i);
    }
  }
};

// the free output derivative of an active interaction in q of the problem
struct assemble_free_output : kernel<active<interactions>, rank, in<"ydot", next>>, assembles {
  static void apply(auto& p, auto row, const auto& ydot_next) {
    set_value(p.vector, row, ydot_next);
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct apply_impact_law
    : stage<"Apply the Newton or the Frémond impact law on the free ydot at step + 1"> {
  static void run(auto& c) {
    if constexpr (std::decay_t<decltype(c.kind)>::nslaw_is_fremond()) {
      c.kind.template run<moreau_jean::apply_fremond_impact_law>(
          c.step, theta_method::stages::parameters(c.kind, c.h, c.integrator.theta()));
    } else {
      c.kind.template run<moreau_jean::apply_newton_impact_law>(c.step, nullptr);
    }
  }
};

// the Moreau-Jean-Frémond scheme: the Frémond impact law whatever the law of
// the interactions (a law without friction has no tangential component, its
// impact law is Newton's)
struct apply_fremond_impact_law
    : stage<"Apply the Frémond impact law on the free ydot at step + 1"> {
  static void run(auto& c) {
    c.kind.template run<moreau_jean::apply_fremond_impact_law>(
        c.step, theta_method::stages::parameters(c.kind, c.h, c.integrator.theta()));
  }
};

// for the schemes without the Frémond impact law
struct apply_newton_impact_law
    : stage<"Apply the Newton impact law on the free ydot at step + 1"> {
  static void run(auto& c) {
    if constexpr (std::decay_t<decltype(c.kind)>::nslaw_is_fremond()) {
      if (c.kind.number_of_interactions() > 0) {
        throw std::runtime_error(
            "the Frémond impact law is implemented in the Moreau-Jean scheme only");
      }
    } else {
      c.kind.template run<moreau_jean::apply_newton_impact_law>(c.step, nullptr);
    }
  }
};

struct assemble_free_velocity_output
    : stage<"Assemble q = H v_free, the free output derivatives"> {
  static void run(auto& c) {
    c.kind.template run<assemble_free_output>(
        c.step, block_vector{c.kind.lambda_view(c.integrator.q_nsp_vector_assembled())});
  }
};

// on a kind of systems

struct update_contact_state : stage<"Update v, q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    auto base = theta_method::stages::parameters(c.kind, c.h, c.integrator.theta());
    auto p = time_step{base.h, base.theta, c.kind.part(),
                       c.kind.velocity_view(c.integrator.velocity_vector_assembled())};

    // only the involved systems: the others keep their free state
    c.kind.template run<fused<update_contact_next_velocity, update_position_involved>>(c.step,
                                                                                       p);

    // the involved systems may have broken them
    c.kind.template run<theta_method::apply_boundary_conditions>(c.step, p);
  }
};

// on the integrator

// W = H M^-1 H^T with a diagonal M: the rigid systems and the interactions
// between them, the first kinds
void delassus_of_rigid_systems(auto osi) {
  auto rigid_systems = std::get<0>(osi.systems());
  auto rigid_interactions = std::get<0>(osi.interactions());

  auto tmp_matrix = rigid_interactions.transposed_h_matrix();
  auto h_mat =
      rigid_interactions.h_matrix_view(osi.h_matrix_assembled(), rigid_systems.ds_offset());
  auto m_mat = rigid_systems.mass_matrix_view(osi.mass_matrix_assembled());
  auto w_mat = rigid_interactions.w_matrix_view(osi.w_matrix_assembled());

  resize(tmp_matrix, size1(h_mat), size0(h_mat));
  if constexpr (osi.with_k_matrix()) {
    solve_linear_system_with_transpose(cast_to_diag_mat(m_mat), h_mat, tmp_matrix);
  } else {
    solve_linear_system_with_transpose(m_mat, h_mat, tmp_matrix);
  }
  prod(h_mat, tmp_matrix, w_mat);
}

// W = H A^-1 H^T
void delassus(auto osi, auto step, auto h) {
  if constexpr (osi.with_k_matrix()) {
    if (algebra::nnz(osi.k_matrix_assembled()) == 0) {
      // no FEM system in contact
      delassus_of_rigid_systems(osi);
    } else {
      algebra::add(1., osi.mass_matrix_assembled(), h * h * osi.theta() * osi.theta(),
                   osi.k_matrix_assembled(), osi.iteration_matrix_assembled());

      // the factorizations of A kept on the FEM systems
      auto blocks = osi.factored_blocks(step, h);
      compute_kkt_matrix(osi.h_matrix_assembled(), osi.iteration_matrix_assembled(),
                         osi.w_matrix_assembled(), blocks);
    }
  } else {
    delassus_of_rigid_systems(osi);
  }
}

// the velocity correction A^-1 H^T lambda (p0 = H^T lambda)
void impulse(auto osi, auto step, auto h) {
  auto& h_matrix = osi.h_matrix_assembled();
  auto& p0 = osi.p0_vector_assembled();
  auto& velocities = osi.velocity_vector_assembled();

  resize(p0, size1(h_matrix));
  resize(velocities, size1(h_matrix));

  transpose(h_matrix);
  prodt1(h_matrix, osi.lambda_vector_assembled(), p0);

  if constexpr (osi.with_k_matrix()) {
    if (algebra::nnz(osi.k_matrix_assembled()) > 0) {
      // block diagonal: the rigid systems are divided, the FEM blocks use the
      // factorizations kept on the systems
      auto blocks = osi.factored_blocks(step, h);
      solve_block_diagonal_system(osi.iteration_matrix_assembled(), p0, velocities, blocks);
      return;
    }
  }
  solve_linear_system(osi.mass_matrix_assembled(), p0, velocities);
}

struct assemble_delassus_operator
    : stage<"Assemble W = H M^-1 H^T, or H (M + h^2 theta^2 K)^-1 H^T">,
      host_code {
  static void run(auto& c) { delassus(c.integrator, c.step, c.h); }
};

struct solve : stage<"Solve the frictional contact problem (W, q, mu): lambda">, host_code {
  static void run(auto& c) {
    using time_stepping_t = typename std::decay_t<decltype(c.time_stepping)>::type;
    c.time_stepping.template solve_nonsmooth_problem<typename time_stepping_t::formulation_t>(
        c.step);
  }
};

struct apply_velocity_correction : stage<"Apply the velocity correction M^-1 H^T lambda">,
                                   host_code {
  static void run(auto& c) { impulse(c.integrator, c.step, c.h); }
};

}  // namespace stages

// the schedules

using storage::schedule::sequence;
using storage::schedule::when;

namespace th = theta_method::stages;
namespace co = contacts::stages;
namespace mj = stages;

// one stage per line: the schedules read as the algorithm
// clang-format off

// a Moreau-Jean time step, with its impact law stage
template <storage::pattern::string_literal Name, typename ImpactLaw>
using time_step_schedule_with = sequence<
    Name,
    for_systems<th::recompute_minv_f, th::compute_free_state>,
    for_interactions<co::sync_relation_sets, co::compute_jacobians, co::compute_outputs, co::compute_free_outputs>,
    co::compute_active_set,
    when<co::with_contacts,
         co::setup_assembly,
         for_interactions<co::compute_time_invariant_jacobians, co::assemble_h_matrix, co::assemble_multipliers,
                          ImpactLaw, mj::assemble_free_velocity_output>,
         for_systems<co::assemble_mass_matrix, co::assemble_stiffness_matrix>,
         mj::assemble_delassus_operator,
         mj::solve,
         mj::apply_velocity_correction,
         for_interactions<co::keep_multipliers>,
         for_systems<mj::update_contact_state>,
         for_interactions<co::update_output_velocities>>>;

// the Moreau-Jean time step: the impact law of the interactions
using time_step_schedule = time_step_schedule_with<"Moreau-Jean time step", mj::apply_impact_law>;

// the Moreau-Jean-Frémond time step: the Frémond impact law on every
// interaction, chosen when the simulation is set up (see selectable_scheme.hpp)
using fremond_time_step_schedule =
    time_step_schedule_with<"Moreau-Jean-Frémond time step", mj::apply_fremond_impact_law>;

// the state before the first step
using initial_schedule =
    sequence<"Moreau-Jean initial state", for_systems<th::set_initial_state>,
             for_interactions<co::set_initial_state>>;

// clang-format on

// the schedules and the storage they need on the kinds (see kinds.hpp)
struct scheme {
  using time_step_schedule = moreau_jean::time_step_schedule;
  using initial_schedule = moreau_jean::initial_schedule;
  template <typename System>
  using system_storage = theta_method::system_storage<System>;
  template <typename Interaction>
  using interaction_storage = contacts::interaction_storage<Interaction>;
};

}  // namespace siconos::simul::moreau_jean
