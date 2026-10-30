#pragma once

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
     q  = H v_free + e ydot(k)         the free output derivatives
     lambda solution of the frictional contact problem (W, q, mu)
     v(k + 1) = v_free + A^-1 H^T lambda
     q(k + 1) = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   The kernels of the scheme, its stages, then the schedules (the time step,
   the initial state). The rest of a step is in contacts.hpp. */

using namespace kernels;
using theta_method::involved_position;
using theta_method::time_step;

// an assembled vector
template <typename Vector>
struct block_vector {
  Vector vector;
};

// the contact impulse of an involved system
struct contact_velocity : kernel<involved<systems>, rank, inout<"velocity", next>> {
  static void apply(auto& p, auto index, auto& v_next) {
    v_next += p.part(p.velocities, index, v_next);
  }
};

// the Newton impact law on the free output derivative
struct impact_law
    : kernel<all<interactions>, in<"ydot">, through<"nslaw", "e">, inout<"ydot", next>> {
  static void apply(const auto&, const auto& ydot, const auto& e, auto& ydot_next) {
    ydot_next += e * ydot;
  }
};

// the free output derivative of an active interaction in q of the problem
struct free_output : kernel<active<interactions>, rank, in<"ydot", next>>, assembles {
  static void apply(auto& p, auto row, const auto& ydot_next) {
    set_value(p.vector, row, ydot_next);
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::host_code;
using storage::schedule::stage;

// on a kind of interactions

struct impact_law : stage<"the Newton impact law on the free ydot at step + 1"> {
  static void run(auto& c) { c.kind.template run<moreau_jean::impact_law>(c.step, nullptr); }
};

struct free_velocity_output : stage<"q = H v_free, the free output derivatives"> {
  static void run(auto& c) {
    c.kind.template run<free_output>(
        c.step, block_vector{c.kind.lambda_view(c.integrator.q_nsp_vector_assembled())});
  }
};

// on a kind of systems

struct contact_update : stage<"v, q at step + 1 of the involved systems"> {
  static void run(auto& c) {
    auto base = theta_method::stages::parameters(c.kind, c.h, c.integrator.theta());
    auto p = time_step{base.h, base.theta, c.kind.part(),
                       c.kind.velocity_view(c.integrator.velocity_vector_assembled())};

    // only the involved systems: the others keep their free state
    c.kind.template run<fused<contact_velocity, involved_position>>(c.step, p);

    // the involved systems may have broken them
    c.kind.template run<theta_method::boundary_conditions>(c.step, p);
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

struct delassus_operator : stage<"W = H M^-1 H^T, or H (M + h^2 theta^2 K)^-1 H^T">,
                           host_code {
  static void run(auto& c) { delassus(c.integrator, c.step, c.h); }
};

struct solve : stage<"the frictional contact problem (W, q, mu): lambda">, host_code {
  static void run(auto& c) {
    using time_stepping_t = typename std::decay_t<decltype(c.time_stepping)>::type;
    c.time_stepping.template solve_nonsmooth_problem<typename time_stepping_t::formulation_t>(
        c.step);
  }
};

struct velocity_correction : stage<"the velocity correction M^-1 H^T lambda">, host_code {
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

// the Moreau-Jean time step
using time_step_schedule = sequence<
    "Moreau-Jean time step",
    for_systems<th::iteration_matrix, th::free_state>,
    for_interactions<co::relation_sets, co::jacobians, co::outputs, co::free_outputs>,
    co::active_set,
    when<co::with_contacts,
         co::assembly_setup,
         for_interactions<co::time_invariant_jacobians, co::h_matrix, co::multipliers,
                          mj::impact_law, mj::free_velocity_output>,
         for_systems<co::mass_matrix, co::stiffness_matrix>,
         mj::delassus_operator,
         mj::solve,
         mj::velocity_correction,
         for_interactions<co::keep_multipliers>,
         for_systems<mj::contact_update>>>;

// the state before the first step
using initial_schedule =
    sequence<"Moreau-Jean initial state", for_systems<th::initial_state>,
             for_interactions<co::initial_state>>;

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
