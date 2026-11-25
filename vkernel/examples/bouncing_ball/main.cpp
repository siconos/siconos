// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

// The whole point of this file: an installed vkernel is used through two
// include strings and nothing from the source tree.
#include <print>

#include "siconos/vkernel.hpp"
#include "siconos/vkernel/config/config_builder.hpp"

namespace siconos::vkernel::config {
struct ball : model::lagrangian_ds {};
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<lcp> {};
struct nslaw : model::newton_impact {};
struct relation : model::lagrangian_r<nslaw::size> {};
struct interaction : simul::interaction<nslaw, relation> {};
struct topo : simul::topology<ball, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

struct data_t : decltype(storage()
                             .with_param<"dof", 3>()
                             .with_items<simulation>()
                             .with_time_invariant<ball, "fext">()
                             .with_diagonal<ball, "mass_matrix">()
                             .with_assembled_diagonal<osi, "mass_matrix_assembled">()
                             .build_runtime_database()){};
}  // namespace siconos::vkernel::config

int main() {
  namespace storage = siconos::vkernel::storage;
  namespace config = siconos::vkernel::config;
  using siconos::vkernel::storage::handle;

  config::data_t data;

  double t0 = 0;               // initial computation time
  double tmax = 10;            // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Ball radius
  double m = 1.;               // Ball mass
  double g = 9.81;             // Gravity

  std::print("====> Model loading ...\n");

  // --------------------------
  // -- The dynamical_system --
  // --------------------------
  handle ball = storage::add<config::ball>(data);

  ball.q() = {position_init, 0, 0};
  ball.velocity() = {velocity_init, 0, 0};
  ball.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

  // -- Set external forces (weight) --
  ball.fext() = {-m * g, 0., 0.};

  // ------------------
  // -- The relation --
  // ------------------

  // -- Lagrangian relation --
  handle relation = storage::add<config::relation>(data);
  relation.h_matrix() = {-1.0, 0., 0.};

  // -- nslaw --
  double e = 0.9;
  handle nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  handle lcp = storage::add<config::lcp>(data);
  lcp.create(SICONOS_LCP_LEMKE);

  // ------------------
  // --- Simulation ---
  // ------------------
  handle simulation = storage::add<config::simulation>(data);

  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;

  simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  handle osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = lcp;

  // -- set the options --
  handle so = storage::add<siconos::vkernel::simul::solver_options>(data);
  so.create();
  osnspb.options() = so;

  // Interaction ball-floor
  handle interaction = simulation.topology().link(ball);
  simulation.topology().bind(interaction, relation);
  interaction.nslaw() = nslaw;

  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  // fix this for constant fext
  simulation.initialize();

  while (simulation.has_next_event()) {
    simulation.compute_one_step();
  }

  std::print("t = {:.6f}  q = {:.15e}  v = {:.15e}\n",
             simulation.current_step() * simulation.time_step(),
             storage::attr<"q">(ball, simulation.current_step())(0),
             storage::attr<"velocity">(ball, simulation.current_step())(0));
  return 0;
}
