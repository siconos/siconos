/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once
// bouncing_ball split in two translation units, to measure the per-TU fixed
// cost of the world type: bouncing_ball_split.cpp sets the model up and
// bouncing_ball_split_sim.cpp instantiates the simulation methods.

#include <cstddef>

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

// defined in bouncing_ball_split_sim.cpp
void simulation_initialize(siconos::vkernel::config::data_t& data, std::size_t simulation);
unsigned simulation_compute_one_step(siconos::vkernel::config::data_t& data,
                                     std::size_t simulation);
bool simulation_has_next_event(siconos::vkernel::config::data_t& data, std::size_t simulation);
