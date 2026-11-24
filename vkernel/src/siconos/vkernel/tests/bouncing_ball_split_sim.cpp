/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "bouncing_ball_split.hpp"

namespace {
auto simulation_handle(siconos::vkernel::config::data_t& data, std::size_t simulation) {
  return siconos::vkernel::storage::make_handle(
      data,
      siconos::vkernel::storage::index<siconos::vkernel::config::simulation, std::size_t>(
          simulation));
}
}  // namespace

void simulation_initialize(siconos::vkernel::config::data_t& data, std::size_t simulation) {
  simulation_handle(data, simulation).initialize();
}

unsigned simulation_compute_one_step(siconos::vkernel::config::data_t& data,
                                     std::size_t simulation) {
  return simulation_handle(data, simulation).compute_one_step();
}

bool simulation_has_next_event(siconos::vkernel::config::data_t& data,
                               std::size_t simulation) {
  return simulation_handle(data, simulation).has_next_event();
}
