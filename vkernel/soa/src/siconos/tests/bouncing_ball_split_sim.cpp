/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "bouncing_ball_split.hpp"

namespace {
auto simulation_handle(siconos::config::data_t& data, std::size_t simulation) {
  return siconos::storage::make_handle(
      data, siconos::storage::index<siconos::config::simulation, std::size_t>(simulation));
}
}  // namespace

void simulation_initialize(siconos::config::data_t& data, std::size_t simulation) {
  simulation_handle(data, simulation).initialize();
}

unsigned simulation_compute_one_step(siconos::config::data_t& data, std::size_t simulation) {
  return simulation_handle(data, simulation).compute_one_step();
}

bool simulation_has_next_event(siconos::config::data_t& data, std::size_t simulation) {
  return simulation_handle(data, simulation).has_next_event();
}
