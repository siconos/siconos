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

#pragma once

#include <cstddef>

#include "siconos/vkernel/simul/anitescu_tasora.hpp"
#include "siconos/vkernel/simul/dem.hpp"
#include "siconos/vkernel/simul/moreau_jean.hpp"
#include "siconos/vkernel/simul/moreau_jean_gosi.hpp"
#include "siconos/vkernel/storage/schedule.hpp"

namespace siconos::vkernel::simul::selectable_scheme {

/* The scheme chosen when the simulation is set up, not when it is compiled:
   the scheme attribute of the integrator (see one_step_integrator.hpp),
   moreau_jean (local problem), moreau_jean_gosi (global problem),
   anitescu_tasora (two solvers), dem, or moreau_jean with the Frémond impact
   law on every interaction (local problem). The schedules are all compiled,
   one runs. */

// the values of the scheme attribute of the integrator
enum : std::size_t {
  local_problem = 0,
  global_problem = 1,
  cone_complementarity = 2,       // Anitescu-Tasora, projected Gauss-Seidel
  cone_complementarity_apgd = 3,  // Anitescu-Tasora, APGD
  explicit_penalty = 4,           // explicit DEM
  local_problem_fremond = 5,      // Moreau-Jean-Frémond
};

using storage::schedule::condition;

struct local : condition<"the scheme is Moreau-Jean on the local problem"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == local_problem;
  }
};

struct global : condition<"the scheme is Moreau-Jean on the global problem"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == global_problem;
  }
};

struct anitescu_tasora : condition<"the scheme is Anitescu-Tasora"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == cone_complementarity;
  }
};

struct anitescu_tasora_apgd : condition<"the scheme is Anitescu-Tasora with APGD"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == cone_complementarity_apgd;
  }
};

struct dem : condition<"the scheme is the explicit DEM"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == explicit_penalty;
  }
};

struct local_fremond
    : condition<
          "the scheme is Moreau-Jean on the local problem, with the Frémond impact law"> {
  static bool holds(const auto& c) {
    auto osi = c.integrator;
    return osi.scheme() == local_problem_fremond;
  }
};

using storage::schedule::sequence;
using storage::schedule::when;

struct scheme {
  using time_step_schedule =
      sequence<"time step", when<local, simul::moreau_jean::time_step_schedule>,
               when<global, simul::moreau_jean_gosi::time_step_schedule>,
               when<anitescu_tasora, simul::anitescu_tasora::time_step_schedule>,
               when<anitescu_tasora_apgd, simul::anitescu_tasora::apgd_time_step_schedule>,
               when<dem, simul::dem::time_step_schedule>,
               when<local_fremond, simul::moreau_jean::fremond_time_step_schedule>>;

  // the same for all
  using initial_schedule = simul::moreau_jean::initial_schedule;

  // the storage of all the schemes: the theta-method on the systems, the
  // contacts and the tangential springs of dem on the interactions
  template <typename System>
  using system_storage = theta_method::system_storage<System>;
  template <typename Interaction>
  using interaction_storage =
      decltype(storage::mp::concat(contacts::interaction_storage<Interaction>{},
                                   gather<simul::dem::tangential_spring<Interaction>>{}));
};

}  // namespace siconos::vkernel::simul::selectable_scheme
