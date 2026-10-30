#pragma once

#include <cstddef>

#include "siconos/simul/anitescu_tasora.hpp"
#include "siconos/simul/moreau_jean.hpp"
#include "siconos/simul/moreau_jean_gosi.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::selectable_scheme {

/* The scheme chosen when the simulation is set up, not when it is compiled:
   the scheme attribute of the integrator (see one_step_integrator.hpp),
   moreau_jean (local problem), moreau_jean_gosi (global problem) or
   anitescu_tasora (the reference for the comparisons). The schedules are
   all compiled, one runs. */

// the values of the scheme attribute of the integrator
enum : std::size_t { local_problem = 0, global_problem = 1, cone_complementarity = 2 };

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

using storage::schedule::sequence;
using storage::schedule::when;

struct scheme {
  using time_step_schedule =
      sequence<"time step", when<local, simul::moreau_jean::time_step_schedule>,
               when<global, simul::moreau_jean_gosi::time_step_schedule>,
               when<anitescu_tasora, simul::anitescu_tasora::time_step_schedule>>;

  // the same for all
  using initial_schedule = simul::moreau_jean::initial_schedule;
};

}  // namespace siconos::simul::selectable_scheme
