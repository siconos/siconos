#pragma once

#include <algorithm>
#include <cmath>

#include "siconos/simul/contacts.hpp"
#include "siconos/simul/kinds.hpp"
#include "siconos/simul/moreau_jean.hpp"
#include "siconos/simul/theta_method.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul::dem {

/* The explicit DEM scheme (Cundall & Strack 1979), here for the
   comparisons with Moreau-Jean: the contacts are penalty springs, their forces explicit (at
   the state of the step), the gap allowed to close (an overlap). For a contact that overlaps,
   gap g <= 0, relative velocity (u_N, u_T):

     F_N = max(0, -k_n g - c_n u_N)                 spring and dashpot
     xi  <- xi + h u_T                              the tangential spring
     F_T = -k_t xi, at most mu F_N (Coulomb): then  xi = -F_T / k_t

   xi is the history of the contact (the tangential_spring of the
   interaction), zero when it does not overlap. k_n, c_n and k_t are the
   normal_stiffness, normal_damping and tangential_stiffness of the
   integrator.

   The contact impulse h F goes through the velocity correction of the
   Moreau-Jean scheme (moreau_jean.hpp), without a problem to solve:

     v(k + 1) = v_free + A^-1 H^T h F(q(k), v(k))
     q(k + 1) = q(k) + h (theta v(k + 1) + (1 - theta) v(k))

   With theta = 1 and rigid systems, it is the symplectic Euler scheme,
   stable for h < 2 sqrt(m / k_n) (the period of a contact spring / pi).
   The FEM systems keep the implicit theta method on their stiffness, the
   contact forces being explicit. */

using namespace kernels;

// the history of a contact: the elongation of its tangential spring
template <typename Interaction>
using tangential_spring =
    storage::attached<Interaction, symbol<"tangential_spring">, some::scalar>;

// the storage of the scheme on the interactions of a kind (the properties of
// the kind, see kinds.hpp): the tangential springs; y and ydot at step only,
// but ydot at two time levels for the contact work of the outputs (io reads
// ydot at step - 1)
template <typename Interaction>
using interaction_storage = gather<tangential_spring<Interaction>,
                                   storage::time_levels<attr_t<Interaction, "ydot">, 2>>;

// the parameters of the contact forces
template <typename H, typename S>
struct penalty {
  H h;
  S normal_stiffness;
  S normal_damping;
  S tangential_stiffness;
};

// an interaction that does not overlap: no force, no tangential spring
struct released : kernel<all<interactions>, in<"activation", by_topology>, inout<"lambda">,
                         inout<"tangential_spring", by_kind>> {
  static void apply(const auto&, bool activation, auto& lambda, auto& spring) {
    if (!activation) {
      algebra::set_zero(lambda);
      spring = 0;
    }
  }
};

// the contact force of an interaction that overlaps (y = gap <= 0, ydot =
// (u_N, u_T) at step), as an impulse h F in lambda
struct contact_force
    : kernel<active<interactions>, in<"y">, in<"ydot">, through<"nslaw", "mu">,
             inout<"tangential_spring", by_kind>, out<"lambda">> {
  static void apply(const auto& p, const auto& y, const auto& ydot, const auto& mu,
                    auto& spring, auto& lambda) {
    const auto normal = std::max(0., -p.normal_stiffness * y(0) - p.normal_damping * ydot(0));
    spring += p.h * ydot(1);
    auto tangent = -p.tangential_stiffness * spring;
    if (std::abs(tangent) > mu * normal) {
      // sliding: the spring holds the Coulomb force
      tangent = std::copysign(mu * normal, tangent);
      spring = -tangent / p.tangential_stiffness;
    }
    lambda(0) = p.h * normal;
    lambda(1) = p.h * tangent;
  }
};

// the stages (see storage/schedule.hpp)

namespace stages {

using storage::schedule::host_code;
using storage::schedule::stage;

// on the integrator

struct overlaps : stage<
                      "the active set: the interactions that overlap at step, the systems "
                      "they involve">,
                  host_code {
  static void run(auto& c) {
    // gamma = 0: the gap at step, threshold 0: an overlap
    c.involved = c.integrator.compute_active_interactions(c.step, c.h, 0., 0.);
  }
};

// on a kind of interactions

struct released : stage<"no force and no tangential spring without overlap"> {
  static void run(auto& c) { c.kind.template run<dem::released>(c.step, nullptr); }
};

struct contact_forces
    : stage<"the penalty forces at step, the tangential springs: lambda = h F"> {
  static void run(auto& c) {
    using scalar = typename decltype(c.kind.env())::scalar;
    auto osi = c.integrator;
    c.kind.template run<contact_force>(
        c.step,
        penalty<decltype(c.h), scalar>{c.h, osi.normal_stiffness(), osi.normal_damping(),
                                       osi.tangential_stiffness()});
  }
};

// the velocity correction of FEM systems needs the assembled iteration
// matrix, formed by the delassus stage in Moreau-Jean
struct iteration_matrix : stage<"A = M + h^2 theta^2 K, if FEM systems are in contact">,
                          host_code {
  static void run(auto& c) {
    auto osi = c.integrator;
    if constexpr (osi.with_k_matrix()) {
      if (algebra::nnz(osi.k_matrix_assembled()) > 0) {
        algebra::add(1., osi.mass_matrix_assembled(), c.h * c.h * osi.theta() * osi.theta(),
                     osi.k_matrix_assembled(), osi.iteration_matrix_assembled());
      }
    }
  }
};

}  // namespace stages

// the schedules

using storage::schedule::sequence;
using storage::schedule::when;

namespace th = theta_method::stages;
namespace co = contacts::stages;
namespace mj = moreau_jean::stages;
namespace de = stages;

// one stage per line: the schedules read as the algorithm
// clang-format off

// the explicit DEM time step
using time_step_schedule = sequence<
    "explicit DEM time step",
    for_systems<th::iteration_matrix, th::free_state>,
    for_interactions<co::relation_sets, co::jacobians, co::outputs>,
    de::overlaps,
    for_interactions<de::released>,
    when<co::with_contacts,
         co::assembly_setup,
         for_interactions<co::time_invariant_jacobians, de::contact_forces, co::h_matrix,
                          co::multipliers>,
         for_systems<co::mass_matrix, co::stiffness_matrix>,
         de::iteration_matrix,
         mj::velocity_correction,
         for_systems<mj::contact_update>>>;

// the state before the first step: as moreau_jean
using initial_schedule = moreau_jean::initial_schedule;

// clang-format on

struct scheme {
  using time_step_schedule = dem::time_step_schedule;
  using initial_schedule = dem::initial_schedule;
  template <typename System>
  using system_storage = theta_method::system_storage<System>;
  template <typename Interaction>
  using interaction_storage = dem::interaction_storage<Interaction>;
};

}  // namespace siconos::simul::dem
