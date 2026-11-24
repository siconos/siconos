// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "RtOsiTest.hpp"

#include "siconos/vkernel/config/config.hpp"
#include "siconos/vkernel/config/environment.hpp"
#include "siconos/vkernel/model/fem.hpp"
#include "siconos/vkernel/model/lagrangian_ds.hpp"
#include "siconos/vkernel/model/nslaws.hpp"
#include "siconos/vkernel/simul/interaction.hpp"
#include "siconos/vkernel/simul/one_step_integrator.hpp"
#include "siconos/vkernel/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/vkernel/simul/time_discretization.hpp"
#include "siconos/vkernel/simul/time_stepping.hpp"
#include "siconos/vkernel/simul/topology.hpp"
#include "siconos/vkernel/storage/mp/mp.hpp"

namespace siconos::vkernel::config {
struct fem_ds : model::lagrangian_ds {};
struct ball : model::lagrangian_ds {};
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_r<nslaw::size>;
using interaction = simul::interaction<nslaw, relation>;
using fem_interaction = simul::rt_ct_interaction<nslaw, relation>;
using topo = simul::topology<ball, interaction, fem_ds, fem_interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using simulation = simul::time_stepping<td, osi, osnspb>;

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>>;
};

struct make
    : storage::make_runtime_database<
          config::env, fem_ds, simulation,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<ball, "fext">>,
              storage::diagonal<storage::attr_t<ball, "mass_matrix">>,
              storage::assembled_diagonal<storage::attr_t<osi, "mass_matrix_assembled">>>> {};

}  // namespace siconos::vkernel::config

namespace config = siconos::vkernel::config;
namespace store = siconos::vkernel::storage;
namespace mp = store::mp;
using siconos::vkernel::algebra::vector;

CPPUNIT_TEST_SUITE_REGISTRATION(RtOsiTest);

void RtOsiTest::setUp() {}

void RtOsiTest::tearDown() {}

// check that the kinds of the integrator exist and that their offsets can be
// modified
void RtOsiTest::testOsi0() {
  auto data = siconos::vkernel::config::make();

  auto osi = store::add<config::osi>(data);

  std::size_t i = 0;
  std::size_t j = 0;
  mp::for_each(osi.systems(), [&](auto kind) { kind.ds_offset() = i++; });
  mp::for_each(osi.interactions(), [&](auto kind) { kind.inter_offset() = j++; });

  CPPUNIT_ASSERT(i == 2);
  CPPUNIT_ASSERT(j == 2);

  i = 0;
  j = 0;
  mp::for_each(osi.systems(), [&](auto kind) { CPPUNIT_ASSERT((kind.ds_offset() == i++)); });
  mp::for_each(osi.interactions(),
               [&](auto kind) { CPPUNIT_ASSERT((kind.inter_offset() == j++)); });
}
