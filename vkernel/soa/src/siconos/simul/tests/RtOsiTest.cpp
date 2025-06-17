
#include "RtOsiTest.hpp"

#include "siconos/config/config.hpp"
#include "siconos/model/fem.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include "siconos/storage/ground/ground.hpp"
#include "siconos/utils/environment.hpp"

namespace siconos::config {
using fem = model::finite_element_linear_tids;
using fem_ds = model::rt_lagrangian_ds;
using material = model::material;
using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_r<nslaw::size>;
using interaction = simul::interaction<nslaw, relation>;
using fem_interaction = simul::rt_interaction<nslaw, relation>;
using topo = simul::topology<ball, interaction, fem_ds, fem_interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using simulation = simul::time_stepping<td, osi, osnspb>;

using params = map<iparam<"dof", 3>>;

struct make
    : storage::make<
          standard_environment<params>, fem_ds, fem, material, simulation,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<ball, "fext">>,
              storage::diagonal<storage::attr_t<ball, "mass_matrix">>,
              storage::assembled_diagonal<storage::attr_t<
                  typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {
};

}  // namespace siconos::config

namespace config = siconos::config;
namespace store = siconos::storage;
namespace ct = store::ground;
using siconos::algebra::vector;

CPPUNIT_TEST_SUITE_REGISTRATION(RtOsiTest);

void RtOsiTest::setUp() {}

void RtOsiTest::tearDown() {}

// check that osi elements exist and that offsets can be modified
void RtOsiTest::testOsi0()
{
  auto data = siconos::config::make();

  auto osi = store::add<config::osi>(data);

  std::size_t i = 0;
  std::size_t j = 0;
  ct::for_each(osi.elements(), [&](auto elem) {
    elem.ds_offset() = i++;
    elem.inter_offset() = j++;
  });

  i = 0;
  j = 0;
  ct::for_each(osi.elements(), [&](auto elem) {
    CPPUNIT_ASSERT((elem.ds_offset() == i++));
    CPPUNIT_ASSERT((elem.inter_offset() == i++));
  });
}
