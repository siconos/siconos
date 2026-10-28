/* End-to-end test of the cxd-backed narrowphase stage.

   Mirrors tests/bouncing_disks.cpp: two bodies fall, collide and separate
   while the simulation runs. The ONLY addition to the bouncing_disks
   harness loop is the np.ensure(step) line, which rebinds the generic
   ds-ds relations to polytopepolytope_r items and fills their cached
   status/h/normal/cp1/cp2 from one batched narrowphase computation.

   The bodies are side-2 axis-aligned squares; body 1 falls on body 2
   (body 2 carries no gravity, so the gap really closes). The contact is
   perfectly aligned, so no torque is expected: the rotation stays ~0 and
   the cached h must match the analytic axis-aligned gap |y2 - y1| - 2
   while the pair is separated. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>
#include <cmath>
#include <limits>
#include <print>
#include <stdexcept>
#include <string>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
/* the relation + stage headers must be visible before the config names
   them; the load-bearing cxd include order lives inside cxd_narrowphase.hpp */
#include "siconos/collision/polytopepolytope_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;

struct polytope : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using diskdisk_r = collision::diskdisk_r;
/* topo is the acyclic marker of the relation/narrowphase pair; defined
   right after the interaction that names it */
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, polytopepolytope_r> {};
struct topo : simul::topology<polytope, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<polytope, collision::empty_shape>;
struct neighborhood : collision::neighborhood<pointd> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t : storage::make_runtime_database<
                    env, simulation, interaction_manager, polytope, polytope_shape, diskdisk_r,
                    nslaw, interaction, topo, polytopepolytope_r, cxd_narrowphase, osi, td,
                    osnspb, neighborhood, space_filter,
                    storage::with_properties<
                        storage::wrapped<polytope, some::unbounded_collection>,
                        storage::wrapped<diskdisk_r, some::unbounded_collection>,
                        storage::wrapped<pointd, some::unbounded_collection>,
                        storage::wrapped<interaction, some::unbounded_collection>,
                        storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
                        storage::attached<polytope, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<polytope_shape>>,
                        storage::time_invariant<storage::pattern::attr_t<polytope, "fext">>,
                        storage::diagonal<storage::pattern::attr_t<polytope, "mass_matrix">>,
                        storage::assembled_diagonal<storage::pattern::attr_t<
                            typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {

class BouncingPolytopesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(BouncingPolytopesTest);
  CPPUNIT_TEST(testBouncingPolytopes);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testBouncingPolytopes() {
    using namespace siconos;
    using namespace storage;

    config::data_t data;

    const double t0 = 0.;
    const double tmax = 2.;  // 400 steps of dt = 0.005
    const double h = 0.005;
    const double theta = 0.5;  // Moreau-Jean integrator parameter
    const double m = 1.;
    const double g = 9.81;
    const double e = 0.9;

    // one square shape shared by both bodies (side 2, CCW)
    auto shape = storage::add<config::polytope_shape>(data);
    const std::array<std::array<double, 2>, 4> square = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
    for (int i = 0; i < 4; ++i) shape.vertices()[i] = {square[i][0], square[i][1]};
    shape.nvertices() = 4;

    // body 1 falls on body 2; body 2 carries no gravity (it is the target,
    // so the gap really closes -- with the same gravity on both the gap
    // would never change)
    auto body1 = storage::add<config::polytope>(data);
    body1.q() = {0, 5, 0};
    body1.velocity() = {0, 0, 0};
    body1.mass_matrix().diagonal() << m, m, 2. / 3.;  // side-2 square: I = (a^2+a^2)/12
    storage::prop<"shape">(body1) = shape;
    body1.fext() = {0., -m * g, 0.};

    auto body2 = storage::add<config::polytope>(data);
    body2.q() = {0, 0, 0};
    body2.velocity() = {0, 0, 0};
    body2.mass_matrix().diagonal() << m, m, 2. / 3.;
    storage::prop<"shape">(body2) = shape;
    body2.fext() = {0., 0., 0.};

    // -- the relation --
    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = e;

    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

    // -- the simulation --
    auto simulation = storage::add<config::simulation>(data);
    simulation.one_step_integrator().theta() = theta;
    simulation.one_step_integrator().constraint_activation_threshold() = 0.;
    simulation.time_discretization().t0() = t0;
    simulation.time_discretization().h() = h;
    simulation.time_discretization().tmax() = tmax;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;

    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_LEMKE);
    osnspb.options() = so;

    // -- the broadphase and the filter --
    auto ngbh = storage::add<config::neighborhood>(data);
    // the radius must exceed the initial centre distance (5) so that the
    // pair is a broadphase neighbour from the very first step on
    ngbh.create(5.5);

    auto dsds = storage::add<config::diskdisk_r>(data);  // generic ds-ds binder
    auto spacef = storage::add<config::space_filter>(data);
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = dsds;
    spacef.insert_nonsmooth_law(nslaw, 0, 0);
    spacef.make_points();
    ngbh.add_point_sets(0);

    // the ds-ds interaction, created the way space_filter::dsds_activation
    // does (the filter binds ds-ds pairs through point-SET pairs, so a world
    // with a single point set -- no ground here -- never produces one): the
    // relation starts on the generic binder, the stage rebinds it.
    // It is created after simulation.initialize(), like the filter creates its
    // interactions inside the loop, so that the step columns know it.
    auto topo_h = spacef.topology();

    // -- the cxd narrowphase stage --
    auto np = storage::add<config::cxd_narrowphase>(data);
    np.dsds_r() = dsds;

    simulation.initialize();

    // one interaction per point of the contact manifold (two in 2D), as
    // space_filter::dsds_activation creates them for a polytope pair; the
    // stage gives them manifold points 0 and 1
    for (int pi = 0; pi < 2; ++pi) {
      auto inter = topo_h.link(body1, body2);
      inter.nslaw() = nslaw;
      topo_h.bind(inter, dsds);
    }

    bool have_prev = false;
    double prev_h = 0.;
    double prev_vrel = 0.;
    double min_gap = std::numeric_limits<double>::infinity();
    double max_dvrel = 0.;
    std::size_t nsteps = 0;

    while (simulation.has_next_event()) {
      const auto step = simulation.current_step();
      // the broadphase is exercised, but space_filter::update_index_set0 is
      // NOT called here: this world has no ground, so the filter has nothing
      // to activate, and its remove_interactions assumes every interaction
      // comes from its own maps (it would dereference an end iterator on the
      // hand-created ds-ds interaction above).
      ngbh.update(step);
      ngbh.search();
      np.ensure(step);  // <- the only addition to the bouncing_disks loop

      const double v1 = storage::attr<"velocity">(body1, step)(1);
      const double v2 = storage::attr<"velocity">(body2, step)(1);
      const double vrel = v1 - v2;
      const double y1 = storage::attr<"q">(body1, step)(1);
      const double y2 = storage::attr<"q">(body2, step)(1);
      const double theta1 = storage::attr<"q">(body1, step)(2);
      const double theta2 = storage::attr<"q">(body2, step)(2);

      for (auto rel : storage::handles<config::polytopepolytope_r>(data, step)) {
        const double hval = rel.h();
        CPPUNIT_ASSERT_MESSAGE("cached h is finite", std::isfinite(hval));
        // a perfectly aligned normal impact must not excite any torque
        CPPUNIT_ASSERT_MESSAGE("no torque on the bodies",
                               std::abs(theta1) <= 1e-6 && std::abs(theta2) <= 1e-6);
        if (have_prev) {
          // the gap can only change with the relative centre velocity
          const double bound = 2. * h * std::max(std::abs(vrel), std::abs(prev_vrel));
          CPPUNIT_ASSERT_MESSAGE("h changes within the velocity bound",
                                 std::abs(hval - prev_h) <= bound + 1e-9);
          // a contact exchanges an impulse: the relative velocity jumps
          max_dvrel = std::max(max_dvrel, std::abs(vrel - prev_vrel));
        }
        // while the pair is separated, the cached h is the analytic gap
        const double ref_gap = std::abs(y2 - y1) - 2.;
        if (ref_gap > 0.)
          CPPUNIT_ASSERT_MESSAGE("h matches the analytic axis-aligned gap",
                                 std::abs(hval - ref_gap) <= 1e-6);
        min_gap = std::min(min_gap, ref_gap);
        prev_h = hval;
        prev_vrel = vrel;
        have_prev = true;
      }

      simulation.compute_one_step();
      ++nsteps;
    }

    CPPUNIT_ASSERT_MESSAGE("the simulation ran for the expected number of steps",
                           nsteps > 300);
    // The contact law is predictive (Moreau-Jean activates the constraint on
    // the predicted gap), so the solver keeps the bodies apart: the narrowphase
    // h stays positive and the status stays "separated" -- the impact itself
    // (the impulse) is what proves the contact was processed.
    CPPUNIT_ASSERT_MESSAGE("the bodies came into contact range", min_gap < 0.05);
    CPPUNIT_ASSERT_MESSAGE("the contact impulse was applied (relative velocity jump)",
                           max_dvrel > 1.);
    CPPUNIT_ASSERT_MESSAGE("no NaN in the final velocities",
                           std::isfinite(storage::attr<"velocity">(body1, nsteps)(1)) &&
                               std::isfinite(storage::attr<"velocity">(body2, nsteps)(1)));

    std::print(
        "bouncing_polytopes: {} steps, contact processed (min gap {:.6f}), h stays "
        "consistent\n",
        nsteps, min_gap);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(BouncingPolytopesTest);

}  // namespace
