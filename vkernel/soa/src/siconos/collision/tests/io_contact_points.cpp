#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cmath>
#include <cstdio>
#include <print>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// Regression for the uninitialized-output column 15 of io::contact_points
// (interaction.y[1], the tangential gap row).
//
// compute_output in moreau_jean_element.hpp writes y[0] only on the nonlinear
// relation branches; y[1] used to be read uninitialized (NaN or 0 depending on
// the allocator) and contact_points exposed that in HDF5 column 15. The fix
// zeros y before the branch, so every row must now read y[1] == 0 and finite,
// on both the relation1 (diskfsegment_r) and relation2 (diskdisk_r) branches.
// This test drives a disks-and-ground-segment world through a real stepper and
// checks y[1] on the active interactions each step.

namespace siconos::config {
namespace some = siconos::storage::some;
struct disk : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
using solver_options = simul::solver_options;
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r> {};
struct topo : simul::topology<disk, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<disk, collision::empty_shape>;
using pointl = collision::point<storage::pattern::empty_item, collision::shape::segment>;
struct neighborhood : collision::neighborhood<pointd, pointl> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t : storage::make_runtime_database<
                    env, simulation, interaction_manager, segment_shape, disk_shape,
                    storage::with_properties<
                        storage::wrapped<disk, some::unbounded_collection>,
                        storage::wrapped<diskdisk_r, some::unbounded_collection>,
                        storage::wrapped<diskfsegment_r, some::unbounded_collection>,
                        storage::wrapped<interaction, some::unbounded_collection>,
                        storage::wrapped<segment_shape, some::unbounded_collection>,
                        storage::wrapped<disk_shape, some::unbounded_collection>,
                        storage::attached<disk, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<disk_shape>>,
                        storage::time_invariant<storage::pattern::attr_t<disk, "fext">>,
                        storage::diagonal<storage::pattern::attr_t<disk, "mass_matrix">>,
                        storage::assembled_diagonal<storage::pattern::attr_t<
                            typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {

class IoContactPointsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(IoContactPointsTest);
  CPPUNIT_TEST(testIoContactPoints);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testIoContactPoints() {
    using namespace siconos;
    using namespace storage;

    config::data_t data;

    double t0 = 0.;
    double tmax = 0.7;
    double h = 0.005;
    double position_init = 1.;
    double theta = 0.5;
    double radius = 0.1;
    double m = 1.;
    double g = 9.81;
    double e = 0.9;

    std::size_t ndisks = 2;

    auto disk_shape = storage::add<config::disk_shape>(data);
    disk_shape.radius() = radius;

    using disk_handle_t = decltype(storage::add<config::disk>(data));
    std::vector<disk_handle_t> disks;
    disks.reserve(ndisks);
    for (std::size_t i = 0; i < ndisks; ++i) disks.push_back(storage::add<config::disk>(data));

    for (std::size_t i = 0; i < disks.size(); ++i) {
      auto& disk = disks[i];
      prop<"id", config::topo>(disk) = static_cast<uint>(i + 1);
      disk.q() = {0, position_init * static_cast<uint>(i + 1), 0};
      disk.velocity() = {0, 0., 0.};
      disk.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
      prop<"shape">(disk) = disk_shape;
      disk.fext() = {0., -m * g, 0.};
    }

    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = e;

    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

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

    auto ngbh = storage::add<config::neighborhood>(data);
    ngbh.create(0.6);

    auto segment = storage::add<config::segment_shape>(data);
    auto diskdisk_r = storage::add<config::diskdisk_r>(data);
    auto ground_r = storage::add<config::diskfsegment_r>(data);
    ground_r.segment() = segment;
    segment.x1() = -10.;
    segment.y1() = 0.;
    segment.x2() = 10.;
    segment.y2() = 0.;
    segment.initialize();
    segment.maxpoints() = 10;

    auto spacef = storage::add<config::space_filter>(data);
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = diskdisk_r;
    spacef.insert_nonsmooth_law(nslaw, 0, 0);
    spacef.insert_diskfsegment_r(ground_r);
    spacef.make_points();
    ngbh.add_point_sets(0);

    simulation.initialize();

    std::size_t checks = 0;

    while (simulation.has_next_event()) {
      ngbh.update(0);
      ngbh.search();
      spacef.update_index_set0(simulation.current_step());

      simulation.compute_one_step();

      // compute_output (called inside compute_one_step) zeroes y before the
      // nonlinear relation branches, so y[1] must be defined (==0, finite) for
      // every active interaction this step.
      auto& ys =
          storage::attr_values<config::interaction, "y">(data, simulation.current_step());
      for (auto [y] : view::zip(ys)) {
        ++checks;
        CPPUNIT_ASSERT_MESSAGE("y[1] (tangential gap row) is zero and finite",
                               std::isfinite(y[0]) && std::isfinite(y[1]) && y[1] == 0.);
      }
    }

    // Require at least one contact to have been produced, so the check is not
    // vacuous: disks fall onto the ground segment within the window.
    CPPUNIT_ASSERT_MESSAGE("at least one interaction was observed", checks > 0);

    std::print("io_contact_points: {} y-rows checked\n", checks);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(IoContactPointsTest);

}  // namespace
