#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <fstream>
#include <print>
#include <stdexcept>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"

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
struct io : siconos::io::io<osi> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t : storage::make_compile_time_database<
                    env, simulation, interaction_manager, io, segment_shape, disk_shape,
                    storage::with_properties<
                        storage::wrapped<disk, some::unbounded_collection>,
                        storage::wrapped<diskdisk_r, some::unbounded_collection>,
                        storage::wrapped<diskfsegment_r, some::unbounded_collection>,
                        storage::wrapped<pointl, some::unbounded_collection>,
                        storage::wrapped<pointd, some::unbounded_collection>,
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

class NslawsMultiGroupTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(NslawsMultiGroupTest);
  CPPUNIT_TEST(testOneInsertOrder);
  CPPUNIT_TEST(testUnregisteredPairThrows);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testOneInsertOrder() {
    using namespace siconos;
    using namespace storage;

    // ============================================================
    // Part A: two collision groups, laws inserted in ONE order only.
    // The mirror insert in insert_nonsmooth_law covers both (a,b) and
    // (b,a), so the non-canonical visit order is not a failure path.
    // ============================================================
    {
      config::data_t data;

      double t0 = 0, tmax = 0.5, h = 0.005, theta = 0.5;
      double radius = 0.1, m = 1., g = 9.81, e = 0.9;

      auto disk_shape = storage::add<config::disk_shape>(data);
      disk_shape.radius() = radius;

      // disk 1: group 0, disk 2: group 1, disk 3: group 1
      auto d1 = storage::add<config::disk>(data);
      auto d2 = storage::add<config::disk>(data);
      auto d3 = storage::add<config::disk>(data);
      d1.collision_group() = 0;
      d2.collision_group() = 1;
      d3.collision_group() = 1;

      d1.q() = {0, 1.5, 0};
      d2.q() = {0, 0.5, 0};
      d3.q() = {2.0, 1.0, 0};
      for (auto d : {std::ref(d1), std::ref(d2), std::ref(d3)}) {
        d.get().velocity() = {0, 0, 0};
        d.get().mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
        prop<"shape">(d.get()) = disk_shape;
        d.get().fext() = {0., -m * g, 0.};
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

      auto diskdisk_r = storage::add<config::diskdisk_r>(data);
      auto ground_r = storage::add<config::diskfsegment_r>(data);
      auto segment = storage::add<config::segment_shape>(data);
      ground_r.segment() = segment;
      segment.x1() = -10.;
      segment.y1() = 0.;
      segment.x2() = 10.;
      segment.y2() = 0.;
      segment.initialize();
      segment.maxpoints() = 10;
      segment.collision_group() = 1;

      auto spacef = storage::add<config::space_filter>(data);
      spacef.neighborhood() = ngbh;
      spacef.diskdisk_r() = diskdisk_r;
      // ONE order only: (0,1) for disk-disk, (1,1) for disk-segment.
      // The mirror insert covers (1,0) and the disk1/segment pair (0,1).
      spacef.insert_nonsmooth_law(nslaw, 0, 1);
      spacef.insert_nonsmooth_law(nslaw, 1, 1);
      spacef.insert_diskfsegment_r(ground_r);
      spacef.make_points();
      ngbh.add_point_sets(0);

      simulation.initialize();

      int steps = 0;
      while (simulation.has_next_event()) {
        ngbh.update(0);
        ngbh.search();
        spacef.update_index_set0(simulation.current_step());
        simulation.compute_one_step();
        ++steps;
      }
      std::print("Part A OK: {} steps, no exception\n", steps);
    }
  }

  void testUnregisteredPairThrows() {
    using namespace siconos;
    using namespace storage;

    // ============================================================
    // Part B: a colliding pair with NO registered law must fail
    // loudly (Issue 4). Disk in group 0 and disk in group 2 overlap;
    // no law for (0,2) or (2,0) is inserted, so update_index_set0
    // throws std::runtime_error from the space filter path.
    // ============================================================
    {
      config::data_t data;

      double radius = 0.1, m = 1., g = 9.81, e = 0.9;

      auto disk_shape = storage::add<config::disk_shape>(data);
      disk_shape.radius() = radius;

      auto d1 = storage::add<config::disk>(data);
      auto d2 = storage::add<config::disk>(data);
      d1.collision_group() = 0;
      d2.collision_group() = 2;  // no law registered for (0,2) or (2,0)

      // overlap on purpose -> guaranteed neighbor pair
      d1.q() = {0, 1.0, 0};
      d2.q() = {0, 1.0, 0};
      d1.velocity() = {0, 0, 0};
      d2.velocity() = {0, 0, 0};
      d1.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
      d2.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;
      prop<"shape">(d1) = disk_shape;
      prop<"shape">(d2) = disk_shape;
      d1.fext() = {0., -m * g, 0.};
      d2.fext() = {0., -m * g, 0.};

      auto nslaw = storage::add<config::nslaw>(data);
      nslaw.e() = e;

      auto fc2d = storage::add<config::fc2d>(data);
      fc2d.create();
      fc2d.instance()->dimension = 2;
      fc2d.instance()->mu = 0;

      auto simulation = storage::add<config::simulation>(data);
      simulation.one_step_integrator().theta() = 0.5;
      simulation.one_step_integrator().constraint_activation_threshold() = 0.;
      simulation.time_discretization().t0() = 0;
      simulation.time_discretization().h() = 0.005;
      simulation.time_discretization().tmax() = 0.1;
      auto osnspb = simulation.one_step_nonsmooth_problem();
      osnspb.problem() = fc2d;
      auto so = storage::add<simul::solver_options>(data);
      so.create(SICONOS_FRICTION_2D_LEMKE);
      osnspb.options() = so;

      auto ngbh = storage::add<config::neighborhood>(data);
      ngbh.create(0.6);

      auto diskdisk_r = storage::add<config::diskdisk_r>(data);

      auto spacef = storage::add<config::space_filter>(data);
      spacef.neighborhood() = ngbh;
      spacef.diskdisk_r() = diskdisk_r;
      // intentionally NO law for (0,2) / (2,0)
      spacef.make_points();
      ngbh.add_point_sets(0);

      simulation.initialize();

      bool threw = false;
      try {
        ngbh.update(0);
        ngbh.search();
        spacef.update_index_set0(simulation.current_step());
      } catch (const std::runtime_error& e) {
        threw = true;
        std::print("Part B PASS: {}\n", e.what());
      }
      CPPUNIT_ASSERT_MESSAGE("no exception thrown for an unregistered pair", threw);
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(NslawsMultiGroupTest);

}  // namespace
