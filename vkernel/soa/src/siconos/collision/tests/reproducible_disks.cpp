#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#ifdef WITH_GPU
#include "cuNSearch.h"
namespace CompactNSearch = cuNSearch;
#else
#include "CompactNSearch/CompactNSearch.h"
#endif

#include <algorithm>
#include <cmath>
#include <format>
#include <print>
#include <string>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

/* Reproducibility with a multithreaded neighbor search (run with several
   OpenMP threads, see vkernel_tests.cmake). The space filter creates the
   contacts in the order of the neighbors and the NSGS solution depends on the
   order of the contacts, so this order must not depend on the scheduling:

   - the neighbor lists of CompactNSearch are sorted, whatever the number of
     threads (the contract, checked on a large point cloud: a scheduling
     dependent order is not sorted, and this check does not depend on the
     scheduling of the run);
   - the same granular flow (a grid of disks with friction sliding on a slope,
     in a box) run twice gives bitwise identical states. */

namespace siconos::config {
namespace some = siconos::storage::some;
struct disk : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
using translated_disk_shape = collision::translated<disk_shape>;

struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction
    : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, diskfdisk_r> {};
struct topo : simul::topology<disk, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<disk, collision::empty_shape>;
using pointl =
    collision::point<storage::pattern::empty_item, collision::shape::segment>;
struct neighborhood : collision::neighborhood<pointd, pointl> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, space_filter, segment_shape, disk_shape,
          storage::with_properties<
              storage::wrapped<disk, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<translated_disk_shape,
                               some::unbounded_collection>,
              storage::attached<disk, storage::pattern::symbol<"shape">,
                                storage::some::item_ref<disk_shape>>,
              storage::time_invariant<storage::pattern::attr_t<disk, "fext">>,
              storage::diagonal<
                  storage::pattern::attr_t<disk, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::pattern::attr_t<osi, "mass_matrix_assembled">>>> {
};
}  // namespace siconos::config

namespace {

// the states of a run: q and velocity of every disk, then the largest number
// of interactions seen during the run
struct run_result {
  std::vector<double> states;
  std::size_t max_interactions = 0;
  std::string dataflow;  // what the stages of the time steps did
};

run_result run_deposit(bool with_dataflow = false)
{
  using namespace siconos;
  using namespace storage;

  config::data_t data;

  constexpr std::size_t n = 20;  // n x n disks
  constexpr double radius = 1.;
  constexpr double spacing = 2.2 * radius;
  constexpr double m = 1.;
  constexpr double g = 9.81;
  constexpr double width = 2. * spacing * n;
  constexpr double slope_height = 0.5 * width;

  auto disk_shape = storage::add<config::disk_shape>(data);
  disk_shape.radius() = radius;

  // columns of disks just above the slope: they slide, pile up and the
  // contacts are created all along the run
  auto slope_y = [&](double x) {
    return x < 0.5 * width ? slope_height * (1. - x / (0.5 * width)) : 0.;
  };
  using disk_handle_t = decltype(storage::add<config::disk>(data));
  std::vector<disk_handle_t> disks;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      auto disk = storage::add<config::disk>(data);
      prop<"id", config::topo>(disk) =
          static_cast<unsigned int>(disks.size() + 1);
      const double x = spacing * (i + 1.);
      disk.q() = {x, slope_y(x) + 1.1 * radius + spacing * j, 0};
      disk.velocity() = {0, 0, 0};
      disk.mass_matrix().diagonal() << m, m, 0.5 * m * radius * radius;
      prop<"shape">(disk) = disk_shape;
      disk.fext() = {0., -m * g, 0.};
      disks.push_back(disk);
    }
  }

  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = 0.;
  nslaw.mu() = 0.5;

  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;

  auto simulation = storage::add<config::simulation>(data);
  simulation.one_step_integrator().theta() = 0.5;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = 0.;
  simulation.time_discretization().h() = 0.005;
  simulation.time_discretization().tmax() = 1.;

  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;

  // an iterative solver stopped at a loose tolerance: its result depends on
  // the order of the contacts
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  so.set_iparam(SICONOS_IPARAM_MAX_ITER, 100);
  so.set_dparam(SICONOS_DPARAM_TOL, 1e-3);
  osnspb.options() = so;

  auto ngbh = storage::add<config::neighborhood>(data);
  ngbh.create(2.1 * radius);

  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = storage::add<config::diskdisk_r>(data);
  spacef.insert_nonsmooth_law(nslaw, 0, 0);

  // the box: a slope, a floor and two walls, with a point every radius
  auto add_wall = [&](double x1, double y1, double x2, double y2) {
    auto segment = storage::add<config::segment_shape>(data);
    segment.x1() = x1;
    segment.y1() = y1;
    segment.x2() = x2;
    segment.y2() = y2;
    segment.initialize();
    // a whole number of points: the points are iota(0, maxpoints)
    segment.maxpoints() = std::ceil(std::hypot(x2 - x1, y2 - y1) / radius);
    auto wall_r = storage::add<config::diskfsegment_r>(data);
    wall_r.segment() = segment;
    spacef.insert_diskfsegment_r(wall_r);
  };
  add_wall(0., slope_height, 0.5 * width, 0.);
  add_wall(0.5 * width, 0., width, 0.);
  add_wall(0., slope_height, 0., 2. * width);
  add_wall(width, 0., width, 2. * width);

  spacef.make_points();
  ngbh.add_point_sets(0);

  simulation.initialize();
  simulation.enable_dataflow(with_dataflow);

  run_result result;
  while (simulation.has_next_event()) {
    ngbh.update(0);
    ngbh.search();
    spacef.update_index_set0(simulation.current_step());
    simulation.compute_one_step();
    result.max_interactions =
        std::max(result.max_interactions,
                 std::size(storage::handles<config::interaction>(data, 0)));
  }

  result.dataflow = simulation.dataflow_report();

  for (auto& disk : disks) {
    auto& q = storage::attr<"q">(disk, simulation.current_step());
    auto& v = storage::attr<"velocity">(disk, simulation.current_step());
    for (int k = 0; k < 3; ++k) result.states.push_back(q(k));
    for (int k = 0; k < 3; ++k) result.states.push_back(v(k));
  }
  return result;
}

class ReproducibleDisksTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ReproducibleDisksTest);
  CPPUNIT_TEST(testSortedNeighbors);
  CPPUNIT_TEST(testSameRunSameStates);
  CPPUNIT_TEST(testDataflow);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testSortedNeighbors()
  {
    using Real = CompactNSearch::Real;

    // a jittered grid of points and a line of points crossing it
    constexpr unsigned int n = 100;
    std::vector<Real> cloud;
    unsigned int seed = 12345;
    auto jitter = [&seed]() {
      seed = 1664525u * seed + 1013904223u;
      return static_cast<Real>(seed >> 8) / static_cast<Real>(1u << 24) -
             Real(0.5);
    };
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j) {
        cloud.insert(cloud.end(), {i + jitter(), j + jitter(), Real(0)});
      }
    }
    std::vector<Real> line;
    for (unsigned int k = 0; k < 4 * n; ++k) {
      line.insert(line.end(), {Real(0.25) * k, Real(0.2) * k, Real(0)});
    }

    CompactNSearch::NeighborhoodSearch search(2.5);
    search.add_point_set(cloud.data(), cloud.size() / 3);
    search.add_point_set(line.data(), line.size() / 3);
    search.find_neighbors();

    std::size_t pairs = 0;
    for (unsigned int ps = 0; ps < 2; ++ps) {
      const auto& points = search.point_set(ps);
      for (unsigned int other = 0; other < 2; ++other) {
        for (unsigned int i = 0; i < points.n_points(); ++i) {
          for (unsigned int k = 1; k < points.n_neighbors(other, i); ++k) {
            if (points.neighbor(other, i, k - 1) >=
                points.neighbor(other, i, k)) {
              CPPUNIT_FAIL(std::format(
                  "neighbors of point {} of set {} in set {} not sorted", i,
                  ps, other));
            }
          }
          pairs += points.n_neighbors(other, i);
        }
      }
    }
    // many neighbors in many cells, not a trivial search
    CPPUNIT_ASSERT(pairs > 100000);
  }

  void testSameRunSameStates()
  {
    auto first = run_deposit();
    auto second = run_deposit();

    std::print("max interactions: {} / {}\n", first.max_interactions,
               second.max_interactions);

    // a flow with many contacts, not a few isolated ones
    CPPUNIT_ASSERT(first.max_interactions > 200);

    CPPUNIT_ASSERT_EQUAL(first.states.size(), second.states.size());
    for (std::size_t k = 0; k < first.states.size(); ++k) {
      if (first.states[k] != second.states[k]) {
        CPPUNIT_FAIL(
            std::format("runs differ: state {} of disk {}: {:.17g} / {:.17g}",
                        k % 6, k / 6, first.states[k], second.states[k]));
      }
    }
  }

  // what the stages of the time steps do, and recording it changes nothing
  void testDataflow()
  {
    auto recorded = run_deposit(true);
    auto plain = run_deposit();
    std::print("{}", recorded.dataflow);

    CPPUNIT_ASSERT(recorded.states == plain.states);
    CPPUNIT_ASSERT(plain.dataflow.empty());

    const auto& text = recorded.dataflow;
    for (const auto* line :
         {"[device] Compute the free state: v, q at step + 1 without the "
          "contact impulses\n",
          "[host] Compute the active set: the interactions whose gap may "
          "close, the systems "
          "they "
          "involve\n",
          "[host] Solve the frictional contact problem (W, q, mu): lambda\n",
          "    compute_free_velocity + update_position<all<systems>>: all "
          "systems (disk)\n",
          "fusion candidates:\n"}) {
      if (text.find(line) == std::string::npos) {
        CPPUNIT_FAIL(std::string("missing in the dataflow: ") + line);
      }
    }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ReproducibleDisksTest);

}  // namespace
