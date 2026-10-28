/* topology::sort_interactions (storage::permute): a pile of disks on a ground
   segment, the space filter creating disk-disk and disk-segment interactions
   in any order, sorted every 3 steps. After each sort:

   - the interactions are the same: the multiset of (ds1, ds2, relation,
     lambda, y) is unchanged;
   - their relation types are in the order of the relation variant, the
     relative order of the interactions of a type unchanged;
   - with the relation sets, the set of each relation type holds exactly the
     interactions of this type, in order.

   Registry backend, or the compile time database with -DSORT_STATIC_DB. The
   output (interactions and a hash of the disks per step) is compared
   without and with the sets by the oracle test. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <print>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

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
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, diskfdisk_r> {};
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

using properties = storage::with_properties<
    storage::wrapped<disk, some::unbounded_collection>,
    storage::wrapped<diskdisk_r, some::unbounded_collection>,
    storage::wrapped<diskfsegment_r, some::unbounded_collection>,
    storage::wrapped<pointl, some::unbounded_collection>,
    storage::wrapped<pointd, some::unbounded_collection>,
    storage::wrapped<interaction, some::unbounded_collection>,
    storage::wrapped<segment_shape, some::unbounded_collection>,
    storage::wrapped<disk_shape, some::unbounded_collection>,
    storage::wrapped<translated_disk_shape, some::unbounded_collection>,
    storage::attached<disk, storage::pattern::symbol<"shape">,
                      storage::some::item_ref<disk_shape>>,
    storage::time_invariant<storage::pattern::attr_t<disk, "fext">>,
    storage::diagonal<storage::pattern::attr_t<disk, "mass_matrix">>,
    storage::assembled_diagonal<
        storage::pattern::attr_t<typename osi::assembled_osi_t, "mass_matrix_assembled">>>;

#ifdef SORT_STATIC_DB
struct data_t : storage::make_compile_time_database<env, simulation, interaction_manager,
                                                    segment_shape, disk_shape, properties> {};
#else
struct data_t : storage::make_runtime_database<env, simulation, interaction_manager,
                                               segment_shape, disk_shape, properties> {};
#endif
}  // namespace siconos::config

namespace {

void fnv(std::uint64_t& h, const void* p, std::size_t n) {
  const auto* b = static_cast<const unsigned char*>(p);
  for (std::size_t i = 0; i < n; ++i) {
    h ^= b[i];
    h *= 1099511628211ull;
  }
}

// (relation type, ds1, ds2, relation item, lambda, y) of each interaction,
// in the order of the interactions
auto snapshot(auto& data, auto step) {
  using namespace siconos;
  using entry_t = std::tuple<std::size_t, std::size_t, std::size_t, std::size_t, double,
                             double, double, double>;
  std::vector<entry_t> out;
  for (auto inter : storage::handles<config::interaction>(data, step)) {
    const auto& rel = inter.relation();
    const std::size_t item = std::visit([](auto& r) { return r.value(); }, rel);
    out.emplace_back(rel.index(), storage::prop<"ds1", config::topo>(inter).value(),
                     storage::prop<"ds2", config::topo>(inter).value(), item,
                     storage::attr<"lambda">(inter)(0), storage::attr<"lambda">(inter)(1),
                     storage::attr<"y">(inter)(0), storage::attr<"y">(inter)(1));
  }
  return out;
}

void check_sort(auto& data, auto topo, auto step) {
  using namespace siconos;
  const auto before = snapshot(data, step);
  topo.sort_interactions(step);
  const auto after = snapshot(data, step);

  // the same interactions
  auto sorted_before = before;
  auto sorted_after = after;
  std::sort(sorted_before.begin(), sorted_before.end());
  std::sort(sorted_after.begin(), sorted_after.end());
  CPPUNIT_ASSERT_MESSAGE("sort_interactions changed the interactions",
                         sorted_before == sorted_after);

  // by relation type, the relative order within a type kept
  auto expected = before;
  std::stable_sort(expected.begin(), expected.end(),
                   [](auto& a, auto& b) { return std::get<0>(a) < std::get<0>(b); });
  CPPUNIT_ASSERT_MESSAGE("sort_interactions: not a stable sort by relation type",
                         expected == after);

  // the relation sets hold the interactions of their type, in order
  if constexpr (config::topo::relation_sets) {
    std::size_t k = 0;
    storage::mp::for_each(typename config::interaction::relations{}, [&]<typename R>(R) {
      auto& set = simul::relation_set<config::interaction, R, config::topo>(data);
      for (std::size_t j = 0; j < set.dense_.size(); ++j, ++k) {
        CPPUNIT_ASSERT_MESSAGE("relation set: an interaction out of order",
                               set.dense_keys_[j] == k);
        CPPUNIT_ASSERT_MESSAGE("relation set: not the relation of the interaction",
                               set.dense_[j].value() == std::get<3>(after[k]));
      }
    });
    CPPUNIT_ASSERT_MESSAGE("relation sets: not every interaction", k == after.size());
  }
}

}  // namespace

namespace {

class SortInteractionsTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(SortInteractionsTest);
  CPPUNIT_TEST(testSortInteractions);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testSortInteractions() {
    using namespace siconos;
    using namespace storage;

    config::data_t data;

    const double r = 0.1, m = 1., g = 9.81, spacing = 0.2005;
    const int ncols = 6, nrows = 6;

    auto dsh = storage::add<config::disk_shape>(data);
    dsh.radius() = r;
    for (int j = 0; j < nrows; ++j) {
      for (int i = 0; i < ncols; ++i) {
        auto d = storage::add<config::disk>(data);
        // shifted every other row, so that the pile collapses
        d.q() = {i * spacing + (j % 2) * 0.05, r + 0.05 + j * spacing, 0};
        d.velocity() = {0, 0, 0};
        d.mass_matrix().diagonal() << m, m, 0.5 * m * r * r;
        d.fext() = {0., -m * g, 0.};
        prop<"shape">(d) = dsh;
      }
    }

    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = 0.;
    nslaw.mu() = 0.3;
    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

    auto simulation = storage::add<config::simulation>(data);
    simulation.one_step_integrator().theta() = 0.5;
    simulation.one_step_integrator().constraint_activation_threshold() = 0.;
    simulation.time_discretization().t0() = 0.;
    simulation.time_discretization().h() = 0.005;
    simulation.time_discretization().tmax() = 0.3;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;
    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_NSGS);
    osnspb.options() = so;

    auto ngbh = storage::add<config::neighborhood>(data);
    ngbh.create(0.3);
    auto segment = storage::add<config::segment_shape>(data);
    auto dsds = storage::add<config::diskdisk_r>(data);
    auto ground_r = storage::add<config::diskfsegment_r>(data);
    ground_r.segment() = segment;
    segment.x1() = -1.;
    segment.y1() = 0.;
    segment.x2() = 2.;
    segment.y2() = 0.;
    segment.initialize();
    segment.maxpoints() = 31;

    auto spacef = storage::add<config::space_filter>(data);
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = dsds;
    spacef.insert_nonsmooth_law(nslaw, 0, 0);
    spacef.insert_diskfsegment_r(ground_r);
    spacef.make_points();
    ngbh.add_point_sets(0);
    simulation.initialize();

    std::size_t nsorts = 0;
    std::vector<bool> types_met(3, false);
    int steps = 0;
    while (simulation.has_next_event()) {
      const auto s = simulation.current_step();
      ngbh.update(s);
      ngbh.search();
      spacef.update_index_set0(s);
      if (steps % 3 == 0) {
        check_sort(data, simulation.topology(), s);
        ++nsorts;
      }
      simulation.compute_one_step();
      ++steps;

      const auto snap = snapshot(data, simulation.current_step());
      for (auto& e : snap) types_met[std::get<0>(e)] = true;

      std::uint64_t h = 14695981039346656037ull;
      for (auto d : storage::handles<config::disk>(data, simulation.current_step())) {
        const Eigen::Vector3d q = storage::attr<"q">(d, simulation.current_step());
        fnv(h, q.data(), 3 * sizeof(double));
      }
      std::print("sort_interactions: step {} interactions {} hash {:016x}\n", steps,
                 snap.size(), h);
    }
    // both relation types met, else the sort was not tested
    CPPUNIT_ASSERT_MESSAGE("a single relation type: the sort is not tested",
                           std::count(types_met.begin(), types_met.end(), true) >= 2);
    std::print("sort_interactions: {} sorts, ok\n", nsorts);
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(SortInteractionsTest);

}  // namespace
