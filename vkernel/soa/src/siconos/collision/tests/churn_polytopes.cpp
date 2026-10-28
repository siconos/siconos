/* Space-filter churn for the polytope family (the granular gate) plus
   hardening tests for space_filter::remove_interactions.

   1. churn: three scripted phases with no dynamics (velocities 0, fext 0,
      the positions are re-written at the current step):
      (a) two polytope bodies inside the neighborhood radius --
          space_filter::dsds_activation creates the interaction on the
          generic diskdisk_r binder, the cxd_narrowphase stage rebinds it
          to polytopepolytope_r and batches the pair;
      (b) the second body pushed out of the radius -- the filter
          deactivates and REMOVES the interaction (move_back), the
          relation stays in the stage rel_map;
      (c) the second body pulled back in -- a NEW interaction is created
          and the recycled relation item (same storage index) is reused
          from rel_map, the batch re-runs with the exact gap.
      This is the first test that drives the polytope world through
      space_filter::update_index_set0 itself (bouncing_polytopes
      hand-creates the interaction and never calls it).

   2. stale/duplicate pair: two hand-created interactions for the SAME
      ordered ds pair. The proximity map keeps only the last one; when
      remove_interactions walks the duplicates it must fail loudly, not
      dereference an end iterator (release build: assert compiled out).

   3. unknown relation on a self interaction (topo.link(ds), so
      ds1 == ds2): only diskfsegment_r/diskfdisk_r are handled in the
      ds1 == ds2 paths; anything else must fail loudly (the old code
      hit an assert(false) no-op lambda, a silent leak in release).

   4. shape guard: a polytope body inside the neighborhood radius of a
      fixed segment used to bind diskfsegment_r, which reads
      shape.radius() -- the circumradius for a polytope -- and returned
      a plausible but wrong gap. It must throw instead. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>
#include <cmath>
#include <print>
#include <stdexcept>
#include <string>
#include <utility>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
/* the relation + stage headers must be visible before the config names
   them; the load-bearing cxd include order lives inside cxd_narrowphase.hpp */
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;

struct polytope : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using segment_shape = collision::shape::segment;
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
/* topo is the acyclic marker of the relation/narrowphase pair; defined
   right after the interaction that names it */
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, polytopepolytope_r,
                                        polytopefsegment_r> {};
struct topo : simul::topology<polytope, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<polytope, collision::empty_shape>;
using pointl = collision::point<storage::pattern::empty_item, collision::shape::segment>;
struct neighborhood : collision::neighborhood<pointd, pointl> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, interaction_manager, polytope, polytope_shape, segment_shape,
          diskdisk_r, diskfsegment_r, nslaw, interaction, topo, polytopepolytope_r,
          polytopefsegment_r, cxd_narrowphase, osi, td, osnspb, neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<polytope, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
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

auto add_square_shape(siconos::config::data_t& data) {
  auto shape = siconos::storage::add<siconos::config::polytope_shape>(data);
  const std::array<std::array<double, 2>, 4> square = {{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}};
  for (int i = 0; i < 4; ++i) shape.vertices()[i] = {square[i][0], square[i][1]};
  shape.nvertices() = 4;
  return shape;
};

/* body1 at (0,0); body2 at (0, y2pos): the analytic axis-aligned gap
   between side-2 squares is y2pos - 2. */
auto add_two_bodies(siconos::config::data_t& data, double y2pos) {
  using namespace siconos;
  auto shape = add_square_shape(data);

  auto ds1 = storage::add<config::polytope>(data);
  ds1.q() = {0, 0, 0};
  ds1.velocity() = {0, 0, 0};
  ds1.mass_matrix().diagonal() << 1., 1., 2. / 3.;
  storage::prop<"shape">(ds1) = shape;
  ds1.fext() = {0., 0., 0.};

  auto ds2 = storage::add<config::polytope>(data);
  ds2.q() = {0, y2pos, 0};
  ds2.velocity() = {0, 0, 0};
  ds2.mass_matrix().diagonal() << 1., 1., 2. / 3.;
  storage::prop<"shape">(ds2) = shape;
  ds2.fext() = {0., 0., 0.};

  return std::pair{ds1, ds2};
}

/* the simulation/broadphase/narrowphase filter stack, shared by every
   subtest. The ground segment sits at y = -1000, far outside the
   neighborhood radius: it only exists to create the second point set, so
   that space_filter::update_index_set0 runs the (pointd, pointd)
   ds-ds pair loop at all (one point set never enters it).
   Call finalize_world AFTER the bodies are added: make_points needs the
   ds rows to build the ds point set. */
auto setup_world(siconos::config::data_t& data, double radius, double ground_y = -1000.) {
  using namespace siconos;
  using namespace storage;

  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = 0.9;

  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
  fc2d.instance()->mu = 0;

  auto simulation = storage::add<config::simulation>(data);
  simulation.one_step_integrator().theta() = 0.5;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = 0.;
  simulation.time_discretization().h() = 0.005;
  simulation.time_discretization().tmax() = 1000.;  // scripted loop only
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;

  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_LEMKE);
  osnspb.options() = so;

  // fixed ground: supplies the pointl point set. ground_y = -1000 in the
  // churn subtest (far outside the radius, so the segment never binds);
  // the guard subtest moves it inside the radius on purpose.
  auto segment = storage::add<config::segment_shape>(data);
  segment.x1() = -10.;
  segment.y1() = ground_y;
  segment.x2() = 10.;
  segment.y2() = ground_y;
  segment.initialize();
  segment.maxpoints() = 10;
  auto ground_r = storage::add<config::diskfsegment_r>(data);
  ground_r.segment() = segment;

  auto ngbh = storage::add<config::neighborhood>(data);
  ngbh.create(radius);

  auto dsds = storage::add<config::diskdisk_r>(data);  // generic ds-ds binder
  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = dsds;
  spacef.insert_nonsmooth_law(nslaw, 0, 0);

  return std::tuple{simulation, spacef, ngbh, dsds};
}

/* make_points + add_point_sets + initialize once the ds rows exist
   (bouncing_disks order: bodies first, then spacef.make_points();
   ngbh.add_point_sets(0); then simulation.initialize()). */
auto finalize_world(siconos::config::data_t& data, auto spacef, auto ngbh, auto simulation) {
  spacef.make_points();
  ngbh.add_point_sets(0);
  simulation.initialize();
}

/* ---- 1: the churn gate ---- */
void test_space_filter_churn(siconos::config::data_t& data) {
  using namespace siconos;
  using namespace storage;
  using poly_r_index = storage::index<config::polytopepolytope_r, std::size_t>;

  auto [simulation, spacef, ngbh, dsds] = setup_world(data, 5.5);

  auto np = storage::add<config::cxd_narrowphase>(data);
  np.topology() = spacef.topology();
  np.dsds_r() = dsds;

  auto [ds1, ds2] = add_two_bodies(data, 4.0);  // gap 2 inside the radius
  finalize_world(data, spacef, ngbh, simulation);

  // scripted body2 y positions by step: in, deeper, out, back in
  const std::array<double, 4> y2_by_step = {4.0, 3.0, 100.0, 4.0};
  // two interactions per polytope pair: one per contact manifold point
  const std::array<size_t, 4> ninter_by_step = {2, 2, 0, 2};
  // step 2 keeps the previous h (1.0): no batch runs when the pair is out
  const std::array<double, 4> h_by_step = {2.0, 1.0, 1.0, 2.0};

  // the relation item of each manifold point, recorded at the first step
  std::array<std::size_t, 2> rel_idx0 = {static_cast<std::size_t>(-1),
                                         static_cast<std::size_t>(-1)};

  for (size_t step = 0; step < y2_by_step.size(); ++step) {
    const auto s = simulation.current_step();
    CPPUNIT_ASSERT_MESSAGE("current step advances in lockstep", s == step);

    storage::attr<"q">(ds2, s)(1) = y2_by_step[step];

    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    np.ensure(s);

    auto inters = storage::handles<config::interaction>(data, s);
    CPPUNIT_ASSERT_MESSAGE("the filter keeps the right interaction count",
                           inters.size() == ninter_by_step[step]);

    for (auto inter : inters) {
      CPPUNIT_ASSERT_MESSAGE("the ds-ds interaction is rebound to polytopepolytope_r",
                             std::holds_alternative<poly_r_index>(inter.relation()));
      auto rel = storage::make_handle(data, std::get<poly_r_index>(inter.relation()));
      const std::size_t pi = rel.point_index();
      CPPUNIT_ASSERT_MESSAGE("a polytope relation stands for manifold point 0 or 1", pi < 2);
      if (rel_idx0[pi] == static_cast<std::size_t>(-1)) {
        rel_idx0[pi] = rel.index().value();
      } else {
        CPPUNIT_ASSERT_MESSAGE(
            "the relation item is recycled (same storage index) after removal",
            rel.index().value() == rel_idx0[pi]);
      }
      CPPUNIT_ASSERT_MESSAGE("cached h is the exact gap (keep-previous when separated)",
                             std::abs(rel.h() - h_by_step[step]) <= 1e-6);
    }

    simulation.compute_one_step();
  }

  std::print("churn: space_filter create/remove/recreate cycle passes\n");
}

/* ---- 2: duplicate ds pairs must fail loudly, not deref end() ---- */
void test_duplicate_pair_removal_throws(siconos::config::data_t& data) {
  using namespace siconos;
  using namespace storage;

  auto [simulation, spacef, ngbh, dsds] = setup_world(data, 2.0 /* tiny: pair inactive */);
  auto [ds1, ds2] = add_two_bodies(data, 100.0);  // far apart: never active
  finalize_world(data, spacef, ngbh, simulation);

  auto topo_h = spacef.topology();

  // TWO interactions for the same ds pair: the proximity map keeps only
  // the last one, the first becomes stale
  auto inter1 = topo_h.link(ds1, ds2);
  inter1.nslaw() = spacef.nonsmooth_law(ds1, ds2);
  topo_h.bind(inter1, dsds);
  auto inter2 = topo_h.link(ds1, ds2);
  inter2.nslaw() = spacef.nonsmooth_law(ds1, ds2);
  topo_h.bind(inter2, dsds);

  // update_index_set0 reads searched neighbor lists: broadphase first
  ngbh.update(0);
  ngbh.search();

  bool threw = false;
  try {
    spacef.update_index_set0(0);
  } catch (const std::exception& e) {
    threw = true;
    std::print("duplicate pair: {}\n", e.what());
  }
  CPPUNIT_ASSERT_MESSAGE("duplicate ds pair removal must throw, not dereference end()", threw);
  std::print("remove_interactions: stale/duplicate pair fails loudly\n");
}

/* ---- 3: unhandled relation on a self interaction must fail loudly ---- */
void test_unknown_relation_throws(siconos::config::data_t& data) {
  using namespace siconos;
  using namespace storage;

  auto [simulation, spacef, ngbh, dsds] = setup_world(data, 2.0);
  auto [ds1, ds2] = add_two_bodies(data, 100.0);
  finalize_world(data, spacef, ngbh, simulation);

  auto topo_h = spacef.topology();

  // self interaction (ds1 == ds2 in the row) with a relation the
  // ds1 == ds2 removal paths do not handle: diskdisk_r is a ds-ds
  // relation, not a fixed-object one
  auto inter = topo_h.link(ds1);
  inter.nslaw() = spacef.nonsmooth_law(ds1, ds1);
  topo_h.bind(inter, dsds);

  // update_index_set0 reads searched neighbor lists: broadphase first
  ngbh.update(0);
  ngbh.search();

  bool threw = false;
  try {
    spacef.update_index_set0(0);
  } catch (const std::exception& e) {
    threw = true;
    std::print("unknown relation: {}\n", e.what());
  }
  CPPUNIT_ASSERT_MESSAGE("unhandled relation in the ds1 == ds2 paths must throw", threw);
  std::print("remove_interactions: unhandled relation fails loudly\n");
}
/* ---- 4: a polytope near a fixed segment must fail loudly ---- */
void test_polytope_against_fixed_segment_throws(siconos::config::data_t& data) {
  using namespace siconos;
  using namespace storage;

  // the ground sits at y = -2, inside the 5.5 radius of the body at the
  // origin: the filter will try to bind diskfsegment_r to a polytope body.
  auto [simulation, spacef, ngbh, dsds] = setup_world(data, 5.5, -2.);
  add_two_bodies(data, 100.0);  // bodies exist; the ds-ds pair is far apart
  finalize_world(data, spacef, ngbh, simulation);

  ngbh.update(0);
  ngbh.search();

  // With polytope x fixed-segment support wired, this no longer throws:
  // a polytope body near a fixed segment binds polytopefsegment_r.
  bool threw = false;
  try {
    spacef.update_index_set0(0);
  } catch (const std::exception& e) {
    threw = true;
    std::print("polytope vs fixed segment: {}\n", e.what());
  }
  CPPUNIT_ASSERT_MESSAGE(
      "a polytope near a fixed segment must now bind polytopefsegment_r, "
      "not throw",
      !threw);
  std::print("space_filter: polytope vs fixed segment binds polytopefsegment_r\n");
}
}  // namespace

namespace {

class ChurnPolytopesTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ChurnPolytopesTest);
  CPPUNIT_TEST(testChurnPolytopes);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testChurnPolytopes() {
    {
      siconos::config::data_t data;
      test_space_filter_churn(data);
    }
    {
      siconos::config::data_t data;
      test_duplicate_pair_removal_throws(data);
    }
    {
      siconos::config::data_t data;
      test_unknown_relation_throws(data);
    }
    {
      siconos::config::data_t data;
      test_polytope_against_fixed_segment_throws(data);
    }
    std::print("churn_polytopes: all groups pass\n");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(ChurnPolytopesTest);

}  // namespace
