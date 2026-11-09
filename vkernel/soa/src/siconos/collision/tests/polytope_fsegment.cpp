/* Unit test for polytopefsegment_r (polytope x fixed segment).

   A side-2 square at the origin and a fixed segment from (0,3) to (0,5):
   the gap is exactly 2 (square top edge y=1 to segment bottom y=3), normal
   (0,1), cp1=(0,1), cp2=(0,3).

   The relation is bound on the self interaction, as
   space_filter::polytopefsegment_activation does; the cxd_narrowphase stage
   handles polytope x polytope pairs only and must leave it untouched. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <cassert>
#include <cmath>
#include <print>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
struct polytope_ds : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using segment_shape = collision::shape::segment;
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, polytopepolytope_r,
                                        polytopefsegment_r> {};
struct topo : simul::topology<polytope_ds, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
struct data_t : storage::make_runtime_database<
                    env, polytope_ds, polytope_shape, segment_shape, diskdisk_r,
                    diskfsegment_r, nslaw, interaction, topo, polytopepolytope_r,
                    polytopefsegment_r, cxd_narrowphase, osi, td, osnspb,
                    storage::with_properties<
                        storage::wrapped<polytope_ds, some::unbounded_collection>,
                        storage::wrapped<interaction, some::unbounded_collection>,
                        storage::attached<polytope_ds, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<polytope_shape>>>> {};
}  // namespace siconos::config

#include <stdexcept>

namespace {
auto near = [](double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; };
}  // namespace

namespace {

class PolytopeFsegmentTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(PolytopeFsegmentTest);
  CPPUNIT_TEST(testPolytopeFsegment);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testPolytopeFsegment() {
    using namespace siconos;
    config::data_t data;

    /* polytope: side-2 square centered at origin, vertices CCW */
    auto sh = storage::add<config::polytope_shape>(data);
    sh.nvertices() = 4;
    sh.vertices()[0] = {1, -1};
    sh.vertices()[1] = {1, 1};
    sh.vertices()[2] = {-1, 1};
    sh.vertices()[3] = {-1, -1};

    auto ds = storage::add<config::polytope_ds>(data);
    storage::prop<"shape", storage::config_declared>(ds) = sh;
    storage::attr<"q">(ds, 0) = Eigen::Vector3d(0, 0, 0);

    /* fixed segment from (0,3) to (0,5) */
    auto seg = storage::add<config::segment_shape>(data);
    seg.set_p1_p2(std::array<double, 6>{0, 3, 0, 0, 5, 0}, 0);
    seg.set_maxpoints(2);
    seg.initialize(0);

    auto topo = storage::add<config::topo>(data);
    auto np = storage::add<config::cxd_narrowphase>(data);
    np.topology() = topo;

    auto inter = topo.link(ds);
    auto rel = storage::add<config::polytopefsegment_r>(data);
    rel.segment() = seg;
    rel.point_index() = 0;
    topo.bind(inter, rel);
    inter.nslaw() = storage::add<config::nslaw>(data);

    /* the stage rebinds polytope x polytope pairs only */
    np.ensure(0);
    using poly_fs_index = storage::index<config::polytopefsegment_r, std::size_t>;
    CPPUNIT_ASSERT_MESSAGE("stage leaves polytopefsegment_r bound",
                           std::holds_alternative<poly_fs_index>(inter.relation()));

    /* h: gap = 2, normal (0,1), cp1=(0,1), cp2=(0,3) */
    double h = rel.compute_h(0, ds);
    std::print("polytopefsegment: h = {}\n", h);
    CPPUNIT_ASSERT_MESSAGE("h == 2 (gap)", near(h, 2));
    const auto ct = rel.contact(0, ds);
    CPPUNIT_ASSERT_MESSAGE("status separated", ct.status == 0);
    CPPUNIT_ASSERT_MESSAGE("normal (0,1)", near(ct.normal[0], 0) && near(ct.normal[1], 1));
    CPPUNIT_ASSERT_MESSAGE("cp1 (0,1)", near(ct.cp1[0], 0) && near(ct.cp1[1], 1));
    CPPUNIT_ASSERT_MESSAGE("cp2 (0,3)", near(ct.cp2[0], 0) && near(ct.cp2[1], 3));

    /* G: nx=0, ny=1, tx=-1, ty=0; arm1 = cp1 - q = (0,1), za1 = (-1,0) */
    Eigen::MatrixXd g1(2, 3);
    rel.compute_jachq(0, ds, g1);
    CPPUNIT_ASSERT_MESSAGE("G1 normal row",
                           near(g1(0, 0), 0) && near(g1(0, 1), -1) && near(g1(0, 2), 0));
    CPPUNIT_ASSERT_MESSAGE("G1 tangent row",
                           near(g1(1, 0), 1) && near(g1(1, 1), 0) && near(g1(1, 2), -1));

    /* the second manifold point: the segment is the incident face, clipped to
       the square's top face, so its far endpoint (0,5) with gap 4 */
    auto rel1 = storage::add<config::polytopefsegment_r>(data);
    rel1.segment() = seg;
    rel1.point_index() = 1;
    const auto ct1 = rel1.contact(0, ds);
    CPPUNIT_ASSERT_MESSAGE("point 1: the far endpoint, gap 4",
                           near(ct1.h, 4) && near(ct1.cp2[1], 5) && near(ct1.cp1[1], 1));

    /* the relation serves every polytope near the segment: a second body
       evaluated through it does not change the first body's contact */
    auto ds2 = storage::add<config::polytope_ds>(data);
    storage::prop<"shape", storage::config_declared>(ds2) = sh;
    storage::attr<"q">(ds2, 0) = Eigen::Vector3d(0, 1.5, 0);
    CPPUNIT_ASSERT_MESSAGE("second body: h == 0.5", near(rel.compute_h(0, ds2), 0.5));
    CPPUNIT_ASSERT_MESSAGE("first body still h == 2", near(rel.compute_h(0, ds), 2));
    CPPUNIT_ASSERT_MESSAGE("first body still cp1 (0,1)", near(rel.contact(0, ds).cp1[1], 1));

    std::print("polytope_fsegment: all checks pass\n");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(PolytopeFsegmentTest);

}  // namespace
