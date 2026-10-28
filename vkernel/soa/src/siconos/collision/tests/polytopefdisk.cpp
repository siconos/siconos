/* polytopefdisk_r: a polytope against a fixed disk (registry backend).

   1. point_polytope_distance on a side-2 square, both orientations: signed
      distance (negative inside), normal, closest boundary point.
   2. the relation's contact outside and with the disk centre inside the
      polygon (h = -(r + distance to the nearest edge)).
   3. the normal row of G against central finite differences of h, for a
      rotated, off-centre square, separated and penetrating.
   4. a square falling on a fixed disk through the space filter and
      Moreau-Jean: the relation is bound and the square comes to rest on top
      of the disk without rotating. */

#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>

#include <array>
#include <cmath>
#include <print>
#include <stdexcept>
#include <string>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopefdisk_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

struct polytope_ds : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using polytope_shape = collision::shape::polytope;
using disk_shape = collision::shape::disk;
using translated_disk_shape = collision::translated<disk_shape>;
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
struct topo;
using polytopefdisk_r = collision::polytopefdisk_r<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfdisk_r, polytopefdisk_r> {};
struct topo : simul::topology<polytope_ds, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<polytope_ds, collision::empty_shape>;
using pointtds = collision::point<storage::pattern::empty_item, translated_disk_shape>;
struct neighborhood : collision::neighborhood<pointd, pointtds> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, interaction_manager, polytope_ds, polytope_shape, disk_shape,
          translated_disk_shape, diskdisk_r, diskfdisk_r, polytopefdisk_r, nslaw, osi, td,
          osnspb, neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<polytope_ds, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<translated_disk_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfdisk_r, some::unbounded_collection>,
              storage::wrapped<polytopefdisk_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointtds, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::attached<polytope_ds, symbol<"shape">, some::item_ref<polytope_shape>>,
              storage::time_invariant<storage::pattern::attr_t<polytope_ds, "fext">>,
              storage::diagonal<storage::pattern::attr_t<polytope_ds, "mass_matrix">>,
              storage::assembled_diagonal<storage::pattern::attr_t<
                  typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {
bool near(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }
bool near2(const auto& a, double x, double y, double tol = 1e-9) {
  return near(a[0], x, tol) && near(a[1], y, tol);
}

const std::array<std::array<double, 2>, 4> square_ccw = {{{1, -1}, {1, 1}, {-1, 1}, {-1, -1}}};

void test_distance() {
  using siconos::collision::point_polytope_distance;
  using square_t = std::array<std::array<double, 2>, 4>;
  square_t square_cw;
  for (int i = 0; i < 4; ++i) square_cw[i] = square_ccw[3 - i];
  for (const square_t* sq : {&square_ccw, static_cast<const square_t*>(&square_cw)}) {
    const std::string o = sq == &square_ccw ? "CCW" : "CW";
    auto out = point_polytope_distance(std::array<double, 2>{3, 0}, *sq, 4);
    CPPUNIT_ASSERT_MESSAGE(
        o + " outside: dist 2, normal (1,0), closest (1,0)",
        near(out.dist, 2) && near2(out.normal, 1, 0) && near2(out.closest, 1, 0));
    auto corner = point_polytope_distance(std::array<double, 2>{2, 2}, *sq, 4);
    const double c = 1 / std::sqrt(2.);
    CPPUNIT_ASSERT_MESSAGE(o + " outside a corner: dist sqrt(2), normal diagonal",
                           near(corner.dist, std::sqrt(2.)) && near2(corner.normal, c, c) &&
                               near2(corner.closest, 1, 1));
    auto in = point_polytope_distance(std::array<double, 2>{0.5, 0.2}, *sq, 4);
    CPPUNIT_ASSERT_MESSAGE(
        o + " inside: dist -0.5, normal of the nearest edge, closest on it",
        near(in.dist, -0.5) && near2(in.normal, 1, 0) && near2(in.closest, 1, 0.2));
  }
}

auto add_square(siconos::config::data_t& data, double x, double y, double theta) {
  using namespace siconos;
  auto sh = storage::add<config::polytope_shape>(data);
  for (int i = 0; i < 4; ++i) sh.vertices()[i] = {square_ccw[i][0], square_ccw[i][1]};
  sh.nvertices() = 4;
  auto ds = storage::add<config::polytope_ds>(data);
  storage::prop<"shape", storage::config_declared>(ds) = sh;
  ds.q() = {x, y, theta};
  ds.velocity() = {0, 0, 0};
  ds.mass_matrix().diagonal() << 1., 1., 8. / 12.;
  ds.fext() = {0., 0., 0.};
  return ds;
}

auto add_fixed_disk(siconos::config::data_t& data, double x, double y, double r) {
  using namespace siconos;
  auto dsh = storage::add<config::disk_shape>(data);
  dsh.radius() = r;
  dsh.maxpoints() = 64;
  auto tds = storage::add<config::translated_disk_shape>(data);
  tds.translated() = dsh;
  tds.translation() = {x, y, 0.};
  return tds;
}

void test_relation() {
  using namespace siconos;
  config::data_t data;
  auto ds = add_square(data, 0, 0, 0);
  auto rel = storage::add<config::polytopefdisk_r>(data);

  rel.translated_disk_shape() = add_fixed_disk(data, 3, 0, 0.5);
  auto out = rel.contact(0, ds);
  CPPUNIT_ASSERT_MESSAGE("outside: h = 2 - 0.5", out.status == 0 && near(out.h, 1.5));
  CPPUNIT_ASSERT_MESSAGE(
      "outside: normal (1,0), cp1 (1,0), cp2 (2.5,0)",
      near2(out.normal, 1, 0) && near2(out.cp1, 1, 0) && near2(out.cp2, 2.5, 0));

  rel.translated_disk_shape() = add_fixed_disk(data, 0.5, 0.2, 0.5);
  auto in = rel.contact(0, ds);
  CPPUNIT_ASSERT_MESSAGE("centre inside: h = -(0.5 + 0.5)",
                         in.status == 1 && near(in.h, -1.0));
  CPPUNIT_ASSERT_MESSAGE(
      "centre inside: normal (1,0), cp1 (1,0.2), cp2 = c - r n",
      near2(in.normal, 1, 0) && near2(in.cp1, 1, 0.2) && near2(in.cp2, 0, 0.2));
}

void test_jacobian_finite_differences() {
  using namespace siconos;
  for (double dx : {2.4, 0.6}) {  // separated, then centre inside
    config::data_t data;
    const std::array<double, 3> q0 = {0.1, -0.2, 0.3};
    auto ds = add_square(data, q0[0], q0[1], q0[2]);
    auto rel = storage::add<config::polytopefdisk_r>(data);
    rel.translated_disk_shape() = add_fixed_disk(data, dx, 0.3, 0.5);
    Eigen::Matrix<double, 2, 3> g1;
    rel.compute_jachq(0, ds, g1);
    const double eps = 1e-6;
    for (int k = 0; k < 3; ++k) {
      auto& q = storage::attr<"q">(ds, 0);
      q(k) = q0[k] + eps;
      const double hp = rel.compute_h(0, ds);
      q(k) = q0[k] - eps;
      const double hm = rel.compute_h(0, ds);
      q(k) = q0[k];
      CPPUNIT_ASSERT_MESSAGE("normal row of G == dh/dq (centre x " + std::to_string(dx) +
                                 ", dof " + std::to_string(k) + ")",
                             near((hp - hm) / (2 * eps), g1(0, k), 1e-6));
    }
  }
}

void test_fall_on_fixed_disk() {
  using namespace siconos;
  config::data_t data;
  auto ds = add_square(data, 0, 2.3, 0);
  ds.fext() = {0., -9.81, 0.};
  add_fixed_disk(data, 0, 0, 1);

  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = 0.;
  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
  fc2d.instance()->mu = 0;
  auto simulation = storage::add<config::simulation>(data);
  simulation.one_step_integrator().theta() = 0.5;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = 0.;
  simulation.time_discretization().h() = 0.005;
  simulation.time_discretization().tmax() = 1.;
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_LEMKE);
  osnspb.options() = so;

  auto ngbh = storage::add<config::neighborhood>(data);
  ngbh.create(1.5);
  auto dsds = storage::add<config::diskdisk_r>(data);
  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = dsds;
  spacef.insert_nonsmooth_law(nslaw, 0, 0);
  spacef.make_points();
  ngbh.add_point_sets(0);
  simulation.initialize();

  using pfd_t = storage::index<config::polytopefdisk_r, std::size_t>;
  bool bound = false;
  while (simulation.has_next_event()) {
    const auto s = simulation.current_step();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    for (auto inter : storage::handles<config::interaction>(data, s))
      bound = bound || std::holds_alternative<pfd_t>(inter.relation());
    simulation.compute_one_step();
  }
  const auto s = simulation.current_step();
  const auto& q = storage::attr<"q">(ds, s);
  std::print(stderr, "square on fixed disk: q = ({:.6f}, {:.6f}, {:.3e})\n", q(0), q(1), q(2));
  CPPUNIT_ASSERT_MESSAGE("a polytopefdisk_r interaction was created", bound);
  CPPUNIT_ASSERT_MESSAGE("the square neither slides nor rotates",
                         near(q(0), 0, 1e-9) && near(q(2), 0, 1e-9));
  CPPUNIT_ASSERT_MESSAGE("the square rests on top of the disk", near(q(1), 2.0, 2e-2));
  CPPUNIT_ASSERT_MESSAGE("the square is at rest",
                         storage::attr<"velocity">(ds, s).norm() < 1e-2);
}
}  // namespace

namespace {

class PolytopeFDiskTest : public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(PolytopeFDiskTest);
  CPPUNIT_TEST(testPolytopeFDisk);
  CPPUNIT_TEST_SUITE_END();

 public:
  void testPolytopeFDisk() {
    test_distance();
    test_relation();
    test_jacobian_finite_differences();
    test_fall_on_fixed_disk();
    std::print("polytopefdisk: all checks pass\n");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(PolytopeFDiskTest);

}  // namespace
