/* polytopedisk_r: a dynamic disk against a dynamic polytope, in a world with
   one rigid body type (per-kind shapes, registry backend).

   1. the contact outside and with the disk centre inside the polygon.
   2. the normal row of G against central finite differences of h, for both
      bodies, a rotated off-centre square, separated and penetrating.
   3. a disk dropped on a square resting on the ground, through the space
      filter and Moreau-Jean, with the square created before and after the
      disk: the mixed pair is linked disk first, bound to polytopedisk_r, and
      the disk comes to rest on top of the square.
   4. the mixed interaction is removed when the disk leaves and created
      again, on the same relation item, when it comes back. */

#include <array>
#include <cmath>
#include <print>
#include <stdexcept>
#include <string>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
#include "siconos/collision/polytopedisk_r.hpp"
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

struct body : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using disk_shape = collision::shape::disk;
using polytope_shape = collision::shape::polytope;
using segment_shape = collision::shape::segment;
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using polytopedisk_r = collision::polytopedisk_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, polytopepolytope_r,
                                        polytopefsegment_r, polytopedisk_r> {};
struct topo : simul::topology<body, interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<body, collision::empty_shape>;
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
          env, simulation, interaction_manager, body, disk_shape, polytope_shape,
          segment_shape, diskdisk_r, diskfsegment_r, nslaw, interaction, topo,
          polytopepolytope_r, polytopefsegment_r, polytopedisk_r, cxd_narrowphase, osi, td,
          osnspb, neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopefsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopedisk_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>,
              storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                         some::item_ref<polytope_shape>>>,
              storage::time_invariant<storage::pattern::attr_t<body, "fext">>,
              storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>,
              storage::assembled_diagonal<storage::pattern::attr_t<
                  typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {
void check(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("CHECK FAILED: " + what);
}
bool near(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol; }
bool near2(const auto& a, double x, double y, double tol = 1e-9) {
  return near(a[0], x, tol) && near(a[1], y, tol);
}

auto add_disk(siconos::config::data_t& data, double r, double x, double y) {
  using namespace siconos;
  auto sh = storage::add<config::disk_shape>(data);
  sh.radius() = r;
  auto b = storage::add<config::body>(data);
  b.q() = {x, y, 0};
  b.velocity() = {0, 0, 0};
  b.mass_matrix().diagonal() << 1., 1., 0.5 * r * r;
  b.fext() = {0., 0., 0.};
  collision::set_shape(b, sh);
  return b;
}

auto add_square(siconos::config::data_t& data, double a, double x, double y, double theta) {
  using namespace siconos;
  auto sh = storage::add<config::polytope_shape>(data);
  const std::array<std::array<double, 2>, 4> v = {{{-a, -a}, {a, -a}, {a, a}, {-a, a}}};
  for (int i = 0; i < 4; ++i) sh.vertices()[i] = {v[i][0], v[i][1]};
  sh.nvertices() = 4;
  auto b = storage::add<config::body>(data);
  b.q() = {x, y, theta};
  b.velocity() = {0, 0, 0};
  b.mass_matrix().diagonal() << 1., 1., (8 * a * a) / 12.;
  b.fext() = {0., 0., 0.};
  collision::set_shape(b, sh);
  return b;
}

void test_contact() {
  using namespace siconos;
  config::data_t data;
  auto sq = add_square(data, 1, 0, 0, 0);
  auto rel = storage::add<config::polytopedisk_r>(data);

  auto d_out = add_disk(data, 0.5, 3, 0);
  auto out = rel.contact(0, d_out, sq);
  check(out.status == 0 && near(out.h, 1.5), "outside: h = 2 - 0.5");
  check(near2(out.normal, -1, 0) && near2(out.cp1, 2.5, 0) && near2(out.cp2, 1, 0),
        "outside: normal disk -> polytope (-1,0), cp1 (2.5,0), cp2 (1,0)");

  auto d_in = add_disk(data, 0.5, 0.5, 0.2);
  auto in = rel.contact(0, d_in, sq);
  check(in.status == 1 && near(in.h, -1.0), "centre inside: h = -(0.5 + 0.5)");
  check(near2(in.normal, -1, 0) && near2(in.cp1, 0, 0.2) && near2(in.cp2, 1, 0.2),
        "centre inside: normal (-1,0), cp1 = c + r n, cp2 on the nearest edge");
}

void test_jacobian_finite_differences() {
  using namespace siconos;
  for (double dx : {2.4, 0.6}) {  // separated, then centre inside
    config::data_t data;
    auto d = add_disk(data, 0.5, dx, 0.3);
    auto sq = add_square(data, 1, 0.1, -0.2, 0.3);
    auto rel = storage::add<config::polytopedisk_r>(data);
    Eigen::Matrix<double, 2, 3> g1, g2;
    rel.compute_jachq(0, d, sq, g1, g2);
    const double eps = 1e-6;
    for (int body = 0; body < 2; ++body) {
      auto& q = storage::attr<"q">(body == 0 ? d : sq, 0);
      for (int k = 0; k < 3; ++k) {
        const double q0 = q(k);
        q(k) = q0 + eps;
        const double hp = rel.compute_h(0, d, sq);
        q(k) = q0 - eps;
        const double hm = rel.compute_h(0, d, sq);
        q(k) = q0;
        const double gk = body == 0 ? g1(0, k) : g2(0, k);
        check(near((hp - hm) / (2 * eps), gk, 1e-6),
              "normal row of G == dh/dq (centre x " + std::to_string(dx) + ", body " +
                  std::to_string(body + 1) + ", dof " + std::to_string(k) + ")");
      }
    }
  }
}

struct world {
  siconos::config::data_t data;
  decltype(siconos::storage::add<siconos::config::simulation>(data)) simulation;
  decltype(siconos::storage::add<siconos::config::space_filter>(data)) spacef;
  decltype(siconos::storage::add<siconos::config::neighborhood>(data)) ngbh;
  decltype(siconos::storage::add<siconos::config::cxd_narrowphase>(data)) np;

  world()
      : simulation(siconos::storage::add<siconos::config::simulation>(data)),
        spacef(siconos::storage::add<siconos::config::space_filter>(data)),
        ngbh(siconos::storage::add<siconos::config::neighborhood>(data)),
        np(siconos::storage::add<siconos::config::cxd_narrowphase>(data)) {
    using namespace siconos;
    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = 0.;
    nslaw.mu() = 0.5;
    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;
    simulation.one_step_integrator().theta() = 0.5;
    simulation.one_step_integrator().constraint_activation_threshold() = 0.;
    simulation.time_discretization().t0() = 0.;
    simulation.time_discretization().h() = 0.005;
    simulation.time_discretization().tmax() = 1.5;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;
    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_LEMKE);
    osnspb.options() = so;

    auto segment = storage::add<config::segment_shape>(data);
    segment.x1() = -10.;
    segment.y1() = 0.;
    segment.x2() = 10.;
    segment.y2() = 0.;
    segment.initialize();
    segment.maxpoints() = 81;

    ngbh.create(1.5);
    auto dsds = storage::add<config::diskdisk_r>(data);
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = dsds;
    spacef.nslaw() = nslaw;
    np.topology() = spacef.topology();
    np.dsds_r() = dsds;
  }

  void finalize() {
    spacef.make_points();
    ngbh.add_point_sets(0);
    simulation.initialize();
  }

  void filter(std::size_t s) {
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    np.ensure(s);
  }
};

using pd_index = siconos::storage::index<siconos::config::polytopedisk_r, std::size_t>;

/* the polytopedisk_r interactions at step s; each must link the disk first */
std::vector<std::size_t> mixed_relations(world& w, std::size_t s) {
  using namespace siconos;
  std::vector<std::size_t> rels;
  for (auto inter : storage::handles<config::interaction>(w.data, s)) {
    if (!std::holds_alternative<pd_index>(inter.relation())) continue;
    auto ds1 = storage::make_handle(w.data, storage::prop<"ds1", config::topo>(inter));
    auto ds2 = storage::make_handle(w.data, storage::prop<"ds2", config::topo>(inter));
    check(collision::has_shape<collision::shape::disk>(ds1) &&
              collision::has_shape<collision::shape::polytope>(ds2),
          "a mixed pair is linked disk first");
    rels.push_back(std::get<pd_index>(inter.relation()).value());
  }
  return rels;
}

void test_disk_on_square(bool square_first) {
  using namespace siconos;
  world w;
  auto make = [&] {
    if (square_first) {
      auto sq = add_square(w.data, 0.5, 0, 0.5, 0);
      auto d = add_disk(w.data, 0.5, 0, 1.8);
      return std::pair{d, sq};
    }
    auto d = add_disk(w.data, 0.5, 0, 1.8);
    auto sq = add_square(w.data, 0.5, 0, 0.5, 0);
    return std::pair{d, sq};
  };
  auto [d, sq] = make();
  d.fext() = {0., -9.81, 0.};
  sq.fext() = {0., -9.81, 0.};
  w.finalize();

  bool mixed = false;
  while (w.simulation.has_next_event()) {
    const auto s = w.simulation.current_step();
    w.filter(s);
    mixed = mixed || !mixed_relations(w, s).empty();
    w.simulation.compute_one_step();
  }
  const auto s = w.simulation.current_step();
  const auto& qd = storage::attr<"q">(d, s);
  const auto& qs = storage::attr<"q">(sq, s);
  std::print(stderr,
             "disk on square ({}): disk ({:.6f}, {:.6f}), square ({:.6f}, {:.6f}, {:.2e})\n",
             square_first ? "square first" : "disk first", qd(0), qd(1), qs(0), qs(1), qs(2));
  check(mixed, "a polytopedisk_r interaction was created");
  // with friction the split of the tangential reaction between the square's
  // two ground points is not unique (the LCP picks one): the symmetric stack
  // drifts by ~1e-4, it does not slide or roll off
  check(near(qd(0), 0, 1e-3) && near(qs(0), 0, 1e-3) && near(qs(2), 0, 1e-3),
        "symmetric stack: no sliding, no rotation beyond the friction indeterminacy");
  // velocity-level drift of the impacts, as for the disk and square columns
  check(near(qd(1), 1.5, 2e-2) && near(qs(1), 0.5, 2e-2), "the disk rests on the square");
  check(storage::attr<"velocity">(d, s).norm() < 1e-2, "the disk is at rest");
}

void test_mixed_churn() {
  using namespace siconos;
  world w;
  auto sq = add_square(w.data, 0.5, 0, 0.5, 0);
  auto d = add_disk(w.data, 0.5, 0, 1.2);
  w.finalize();
  const std::array<double, 3> y_by_step = {1.2, 10., 1.2};
  const std::array<std::size_t, 3> n_by_step = {1, 0, 1};
  std::size_t rel0 = 0;
  for (std::size_t k = 0; k < y_by_step.size(); ++k) {
    const auto s = w.simulation.current_step();
    storage::attr<"q">(d, s)(1) = y_by_step[k];
    w.filter(s);
    const auto rels = mixed_relations(w, s);
    check(rels.size() == n_by_step[k],
          "mixed interaction count at scripted step " + std::to_string(k));
    if (k == 0) rel0 = rels.front();
    if (k == 2) check(rels.front() == rel0, "the stateless relation item is reused");
    w.simulation.compute_one_step();
  }
  (void)sq;
}
}  // namespace

int main() {
  test_contact();
  test_jacobian_finite_differences();
  test_disk_on_square(true);
  test_disk_on_square(false);
  test_mixed_churn();
  std::print("polytopedisk: all checks pass\n");
  return 0;
}
