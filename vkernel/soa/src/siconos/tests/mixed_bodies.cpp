/* Disks and polytopes as ONE rigid body type, the shape kind per body in
   per-kind sparse sets (collision/shape_of.hpp), on the registry backend,
   driven by the space filter, the cxd narrowphase stage and Moreau-Jean.

   1. two columns on a ground segment, 10 apart (never neighbors):
      two disks (radius 0.5) at x = -5 and two side-1 squares at x = 5,
      dropped from rest. Every interaction must carry the relation of its
      shape pair (diskdisk_r, diskfsegment_r, polytopepolytope_r,
      polytopefsegment_r), both columns must come to rest stacked (centres
      near y = 0.5 and 1.5) and the two columns must agree to 1e-8.

   2. a disk and a square inside the neighborhood radius: there is no
      disk x polytope relation yet, so the space filter must throw instead
      of binding the radius-based diskdisk_r to a polytope. */

#include <array>
#include <cmath>
#include <print>
#include <stdexcept>
#include <string>

#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
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
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, polytopepolytope_r,
                                        polytopefsegment_r> {};
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
          polytopepolytope_r, polytopefsegment_r, cxd_narrowphase, osi, td, osnspb,
          neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopefsegment_r, some::unbounded_collection>,
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

constexpr double m = 1.;
constexpr double g = 9.81;
constexpr double r = 0.5;  // disk radius, square half side

struct world {
  siconos::config::data_t data;
  decltype(siconos::storage::add<siconos::config::simulation>(data)) simulation;
  decltype(siconos::storage::add<siconos::config::space_filter>(data)) spacef;
  decltype(siconos::storage::add<siconos::config::neighborhood>(data)) ngbh;
  decltype(siconos::storage::add<siconos::config::cxd_narrowphase>(data)) np;
  decltype(siconos::storage::add<siconos::config::disk_shape>(data)) dsh;
  decltype(siconos::storage::add<siconos::config::polytope_shape>(data)) psh;

  world()
      : simulation(siconos::storage::add<siconos::config::simulation>(data)),
        spacef(siconos::storage::add<siconos::config::space_filter>(data)),
        ngbh(siconos::storage::add<siconos::config::neighborhood>(data)),
        np(siconos::storage::add<siconos::config::cxd_narrowphase>(data)),
        dsh(siconos::storage::add<siconos::config::disk_shape>(data)),
        psh(siconos::storage::add<siconos::config::polytope_shape>(data)) {
    using namespace siconos;
    dsh.radius() = r;
    const std::array<std::array<double, 2>, 4> square = {{{-r, -r}, {r, -r}, {r, r}, {-r, r}}};
    for (int i = 0; i < 4; ++i) psh.vertices()[i] = {square[i][0], square[i][1]};
    psh.nvertices() = 4;

    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = 0.;

    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

    simulation.one_step_integrator().theta() = 0.5;
    simulation.one_step_integrator().constraint_activation_threshold() = 0.;
    simulation.time_discretization().t0() = 0.;
    simulation.time_discretization().h() = 0.005;
    simulation.time_discretization().tmax() = 2.;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;
    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_LEMKE);
    osnspb.options() = so;

    // ground y = 0, sampled every 0.25 so that a body centre at rest
    // (y = 0.5) always has a sample point inside the radius
    auto segment = storage::add<config::segment_shape>(data);
    segment.x1() = -20.;
    segment.y1() = 0.;
    segment.x2() = 20.;
    segment.y2() = 0.;
    segment.initialize();
    segment.maxpoints() = 161;

    ngbh.create(1.5);
    auto dsds = storage::add<config::diskdisk_r>(data);  // generic ds-ds binder
    spacef.neighborhood() = ngbh;
    spacef.diskdisk_r() = dsds;
    spacef.nslaw() = nslaw;
    np.topology() = spacef.topology();
    np.dsds_r() = dsds;
  }

  auto add_body(auto shape, double x, double y, double inertia) {
    using namespace siconos;
    auto b = storage::add<config::body>(data);
    b.q() = {x, y, 0};
    b.velocity() = {0, 0, 0};
    b.mass_matrix().diagonal() << m, m, inertia;
    b.fext() = {0., -m * g, 0.};
    collision::set_shape(b, shape);
    return b;
  }
  auto add_disk(double x, double y) { return add_body(dsh, x, y, 0.5 * m * r * r); }
  auto add_square(double x, double y) {
    return add_body(psh, x, y, m * (4 * r * r + 4 * r * r) / 12.);
  }

  void finalize() {
    spacef.make_points();
    ngbh.add_point_sets(0);
    simulation.initialize();
  }

  // the filter, then the stage, then the step (the stage rebinds the
  // polytope pairs the filter just created on the generic binder)
  void step() {
    const auto s = simulation.current_step();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    np.ensure(s);
    simulation.compute_one_step();
  }
};

void test_two_columns() {
  using namespace siconos;
  using collision::shape::disk;
  using collision::shape::polytope;

  world w;
  auto d1 = w.add_disk(-5., 0.7);
  auto d2 = w.add_disk(-5., 2.0);
  auto p1 = w.add_square(5., 0.7);
  auto p2 = w.add_square(5., 2.0);
  w.finalize();

  std::array<std::size_t, 4> seen{};  // dd, df, pp, pf
  while (w.simulation.has_next_event()) {
    w.step();
    const auto s = w.simulation.current_step();
    for (auto inter : storage::handles<config::interaction>(w.data, s)) {
      auto ds1 = storage::make_handle(w.data, storage::prop<"ds1", config::topo>(inter));
      auto ds2 = storage::make_handle(w.data, storage::prop<"ds2", config::topo>(inter));
      const bool self = ds1.index().value() == ds2.index().value();
      const bool poly = collision::has_shape<polytope>(ds1);
      check(poly == collision::has_shape<polytope>(ds2), "no mixed pair");
      auto& rel = inter.relation();
      using dd_t = storage::index<config::diskdisk_r, std::size_t>;
      using df_t = storage::index<config::diskfsegment_r, std::size_t>;
      using pp_t = storage::index<config::polytopepolytope_r, std::size_t>;
      using pf_t = storage::index<config::polytopefsegment_r, std::size_t>;
      if (!poly && !self) check(std::holds_alternative<dd_t>(rel), "disk pair: diskdisk_r");
      if (!poly && self)
        check(std::holds_alternative<df_t>(rel), "disk on ground: diskfsegment_r");
      if (poly && !self)
        check(std::holds_alternative<pp_t>(rel), "polytope pair: polytopepolytope_r");
      if (poly && self)
        check(std::holds_alternative<pf_t>(rel), "polytope on ground: polytopefsegment_r");
      seen[(poly ? 2 : 0) + (self ? 1 : 0)] += 1;
    }
  }
  check(seen[0] > 0 && seen[1] > 0 && seen[2] > 0 && seen[3] > 0,
        "every shape pair produced interactions");

  const auto s = w.simulation.current_step();
  auto y = [&](auto& b) { return storage::attr<"q">(b, s)(1); };
  auto v = [&](auto& b) { return storage::attr<"velocity">(b, s).norm(); };
  std::print(stderr, "disks y = {:.6f} {:.6f}, squares y = {:.6f} {:.6f}\n", y(d1), y(d2),
             y(p1), y(p2));
  // Along the normal a disk and an axis-aligned square of half side r have
  // the same gaps and the same G, so the two columns must evolve alike: this
  // compares the polytope relations with the analytic disk ones.
  check(std::abs(y(d1) - y(p1)) < 1e-8 && std::abs(y(d2) - y(p2)) < 1e-8,
        "the disk and square columns evolve alike");
  // Moreau-Jean closes the gap at the velocity level: after a ~4 m/s impact
  // with h = 0.005 the stack keeps a residual interpenetration of ~1e-2
  for (auto* b : {&d1, &p1})
    check(std::abs(y(*b) - 0.5) < 2e-2, "lower body rests on the ground");
  for (auto* b : {&d2, &p2})
    check(std::abs(y(*b) - 1.5) < 2e-2, "upper body rests on the lower");
  for (auto* b : {&d1, &d2, &p1, &p2}) check(v(*b) < 1e-2, "the bodies are at rest");
  for (auto* b : {&p1, &p2})
    check(std::abs(storage::attr<"q">(*b, s)(2)) < 1e-6, "no torque on the squares");
}

void test_disk_polytope_neighbors_throw() {
  world w;
  w.add_disk(0., 0.7);
  w.add_square(1.2, 0.7);
  w.finalize();
  bool thrown = false;
  try {
    w.step();
  } catch (const std::runtime_error& e) {
    thrown = std::string(e.what()).find("disk x polytope") != std::string::npos;
  }
  check(thrown, "a disk and a polytope as neighbors throw");
}
}  // namespace

int main() {
  test_two_columns();
  test_disk_polytope_neighbors_throw();
  std::print("mixed_bodies: all checks pass\n");
  return 0;
}
