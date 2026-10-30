#pragma once

/* The scenes of the NSCD test cases and of their comparisons with the
   Anitescu-Tasora and explicit DEM schemes (collision/tests/nscd_cases.cpp,
   tests/nscd_comparisons.cpp): rigid disks and unit squares on a ground
   segment y = 0, the space filter, the cxd narrowphase stage, the
   selectable integrator (the dem integrator with NSCD_DEM_ONLY), NSGS at a
   tolerance of 1e-12 by default. Fixed disks are obstacles without dofs
   (diskfdisk_r). An incline of angle alpha is the ground with the
   gravity tilted by alpha. */

#include <array>
#include <cmath>
#include <numbers>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
#include "siconos/collision/polytopedisk_r.hpp"
#include "siconos/collision/polytopefdisk_r.hpp"
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
using diskfdisk_r = collision::diskfdisk_r;
using translated_disk_shape = collision::translated<disk_shape>;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using polytopedisk_r = collision::polytopedisk_r<topo>;
using polytopefdisk_r = collision::polytopefdisk_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct interaction
    : simul::interaction<nslaw, diskdisk_r, diskfsegment_r, diskfdisk_r, polytopepolytope_r,
                         polytopefsegment_r, polytopedisk_r, polytopefdisk_r> {};
struct topo : simul::topology<body, interaction> {};
#ifdef NSCD_DEM_ONLY
// the DEM cases on the dem integrator: its own storage (y at one time level)
struct osi : simul::one_step_integrator<topo>::dem {};
#else
struct osi : simul::one_step_integrator<topo>::selectable {};
#endif
struct td : simul::time_discretization<> {};
using pointd = collision::point<body, collision::empty_shape>;
using pointl = collision::point<storage::pattern::empty_item, collision::shape::segment>;
using pointtds = collision::point<storage::pattern::empty_item, translated_disk_shape>;
struct neighborhood : collision::neighborhood<pointd, pointl, pointtds> {};
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
          segment_shape, translated_disk_shape, diskdisk_r, diskfsegment_r, diskfdisk_r, nslaw,
          interaction, topo, polytopepolytope_r, polytopefsegment_r, cxd_narrowphase, osi, td,
          osnspb, neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<diskfdisk_r, some::unbounded_collection>,
              storage::wrapped<translated_disk_shape, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopefsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopedisk_r, some::unbounded_collection>,
              storage::wrapped<polytopefdisk_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
              storage::wrapped<pointtds, some::unbounded_collection>,
              storage::wrapped<interaction, some::unbounded_collection>,
              storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>,
              storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                         some::item_ref<polytope_shape>>>,
              storage::time_invariant<storage::pattern::attr_t<body, "fext">>,
              storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::pattern::attr_t<osi, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace nscd {

constexpr double m = 1.;
constexpr double g = 9.81;
constexpr double r = 0.5;  // disk radius, square half side
constexpr double deg = std::numbers::pi / 180.;

namespace schemes = siconos::simul::selectable_scheme;

struct scene {
  double mu = 0.;
  double h = 0.005;
  double tmax = 1.;
  double alpha = 0.;  // incline angle: the gravity tilted towards +x
  // a contact resting at gap 0 has y and ydot of the order of the rounding
  // errors: with a threshold of 0 its activation is decided by their sign,
  // and a dropped row lets the body fall by g h^2 on the other contact
  double threshold = 1e-10;
  std::size_t scheme = schemes::local_problem;
  double theta = 0.5;
  // the penalty law of the explicit DEM scheme
  double kn = 0.;
  double cn = 0.;
  double kt = 0.;
  // the nonsmooth solver: NSGS for Moreau-Jean, the sweeps (iterations) of
  // Anitescu-Tasora
  double tolerance = 1e-12;
  int max_iterations = 100000;
};

struct world {
  siconos::config::data_t data;
  decltype(siconos::storage::add<siconos::config::simulation>(data)) simulation;
  decltype(siconos::storage::add<siconos::config::space_filter>(data)) spacef;
  decltype(siconos::storage::add<siconos::config::neighborhood>(data)) ngbh;
  decltype(siconos::storage::add<siconos::config::cxd_narrowphase>(data)) np;
  decltype(siconos::storage::add<siconos::config::disk_shape>(data)) dsh;
  decltype(siconos::storage::add<siconos::config::polytope_shape>(data)) psh;
  scene sc;

  explicit world(scene sc_)
      : simulation(siconos::storage::add<siconos::config::simulation>(data)),
        spacef(siconos::storage::add<siconos::config::space_filter>(data)),
        ngbh(siconos::storage::add<siconos::config::neighborhood>(data)),
        np(siconos::storage::add<siconos::config::cxd_narrowphase>(data)),
        dsh(siconos::storage::add<siconos::config::disk_shape>(data)),
        psh(siconos::storage::add<siconos::config::polytope_shape>(data)),
        sc(sc_) {
    using namespace siconos;
    dsh.radius() = r;
    const std::array<std::array<double, 2>, 4> square = {{{-r, -r}, {r, -r}, {r, r}, {-r, r}}};
    for (int i = 0; i < 4; ++i) psh.vertices()[i] = {square[i][0], square[i][1]};
    psh.nvertices() = 4;

    auto nslaw = storage::add<config::nslaw>(data);
    nslaw.e() = 0.;
    nslaw.mu() = sc.mu;

    auto fc2d = storage::add<config::fc2d>(data);
    fc2d.create();
    fc2d.instance()->dimension = 2;
    fc2d.instance()->mu = 0;

    simulation.one_step_integrator().scheme() = sc.scheme;
    simulation.one_step_integrator().theta() = sc.theta;
    simulation.one_step_integrator().normal_stiffness() = sc.kn;
    simulation.one_step_integrator().normal_damping() = sc.cn;
    simulation.one_step_integrator().tangential_stiffness() = sc.kt;
    simulation.one_step_integrator().constraint_activation_threshold() = sc.threshold;
    simulation.time_discretization().t0() = 0.;
    simulation.time_discretization().h() = sc.h;
    simulation.time_discretization().tmax() = sc.tmax;
    auto osnspb = simulation.one_step_nonsmooth_problem();
    osnspb.problem() = fc2d;
    auto so = storage::add<simul::solver_options>(data);
    so.create(SICONOS_FRICTION_2D_NSGS);
    so.set_iparam(SICONOS_IPARAM_MAX_ITER, sc.max_iterations);
    so.set_dparam(SICONOS_DPARAM_TOL, sc.tolerance);
    osnspb.options() = so;

    // ground y = 0, sampled every 0.25 so that a body centre in contact
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
    spacef.insert_nonsmooth_law(nslaw, 0, 0);
    np.topology() = spacef.topology();
    np.dsds_r() = dsds;
  }

  auto add_body(auto shape, double x, double y, double vx, double inertia) {
    using namespace siconos;
    auto b = storage::add<config::body>(data);
    b.q() = {x, y, 0};
    b.velocity() = {vx, 0, 0};
    b.mass_matrix().diagonal() << m, m, inertia;
    b.fext() = {m * g * std::sin(sc.alpha), -m * g * std::cos(sc.alpha), 0.};
    collision::set_shape(b, shape);
    return b;
  }
  auto add_disk(double x, double y) { return add_body(dsh, x, y, 0., 0.5 * m * r * r); }
  // a fixed disk of radius r at (x, y): an obstacle without dofs, its points
  // for the neighborhood sampled on its circle
  void add_fixed_disk(double x, double y) {
    using namespace siconos;
    auto shape = storage::add<config::disk_shape>(data);
    shape.radius() = r;
    shape.maxpoints() = 8;
    auto tds = storage::add<config::translated_disk_shape>(data);
    tds.translated() = shape;
    tds.translation() = {x, y, 0.};
    auto rel = storage::add<config::diskfdisk_r>(data);
    rel.translated_disk_shape() = tds;
    spacef.insert_diskfdisk_r(rel);
  }

  auto add_square(double x, double y, double vx = 0.) {
    return add_body(psh, x, y, vx, m * (4 * r * r + 4 * r * r) / 12.);
  }

  void finalize() {
    spacef.make_points();
    ngbh.add_point_sets(0);
    simulation.initialize();
  }

  void step() {
    const auto s = simulation.current_step();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    np.ensure(s);
    simulation.compute_one_step();
  }
};

}  // namespace nscd
