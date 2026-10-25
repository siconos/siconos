/* Cost of per-kind shape sets against a static "shape" declaration, on a
   disks-only world (registry backend).

   Built twice from this source: disks_shape_static (the shape kind is part
   of the ds type) and disks_shape_per_kind (-DPER_KIND: one rigid ds type,
   the disk set of collision/shape_of.hpp). Both worlds hold the same disks,
   so the trajectories must be bit-identical: the run prints an FNV hash of
   the final q and velocity of every disk.

   Usage: disks_shape_<variant> [ncols nrows nsteps nreps]
     a pile of ncols x nrows disks (radius 0.5, spacing 1.001) on a ground
     segment, nsteps Moreau-Jean steps (NSGS), then nreps evaluations of
     compute_h and compute_jachq over every disk-disk interaction (the only
     place where the shape declaration is read on the hot path).

   Output: per-step mean of the filter (neighbor search + update_index_set0)
   and of compute_one_step, the ns per relation evaluation, the hash.
   Run with OMP_NUM_THREADS=1: the neighbor order must be reproducible. */

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <print>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

struct disk : model::lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using disk_shape = collision::shape::disk;
using segment_shape = collision::shape::segment;
using diskdisk_r = collision::diskdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
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

#ifdef PER_KIND
using shape_decl = storage::dynamic_storage<
    storage::attached<disk, symbol<"disk_shape">, some::item_ref<disk_shape>>>;
#else
using shape_decl = storage::attached<disk, symbol<"shape">, some::item_ref<disk_shape>>;
#endif

struct data_t : storage::make_runtime_database<
                    env, simulation, interaction_manager, disk, disk_shape, segment_shape,
                    diskdisk_r, diskfsegment_r, nslaw, interaction, topo, osi, td, osnspb,
                    neighborhood, space_filter,
                    storage::with_properties<
                        storage::wrapped<disk, some::unbounded_collection>,
                        storage::wrapped<disk_shape, some::unbounded_collection>,
                        storage::wrapped<segment_shape, some::unbounded_collection>,
                        storage::wrapped<diskdisk_r, some::unbounded_collection>,
                        storage::wrapped<diskfsegment_r, some::unbounded_collection>,
                        storage::wrapped<pointd, some::unbounded_collection>,
                        storage::wrapped<pointl, some::unbounded_collection>,
                        storage::wrapped<interaction, some::unbounded_collection>, shape_decl,
                        storage::time_invariant<storage::pattern::attr_t<disk, "fext">>,
                        storage::diagonal<storage::pattern::attr_t<disk, "mass_matrix">>,
                        storage::assembled_diagonal<storage::pattern::attr_t<
                            typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {
using clock_t_ = std::chrono::steady_clock;
double seconds(clock_t_::time_point t0) {
  return std::chrono::duration<double>(clock_t_::now() - t0).count();
}
void fnv(std::uint64_t& h, const void* p, std::size_t n) {
  const auto* b = static_cast<const unsigned char*>(p);
  for (std::size_t i = 0; i < n; ++i) {
    h ^= b[i];
    h *= 1099511628211ull;
  }
}
}  // namespace

int main(int argc, char* argv[]) {
  using namespace siconos;
  const int ncols = argc > 1 ? std::atoi(argv[1]) : 50;
  const int nrows = argc > 2 ? std::atoi(argv[2]) : 50;
  const int nsteps = argc > 3 ? std::atoi(argv[3]) : 100;
  const int nreps = argc > 4 ? std::atoi(argv[4]) : 100;
  const double r = 0.5, m = 1., g = 9.81, spacing = 1.001;

  config::data_t data;

  auto dsh = storage::add<config::disk_shape>(data);
  dsh.radius() = r;
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < ncols; ++i) {
      auto d = storage::add<config::disk>(data);
      d.q() = {i * spacing, r + 0.001 + j * spacing, 0};
      d.velocity() = {0, 0, 0};
      d.mass_matrix().diagonal() << m, m, 0.5 * m * r * r;
      d.fext() = {0., -m * g, 0.};
      collision::set_shape(d, dsh);
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
  simulation.time_discretization().tmax() = 0.005 * nsteps;
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_NSGS);
  osnspb.options() = so;

  const double xmax = ncols * spacing;
  auto segment = storage::add<config::segment_shape>(data);
  segment.x1() = -1.;
  segment.y1() = 0.;
  segment.x2() = xmax + 1.;
  segment.y2() = 0.;
  segment.initialize();
  segment.maxpoints() = static_cast<std::size_t>(4 * (xmax + 2.)) + 1;  // every 0.25

  auto ngbh = storage::add<config::neighborhood>(data);
  ngbh.create(1.5);
  auto dsds = storage::add<config::diskdisk_r>(data);
  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = dsds;
  spacef.nslaw() = nslaw;
  spacef.make_points();
  ngbh.add_point_sets(0);
  simulation.initialize();

  double t_filter = 0., t_step = 0.;
  int steps = 0;
  while (simulation.has_next_event()) {
    const auto s = simulation.current_step();
    auto t0 = clock_t_::now();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    t_filter += seconds(t0);
    t0 = clock_t_::now();
    simulation.compute_one_step();
    t_step += seconds(t0);
    ++steps;
  }
  const auto s = simulation.current_step();

  // relation evaluations over every disk-disk interaction
  using dd_t = storage::index<config::diskdisk_r, std::size_t>;
  std::size_t npairs = 0;
  double acc = 0.;
  Eigen::Matrix<double, 2, 3> g1, g2;
  auto t0 = clock_t_::now();
  for (int rep = 0; rep < nreps; ++rep) {
    npairs = 0;
    for (auto inter : storage::handles<config::interaction>(data, s)) {
      if (!std::holds_alternative<dd_t>(inter.relation())) continue;
      auto rel = storage::make_handle(data, std::get<dd_t>(inter.relation()));
      auto ds1 = storage::make_handle(data, storage::prop<"ds1", config::topo>(inter));
      auto ds2 = storage::make_handle(data, storage::prop<"ds2", config::topo>(inter));
      acc += rel.compute_h(s, ds1, ds2);
      ++npairs;
    }
  }
  const double t_h = seconds(t0);
  t0 = clock_t_::now();
  for (int rep = 0; rep < nreps; ++rep) {
    for (auto inter : storage::handles<config::interaction>(data, s)) {
      if (!std::holds_alternative<dd_t>(inter.relation())) continue;
      auto rel = storage::make_handle(data, std::get<dd_t>(inter.relation()));
      auto ds1 = storage::make_handle(data, storage::prop<"ds1", config::topo>(inter));
      auto ds2 = storage::make_handle(data, storage::prop<"ds2", config::topo>(inter));
      rel.compute_jachq(s, ds1, ds2, g1, g2);
      acc += g1(0, 2) + g2(1, 2);
    }
  }
  const double t_jac = seconds(t0);

  std::uint64_t hash = 14695981039346656037ull;
  for (auto d : storage::handles<config::disk>(data, s)) {
    const Eigen::Vector3d q = storage::attr<"q">(d, s);
    const Eigen::Vector3d v = storage::attr<"velocity">(d, s);
    fnv(hash, q.data(), 3 * sizeof(double));
    fnv(hash, v.data(), 3 * sizeof(double));
  }

  const double per_eval = npairs * static_cast<double>(nreps);
#ifdef PER_KIND
  const char* variant = "per_kind";
#else
  const char* variant = "static";
#endif
  std::print(stderr,
             "variant={} disks={} steps={} pairs={} filter_ms={:.3f} step_ms={:.3f} "
             "h_ns={:.2f} jachq_ns={:.2f} hash={:016x} acc={:.6e}\n",
             variant, ncols * nrows, steps, npairs, 1e3 * t_filter / steps,
             1e3 * t_step / steps, 1e9 * t_h / per_eval, 1e9 * t_jac / per_eval, hash, acc);
  return 0;
}
