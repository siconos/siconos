#include <cmath>
#include <iterator>
#include <print>

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

// Disks bigger than the neighborhood search radius. With one point at
// each center, a pair is seen only when r1 + r2 <= search radius, so
// here every contact would be missed and the disks would go through
// each other. make_points puts points on the circle of the big disks.
int main(int argc, char* argv[]) {
  using namespace siconos;
  using namespace storage;

  config::data_t data;

  double t0 = 0, tmax = 0.5, h = 0.005, theta = 0.5, e = 0.9;
  double search_radius = 0.6;
  double big_r = 1.0, small_r = 0.1;
  double big_m = 100., small_m = 1.;

  auto big_shape = storage::add<config::disk_shape>(data);
  big_shape.radius() = big_r;
  auto small_shape = storage::add<config::disk_shape>(data);
  small_shape.radius() = small_r;

  // big disk at rest, a small disk coming from the right, another big
  // disk coming from above, both 0.2 away from contact.
  auto big0 = storage::add<config::disk>(data);
  auto small = storage::add<config::disk>(data);
  auto big1 = storage::add<config::disk>(data);

  big0.q() = {0, 0, 0};
  big0.velocity() = {0, 0, 0};
  small.q() = {big_r + small_r + 0.2, 0, 0};
  small.velocity() = {-1, 0, 0};
  big1.q() = {0, 2 * big_r + 0.2, 0};
  big1.velocity() = {0, -1, 0};

  big0.mass_matrix().diagonal() << big_m, big_m, big_m * big_r * big_r / 2;
  big1.mass_matrix().diagonal() << big_m, big_m, big_m * big_r * big_r / 2;
  small.mass_matrix().diagonal() << small_m, small_m, small_m * small_r * small_r / 2;
  prop<"shape">(big0) = big_shape;
  prop<"shape">(big1) = big_shape;
  prop<"shape">(small) = small_shape;
  for (auto d : {std::ref(big0), std::ref(small), std::ref(big1)}) {
    d.get().fext() = {0., 0., 0.};
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
  ngbh.create(search_radius);

  auto diskdisk_r = storage::add<config::diskdisk_r>(data);

  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = diskdisk_r;
  spacef.insert_nonsmooth_law(nslaw, 0, 0);
  spacef.make_points();
  ngbh.add_point_sets(0);

  // points spaced by search_radius / 2 on each big circle, one point
  // at the center of the small disk
  auto npoints_big =
      static_cast<std::size_t>(std::ceil(2 * std::numbers::pi * big_r / (search_radius / 2)));
  auto npoints = std::ranges::distance(storage::handles<config::pointd>(data, 0));
  if (npoints != static_cast<std::ptrdiff_t>(2 * npoints_big + 1)) {
    std::print("FAIL: {} points, expected {}\n", npoints, 2 * npoints_big + 1);
    return 1;
  }

  simulation.initialize();

  // both pairs are within the search radius: two interactions,
  // small/big0 and big1/big0. The points of a big disk are neighbors
  // of each other: this must not link a disk with itself.
  ngbh.update(0);
  ngbh.search();
  spacef.update_index_set0(simulation.current_step());
  auto interactions = storage::handles<config::interaction>(data, 0);
  auto ninter0 = std::ranges::distance(interactions);
  for (auto inter : interactions) {
    if (prop<"ds1", config::topo>(inter) == prop<"ds2", config::topo>(inter)) {
      std::print("FAIL: a disk is linked with itself\n");
      return 1;
    }
  }
  if (ninter0 != 2) {
    std::print("FAIL: {} interactions, expected 2\n", ninter0);
    return 1;
  }

  auto dist = [](auto& d1, auto& d2) {
    auto dx = d1.q()(0) - d2.q()(0);
    auto dy = d1.q()(1) - d2.q()(1);
    return std::sqrt(dx * dx + dy * dy);
  };

  double min_gap_small = 1e10, min_gap_big = 1e10;
  std::ptrdiff_t max_inter = 0;
  int steps = 0;
  while (simulation.has_next_event()) {
    ngbh.update(0);
    ngbh.search();
    spacef.update_index_set0(simulation.current_step());
    max_inter = std::max<std::ptrdiff_t>(
        max_inter, std::ranges::distance(storage::handles<config::interaction>(data, 0)));
    simulation.compute_one_step();
    min_gap_small = std::min(min_gap_small, dist(big0, small) - (big_r + small_r));
    min_gap_big = std::min(min_gap_big, dist(big0, big1) - 2 * big_r);
    ++steps;
  }

  std::print("{} steps, max interactions {}, min gaps: small/big {}, big/big {}\n", steps,
             max_inter, min_gap_small, min_gap_big);

  // one step at 1 m/s is 0.005: a missed contact gives a penetration of
  // the order of the whole run, 0.3
  double tol = 2 * h;
  if (max_inter < 2 || min_gap_small < -tol || min_gap_big < -tol) {
    std::print("FAIL: contacts with big disks missed\n");
    return 1;
  }
  // after the impacts, both pairs separate
  if (small.velocity()(0) - big0.velocity()(0) <= 0 ||
      big1.velocity()(1) - big0.velocity()(1) <= 0) {
    std::print("FAIL: no bounce\n");
    return 1;
  }

  std::print("ALL TESTS PASSED\n");
  return 0;
}
