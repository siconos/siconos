/* Unit tests for polytope x fixed-segment narrowphase (item 4). */

#include <cassert>
#include <cmath>
#include <print>

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
auto check = [](bool ok, auto&& what) {
  if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
};
}  // namespace

int main() {
  using namespace siconos;
  config::data_t data;

  // a square polytope
  auto sh = storage::add<config::polytope_shape>(data);
  sh.nvertices() = 4;
  sh.vertices()[0] = {1, -1};
  sh.vertices()[1] = {1, 1};
  sh.vertices()[2] = {-1, 1};
  sh.vertices()[3] = {-1, -1};

  // a ds with that shape
  auto ds = storage::add<config::polytope_ds>(data);
  storage::prop<"shape", storage::config_declared>(ds) = sh;
  storage::attr<"q">(ds, 0) = Eigen::Vector3d(0, 0, 0);

  // a fixed segment from (0,3) to (0,5) -- separated by 2 along y
  auto seg = storage::add<config::segment_shape>(data);
  seg.set_p1_p2(std::array<double, 6>{0, 3, 0, 0, 5, 0}, 0);
  seg.set_maxpoints(2);
  seg.initialize(0);

  // topology + narrowphase
  auto topo = storage::add<config::topo>(data);
  auto np = storage::add<config::cxd_narrowphase>(data);
  np.topology() = topo;
  // np has no dsds_r set (not needed for polytopefsegment); leave default

  // interaction: topo.link(ds) (self), relation = diskfsegment_r with segment
  auto inter = topo.link(ds);
  auto dsr = storage::add<config::diskfsegment_r>(data);
  dsr.segment() = seg;
  inter.relation() = dsr;
  inter.nslaw() = storage::add<config::nslaw>(data);

  // run the stage
  np.ensure(0);

  // relation should now be rebound to polytopefsegment_r
  using poly_fs_index = storage::index<config::polytopefsegment_r, std::size_t>;
  check(std::holds_alternative<poly_fs_index>(inter.relation()),
        "relation rebound to polytopefsegment_r");

  auto rel = storage::make_handle(data, std::get<poly_fs_index>(inter.relation()));

  // h should be finite (separated: distance from square to segment ~ 2)
  std::print("polytopefsegment: h = {}\n", rel.h());
  check(std::isfinite(rel.h()), "h finite");

  std::print("polytope_fsegment: all checks pass\n");
  return 0;
}
