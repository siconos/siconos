/* Unit test for polytopemesh_r (item 2: polytope x mesh edge).

   Builds a side-2 square polytope at the origin and a single mesh segment
   from (0,3) to (0,5): the gap is exactly 2 (square top edge y=1 to segment
   bottom y=3), contact normal (0,1), cp1=(0,1), cp2=(0,3).

   The relation carries its own mesh/item_ref, so no topology/simulation is
   needed for this direct unit test. */

#include <cassert>
#include <cmath>
#include <print>
#include <stdexcept>

#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopemesh_r.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
struct polytope_ds : model::lagrangian_ds {};
struct mesh_ds : model::elastic_lagrangian_ds {};
using polytope_shape = collision::shape::polytope;
using mesh_shape = collision::shape::mesh;
struct topo_tag {};
using polytopemesh_r = collision::polytopemesh_r<topo_tag>;
template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
struct data_t : storage::make_runtime_database<
                    env, polytope_ds, polytope_shape, mesh_ds, mesh_shape, polytopemesh_r,
                    storage::with_properties<
                        storage::wrapped<polytope_ds, some::unbounded_collection>,
                        storage::wrapped<mesh_ds, some::unbounded_collection>,
                        storage::attached<polytope_ds, storage::pattern::symbol<"shape">,
                                          storage::some::item_ref<polytope_shape>>>> {};
}  // namespace siconos::config

#include <stdexcept>

namespace {
auto check = [](bool ok, auto&& what) {
  if (!ok) throw std::runtime_error(std::string("CHECK FAILED: ") + what);
};
auto near = [](double a, double b, double tol = 1e-9) { return std::fabs(a - b) <= tol; };
}  // namespace

int main() {
  using namespace siconos;
  config::data_t data;

  /* polytope: side-2 square centered at origin, vertices CCW */
  auto sh = storage::add<config::polytope_shape>(data);
  sh.nvertices() = 4;
  sh.vertices()[0] = {1, -1};
  sh.vertices()[1] = {1, 1};
  sh.vertices()[2] = {-1, 1};
  sh.vertices()[3] = {-1, -1};

  auto pd = storage::add<config::polytope_ds>(data);
  storage::prop<"shape", storage::config_declared>(pd) = sh;
  storage::attr<"q">(pd, 0) = Eigen::Vector3d(0, 0, 0);

  /* mesh: one segment from (0,3) to (0,5) */
  auto msh = storage::add<config::mesh_shape>(data);
  auto seg = msh.segments();
  seg.insert(0.0, 3.0, 0.0);
  seg.insert(0.0, 5.0, 0.0);
  check(near(seg.x1(0), 0) && near(seg.y1(0), 3), "mesh segment p1");
  check(near(seg.x2(0), 0) && near(seg.y2(0), 5), "mesh segment p2");

  auto md = storage::add<config::mesh_ds>(data);

  /* bind polytopemesh_r relation on (mesh, contact_index 0) */
  auto rel = storage::add<config::polytopemesh_r>(data);
  rel.mesh() = msh;
  rel.contact_index() = 0;

  /* h: gap = 2, normal (0,1), cp1=(0,1), cp2=(0,3) */
  double h = rel.compute_h(0, pd, md);
  std::print("polytopemesh: h = {}\n", h);
  check(near(h, 2), "h == 2 (gap)");
  check(near(rel.status(), 0), "status separated");
  check(near(rel.normal()[0], 0) && near(rel.normal()[1], 1), "normal (0,1)");
  check(near(rel.cp1()[0], 0) && near(rel.cp1()[1], 1), "cp1 (0,1)");
  check(near(rel.cp2()[0], 0) && near(rel.cp2()[1], 3), "cp2 (0,3)");

  /* G matrices */
  Eigen::MatrixXd g1(2, 3), g2(2, 4);
  rel.compute_jachq(0, pd, md, g1, g2);

  /* nx=0,ny=1,tx=-1,ty=0; arm1=cp1-q1=(0,1); za1=(-1,0) */
  check(near(g1(0, 0), 0) && near(g1(0, 1), -1) && near(g1(0, 2), 0), "G1 normal row");
  check(near(g1(1, 0), 1) && near(g1(1, 1), 0) && near(g1(1, 2), -1), "G1 tangent row");
  /* t=0 -> g2 normal row: x1=0,y1=1,x2=0,y2=0 ; tangent row: x1=-1,y1=0,x2=0,y2=0 */
  check(near(g2(0, 0), 0) && near(g2(0, 1), 1) && near(g2(0, 2), 0) && near(g2(0, 3), 0),
        "G2 normal row");
  check(near(g2(1, 0), -1) && near(g2(1, 1), 0) && near(g2(1, 2), 0) && near(g2(1, 3), 0),
        "G2 tangent row");

  std::print("polytopemesh: all checks pass\n");
  return 0;
}
