/* Cost of the relation and index sets on a world with many relation types:
   the configuration of the python module (py/py_head.hpp, registry backend),
   seven relation types between rigid bodies and two between a rigid body and
   a mesh.

   A pile of ncols x nrows rigid bodies, disks and squares in a checkerboard
   (one body type, the shape kind per body), between two wall segments, on a
   deformable floor: a mass-spring strip (bottom nodes fixed by
   bc_velocities_0) whose top edge is the contact mesh. The space filter and
   the cxd narrowphase create disk-disk, square-square, disk-square,
   disk-wall, square-wall, disk-mesh and square-mesh interactions, in the
   order the bodies meet.

   Usage: mixed_relations [ncols nrows nsteps sort]
     nsteps Moreau-Jean steps (NSGS). With sort > 0 the interactions are
     sorted by relation type (topology::sort_interactions) every sort steps,
     after the filter and the narrowphase: the trajectory then differs from
     the unsorted one (NSGS order).

   Output: per-step mean of the filter (neighbor search + update_index_set0),
   of the narrowphase, of the sort and of compute_one_step, the number of
   interactions of each relation type at the end, a hash of the bodies.
   Built with VKERNEL_RELATION_SETS=0 VKERNEL_INDEX_SETS=0 or not, the hash
   must be the same. Run with OMP_NUM_THREADS=1: the neighbor order must be
   reproducible. */

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <print>
#include <vector>

#include "siconos/collision/space_filter.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopefsegment_r.hpp"
#include "siconos/collision/polytopefdisk_r.hpp"
#include "siconos/collision/polytopedisk_r.hpp"
#include "siconos/collision/polytopemesh_r.hpp"
#include "siconos/collision/rigid_body.hpp"
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

// as py/py_head.hpp
struct body : collision::rigid_body {};
struct fem : model::elastic_lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
using diskmesh_r = collision::diskmesh_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
using mesh_shape = collision::shape::mesh;
using translated_disk_shape = collision::translated<disk_shape>;
using polytope_shape = collision::shape::polytope;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopefsegment_r = collision::polytopefsegment_r<topo>;
using polytopefdisk_r = collision::polytopefdisk_r<topo>;
using polytopedisk_r = collision::polytopedisk_r<topo>;
using polytopemesh_r = collision::polytopemesh_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct ct_interaction
    : simul::interaction<nslaw, diskdisk_r, diskfdisk_r, diskfsegment_r, polytopepolytope_r,
                         polytopefsegment_r, polytopefdisk_r, polytopedisk_r> {};
struct rt_interaction : simul::rt_ct_interaction<nslaw, diskmesh_r, polytopemesh_r> {};
struct topo : simul::topology<body, ct_interaction, fem, rt_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<body, collision::empty_shape>;
using pointf = collision::point<fem, mesh_shape>;
using pointl = collision::point<storage::pattern::empty_item, segment_shape>;
using pointtds = collision::point<storage::pattern::empty_item, translated_disk_shape>;
struct neighborhood : collision::neighborhood<pointd, pointf, pointl, pointtds> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, interaction_manager, segment_shape, disk_shape, mesh_shape,
          polytope_shape, cxd_narrowphase,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopefsegment_r, some::unbounded_collection>,
              storage::wrapped<polytopefdisk_r, some::unbounded_collection>,
              storage::wrapped<polytopedisk_r, some::unbounded_collection>,
              storage::wrapped<polytopemesh_r, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskfsegment_r, some::unbounded_collection>,
              storage::wrapped<diskmesh_r, some::unbounded_collection>,
              storage::wrapped<diskfdisk_r, some::unbounded_collection>,
              storage::wrapped<pointl, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointf, some::unbounded_collection>,
              storage::wrapped<pointtds, some::unbounded_collection>,
              storage::wrapped<ct_interaction, some::unbounded_collection>,
              storage::wrapped<rt_interaction, some::unbounded_collection>,
              storage::wrapped<segment_shape, some::unbounded_collection>,
              storage::wrapped<mesh_shape, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<translated_disk_shape, some::unbounded_collection>,
              storage::without_binding<storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>>,
              storage::without_binding<storage::dynamic_storage<storage::attached<
                  body, symbol<"polytope_shape">, some::item_ref<polytope_shape>>>>,
              storage::attached<fem, symbol<"shape">, some::item_ref<mesh_shape>>,
              storage::sparse<storage::pattern::attr_t<fem, "mass_matrix">>,
              storage::sparse<storage::pattern::attr_t<fem, "k_matrix">>,
              storage::unbounded<storage::pattern::attr_t<fem, "q">>,
              storage::unbounded<storage::pattern::attr_t<fem, "velocity">>,
              storage::unbounded<storage::pattern::attr_t<fem, "fext">>,
              storage::time_invariant<storage::pattern::attr_t<body, "fext">>,
              storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>,
              storage::assembled_diagonal<
                  storage::pattern::attr_t<osi, "mass_matrix_assembled">>>> {};
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
  const int ncols = argc > 1 ? std::atoi(argv[1]) : 30;
  const int nrows = argc > 2 ? std::atoi(argv[2]) : 30;
  const int nsteps = argc > 3 ? std::atoi(argv[3]) : 100;
  const int sort_period = argc > 4 ? std::atoi(argv[4]) : 0;
  const double r = 0.5, m = 1., g = 9.81, spacing = 1.05;

  config::data_t data;

  /* the bodies: disks and squares in a checkerboard */
  auto dsh = storage::add<config::disk_shape>(data);
  dsh.radius() = r;
  auto psh = storage::add<config::polytope_shape>(data);
  const std::array<std::array<double, 2>, 4> square = {{{-r, -r}, {r, -r}, {r, r}, {-r, r}}};
  for (int i = 0; i < 4; ++i) psh.vertices()[i] = {square[i][0], square[i][1]};
  psh.nvertices() = 4;

  const double x0 = 0.5, width = x0 + ncols * spacing + 0.5;
  for (int j = 0; j < nrows; ++j) {
    for (int i = 0; i < ncols; ++i) {
      auto b = storage::add<config::body>(data);
      // odd rows shifted, so that the pile does not stand still
      b.q() = {x0 + r + i * spacing + (j % 2) * 0.1, r + 0.05 + j * spacing, 0};
      b.velocity() = {0, 0, 0};
      b.fext() = {0., -m * g, 0.};
      if ((i + j) % 2 == 0) {
        b.mass_matrix().diagonal() << m, m, 0.5 * m * r * r;
        collision::set_shape(b, dsh);
      } else {
        b.mass_matrix().diagonal() << m, m, m * (4 * r * r + 4 * r * r) / 12.;
        collision::set_shape(b, psh);
      }
    }
  }

  /* the floor: nx + 1 bottom nodes (fixed) and nx + 1 top nodes, at y = -1
     and y = 0, x every dx; dof 2n, 2n + 1 of node n; bottom nodes 0..nx, top
     nodes nx + 1..2 nx + 1; springs along the verticals and the top row */
  const double dx = 0.5;
  const int nx = static_cast<int>(std::ceil(width / dx));
  const int nnodes = 2 * (nx + 1), ndofs = 2 * nnodes;
  auto top = [nx](int j) { return nx + 1 + j; };
  auto f = storage::add<config::fem>(data);
  f.q() = Eigen::VectorXd::Zero(ndofs);
  f.velocity() = Eigen::VectorXd::Zero(ndofs);
  f.fext() = Eigen::VectorXd::Zero(ndofs);
  {
    Eigen::SparseMatrix<double> M(ndofs, ndofs), K(ndofs, ndofs);
    std::vector<Eigen::Triplet<double>> mt, kt;
    for (int i = 0; i < ndofs; ++i) mt.emplace_back(i, i, 0.25);
    const double k = 1e5;
    auto spring = [&kt, k](int a, int c) {
      kt.emplace_back(a, a, k);
      kt.emplace_back(c, c, k);
      kt.emplace_back(a, c, -k);
      kt.emplace_back(c, a, -k);
    };
    for (int j = 0; j <= nx; ++j) {
      spring(2 * j, 2 * top(j));          // x, vertical
      spring(2 * j + 1, 2 * top(j) + 1);  // y, vertical
      if (j < nx) {
        spring(2 * top(j), 2 * top(j + 1));          // x, top row
        spring(2 * top(j) + 1, 2 * top(j + 1) + 1);  // y, top row
      }
    }
    M.setFromTriplets(mt.begin(), mt.end());
    K.setFromTriplets(kt.begin(), kt.end());
    f.mass_matrix() = M;
    f.k_matrix() = K;
  }
  auto& bc = storage::prop<"bc_velocities_0", config::topo>(f);
  bc.resize(2 * (nx + 1));
  for (int i = 0; i < 2 * (nx + 1); ++i) bc[i] = i;  // bottom nodes fixed

  // the top edge, segment j from top node j to top node j + 1
  auto msh = storage::add<config::mesh_shape>(data);
  Eigen::MatrixXd nodes(nx + 1, 3);
  for (int j = 0; j <= nx; ++j) nodes.row(j) << j * dx, 0., 0.;
  msh.set_nodes(nodes);
  msh.segments().maxpoints() = static_cast<std::size_t>(nx * dx / 0.25) + 1;
  msh.global_indices().resize(4 * nx);
  for (int j = 0; j < nx; ++j) {
    msh.global_indices()(4 * j) = 2 * top(j);
    msh.global_indices()(4 * j + 1) = 2 * top(j) + 1;
    msh.global_indices()(4 * j + 2) = 2 * top(j + 1);
    msh.global_indices()(4 * j + 3) = 2 * top(j + 1) + 1;
  }
  storage::prop<"shape", storage::config_declared>(f) = msh;

  /* the walls */
  const double height = nrows * spacing + 2.;
  for (double xw : {0., nx * dx}) {
    auto wall = storage::add<config::segment_shape>(data);
    wall.x1() = xw;
    wall.y1() = 0.;
    wall.x2() = xw;
    wall.y2() = height;
    wall.initialize();
    wall.maxpoints() = static_cast<std::size_t>(height / 0.25) + 1;
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

  auto ngbh = storage::add<config::neighborhood>(data);
  ngbh.create(1.5);
  auto dsds = storage::add<config::diskdisk_r>(data);  // generic ds-ds binder
  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = dsds;
  spacef.insert_nonsmooth_law(nslaw, 0, 0);
  auto np = storage::add<config::cxd_narrowphase>(data);
  np.topology() = spacef.topology();
  np.dsds_r() = dsds;
  spacef.make_points();
  ngbh.add_point_sets(0);
  simulation.initialize();

  double t_filter = 0., t_np = 0., t_sort = 0., t_step = 0.;
  int steps = 0;
  while (simulation.has_next_event()) {
    const auto s = simulation.current_step();
    auto t0 = clock_t_::now();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    t_filter += seconds(t0);
    t0 = clock_t_::now();
    np.ensure(s);
    t_np += seconds(t0);
    if (sort_period > 0 && steps % sort_period == 0) {
      t0 = clock_t_::now();
      simulation.topology().sort_interactions(s);
      t_sort += seconds(t0);
    }
    t0 = clock_t_::now();
    simulation.compute_one_step();
    t_step += seconds(t0);
    ++steps;
  }
  const auto s = simulation.current_step();

  // interactions per relation type: the ct ones, then the rt (mesh) ones
  std::array<std::size_t, 9> count{};
  for (auto inter : storage::handles<config::ct_interaction>(data, s)) {
    ++count[inter.relation().index()];
  }
  for (auto inter : storage::handles<config::rt_interaction>(data, s)) {
    ++count[7 + inter.relation().index()];
  }

  std::uint64_t hash = 14695981039346656037ull;
  for (auto b : storage::handles<config::body>(data, s)) {
    const Eigen::Vector3d q = storage::attr<"q">(b, s);
    const Eigen::Vector3d v = storage::attr<"velocity">(b, s);
    fnv(hash, q.data(), 3 * sizeof(double));
    fnv(hash, v.data(), 3 * sizeof(double));
  }
  const Eigen::VectorXd fq = storage::attr<"q">(f, s);
  fnv(hash, fq.data(), fq.size() * sizeof(double));

  std::print(stderr,
             "bodies={} floor_nodes={} steps={} sort={} filter_ms={:.3f} np_ms={:.3f} "
             "sort_ms={:.3f} step_ms={:.3f} dd={} dfd={} dfs={} pp={} pfs={} pfd={} pd={} "
             "dm={} pm={} hash={:016x}\n",
             ncols * nrows, nnodes, steps, sort_period, 1e3 * t_filter / steps,
             1e3 * t_np / steps, 1e3 * t_sort / steps, 1e3 * t_step / steps, count[0],
             count[1], count[2], count[3], count[4], count[5], count[6], count[7], count[8],
             hash);
  return 0;
}
