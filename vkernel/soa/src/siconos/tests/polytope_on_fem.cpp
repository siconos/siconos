/* A polytope resting on a deformable body: the reaction must reach the
   mesh nodes.

   The deformable body is a hand-built 4-node elastic block (8 dofs, bottom
   nodes fixed by bc_velocities_0, springs as stiffness) whose top edge
   (0,1)-(1,1) is the contact mesh. A rigid body falls on the middle of that
   edge, once as a disk of radius 0.1 and once as a square of half side
   0.1, in a world with one rigid body type (per-kind shapes, registry).

   Along the normal the square's face contact is the midpoint of the overlap
   window, the point below its centre: it has the disk's gap, G and nodal
   split. So the two runs must agree, rigid body and block alike, and the
   block must deflect. This fails when polytopemesh_r is not assembled on
   the mesh side (one_step_integrator, moreau_jean_element). */

#include <array>
#include <cmath>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include "siconos/siconos.hpp"
#include "siconos/storage/registry.hpp"

// clang-format off
#include "siconos/collision/polytopepolytope_r.hpp"
#include "siconos/collision/polytopemesh_r.hpp"
#include "siconos/collision/shape_of.hpp"
// clang-format on

namespace siconos::config {
namespace some = siconos::storage::some;
using storage::pattern::symbol;

struct body : model::lagrangian_ds {};
struct fem : model::elastic_lagrangian_ds {};
struct nslaw : model::newton_impact_friction {};
using disk_shape = collision::shape::disk;
using polytope_shape = collision::shape::polytope;
using mesh_shape = collision::shape::mesh;
using diskdisk_r = collision::diskdisk_r;
using diskmesh_r = collision::diskmesh_r;
struct topo;
using polytopepolytope_r = collision::polytopepolytope_r<topo>;
using polytopemesh_r = collision::polytopemesh_r<topo>;
using cxd_narrowphase = collision::cxd_narrowphase<topo>;
struct fc2d : simul::nonsmooth_problem<FrictionContactProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<fc2d> {};
struct ct_interaction : simul::interaction<nslaw, diskdisk_r, polytopepolytope_r> {};
struct rt_interaction : simul::rt_ct_interaction<nslaw, diskmesh_r, polytopemesh_r> {};
struct topo : simul::topology<body, ct_interaction, fem, rt_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
using pointd = collision::point<body, collision::empty_shape>;
using pointf = collision::point<fem, mesh_shape>;
struct neighborhood : collision::neighborhood<pointd, pointf> {};
struct space_filter : collision::space_filter<topo, neighborhood> {};
struct interaction_manager : simul::interaction_manager<space_filter> {};
struct simulation : simul::time_stepping<td, osi, osnspb, topo> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};

struct data_t
    : storage::make_runtime_database<
          env, simulation, interaction_manager, body, fem, disk_shape, polytope_shape,
          mesh_shape, diskdisk_r, diskmesh_r, polytopepolytope_r, polytopemesh_r,
          cxd_narrowphase, nslaw, osi, td, osnspb, neighborhood, space_filter,
          storage::with_properties<
              storage::wrapped<body, some::unbounded_collection>,
              storage::wrapped<disk_shape, some::unbounded_collection>,
              storage::wrapped<polytope_shape, some::unbounded_collection>,
              storage::wrapped<mesh_shape, some::unbounded_collection>,
              storage::wrapped<diskdisk_r, some::unbounded_collection>,
              storage::wrapped<diskmesh_r, some::unbounded_collection>,
              storage::wrapped<polytopepolytope_r, some::unbounded_collection>,
              storage::wrapped<polytopemesh_r, some::unbounded_collection>,
              storage::wrapped<pointd, some::unbounded_collection>,
              storage::wrapped<pointf, some::unbounded_collection>,
              storage::wrapped<ct_interaction, some::unbounded_collection>,
              storage::wrapped<rt_interaction, some::unbounded_collection>,
              storage::dynamic_storage<
                  storage::attached<body, symbol<"disk_shape">, some::item_ref<disk_shape>>>,
              storage::dynamic_storage<storage::attached<body, symbol<"polytope_shape">,
                                                         some::item_ref<polytope_shape>>>,
              storage::attached<fem, symbol<"shape">, some::item_ref<mesh_shape>>,
              storage::sparse<storage::pattern::attr_t<fem, "mass_matrix">>,
              storage::sparse<storage::pattern::attr_t<fem, "k_matrix">>,
              storage::unbounded<storage::pattern::attr_t<fem, "q">>,
              storage::unbounded<storage::pattern::attr_t<fem, "velocity">>,
              storage::unbounded<storage::pattern::attr_t<fem, "fext">>,
              storage::time_invariant<storage::pattern::attr_t<body, "fext">>,
              storage::diagonal<storage::pattern::attr_t<body, "mass_matrix">>,
              storage::assembled_diagonal<storage::pattern::attr_t<
                  typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};
}  // namespace siconos::config

namespace {
void check(bool ok, const std::string& what) {
  if (!ok) throw std::runtime_error("CHECK FAILED: " + what);
}

constexpr double r = 0.1;  // disk radius, square half side

struct result {
  std::vector<double> body_y;  // per step
  std::vector<double> top_y;   // per step, y displacement of node 2
  std::vector<double> fem_q;   // final
  std::size_t n_mesh_contacts;
};

result run(bool square) {
  using namespace siconos;
  config::data_t data;

  /* rigid body above the middle of the top edge */
  auto b = storage::add<config::body>(data);
  b.q() = {0.5, 1. + r + 0.02, 0};
  b.velocity() = {0, 0, 0};
  b.fext() = {0., -9.81, 0.};
  if (square) {
    auto psh = storage::add<config::polytope_shape>(data);
    const std::array<std::array<double, 2>, 4> v = {{{-r, -r}, {r, -r}, {r, r}, {-r, r}}};
    for (int i = 0; i < 4; ++i) psh.vertices()[i] = {v[i][0], v[i][1]};
    psh.nvertices() = 4;
    b.mass_matrix().diagonal() << 1., 1., (4 * r * r + 4 * r * r) / 12.;
    collision::set_shape(b, psh);
  } else {
    auto dsh = storage::add<config::disk_shape>(data);
    dsh.radius() = r;
    b.mass_matrix().diagonal() << 1., 1., 0.5 * r * r;
    collision::set_shape(b, dsh);
  }

  /* 4-node block: nodes 0 (0,0), 1 (1,0), 2 (1,1), 3 (0,1); dof 2i = x_i,
     2i+1 = y_i; q is the displacement */
  auto f = storage::add<config::fem>(data);
  const int n = 8;
  f.q() = Eigen::VectorXd::Zero(n);
  f.velocity() = Eigen::VectorXd::Zero(n);
  f.fext() = Eigen::VectorXd::Zero(n);
  Eigen::SparseMatrix<double> M(n, n), K(n, n);
  std::vector<Eigen::Triplet<double>> mt, kt;
  for (int i = 0; i < n; ++i) mt.emplace_back(i, i, 0.25);
  const double k = 1e4;
  for (auto [a, c] :
       std::vector<std::pair<int, int>>{{1, 7}, {3, 5}, {0, 6}, {2, 4}, {6, 4}}) {
    kt.emplace_back(a, a, k);
    kt.emplace_back(c, c, k);
    kt.emplace_back(a, c, -k);
    kt.emplace_back(c, a, -k);
  }
  M.setFromTriplets(mt.begin(), mt.end());
  K.setFromTriplets(kt.begin(), kt.end());
  f.mass_matrix() = M;
  f.k_matrix() = K;
  auto& bc = storage::prop<"bc_velocities_0", config::topo>(f);
  bc.resize(4);
  for (int i = 0; i < 4; ++i) bc[i] = i;  // bottom nodes fixed

  auto msh = storage::add<config::mesh_shape>(data);
  Eigen::MatrixXd nodes(2, 3);
  nodes << 0., 1., 0., 1., 1., 0.;  // node 3 then node 2
  msh.set_nodes(nodes);
  msh.segments().maxpoints() = 11;
  msh.global_indices().resize(4);
  msh.global_indices() << 6, 7, 4, 5;
  storage::prop<"shape", storage::config_declared>(f) = msh;

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
  ngbh.create(0.5);
  auto dsds = storage::add<config::diskdisk_r>(data);
  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = dsds;
  spacef.nslaw() = nslaw;
  spacef.make_points();
  ngbh.add_point_sets(0);
  simulation.initialize();

  result res{};
  std::size_t& n_mesh_contacts = res.n_mesh_contacts;
  while (simulation.has_next_event()) {
    const auto s = simulation.current_step();
    ngbh.update(s);
    ngbh.search();
    spacef.update_index_set0(s);
    n_mesh_contacts = std::max(n_mesh_contacts,
                               std::size(storage::handles<config::rt_interaction>(data, s)));
    simulation.compute_one_step();
    res.body_y.push_back(storage::attr<"q">(b, simulation.current_step())(1));
    res.top_y.push_back(storage::attr<"q">(f, simulation.current_step())(5));
  }
  const auto s = simulation.current_step();
  const auto& fq = storage::attr<"q">(f, s);
  for (int i = 0; i < n; ++i) res.fem_q.push_back(fq(i));
  return res;
}
}  // namespace

int main() {
  const auto disk = run(false);
  const auto square = run(true);
  double max_dy = 0., max_dtop = 0., max_top = 0.;
  for (std::size_t k = 0; k < disk.body_y.size(); ++k) {
    max_dy = std::max(max_dy, std::abs(square.body_y[k] - disk.body_y[k]));
    max_dtop = std::max(max_dtop, std::abs(square.top_y[k] - disk.top_y[k]));
    max_top = std::max(max_top, std::abs(disk.top_y[k]));
  }
  std::print(stderr,
             "disk: final y = {:.9f}, max |top y| = {:.3e}; square: final y = {:.9f}; "
             "max |dy| = {:.3e}, max |dtop| = {:.3e}\n",
             disk.body_y.back(), max_top, square.body_y.back(), max_dy, max_dtop);
  check(disk.n_mesh_contacts > 0 && square.n_mesh_contacts > 0, "both bodies touch the mesh");
  check(square.body_y.size() == disk.body_y.size(), "same number of steps");
  // the block is undamped: the top oscillates around the static deflection
  // (weight / 2 nodes / k ~ 4.9e-4), so test the amplitude, not the sign
  check(max_top > 1e-4, "the disk moves the block");
  check(std::abs(disk.body_y.back() - 1.1) < 1e-2, "the disk rests on the block");
  check(max_dy < 1e-8, "the square follows the disk's trajectory");
  check(max_dtop < 1e-8, "the square moves the block like the disk");
  for (int i = 0; i < 8; ++i)
    check(std::abs(square.fem_q[i] - disk.fem_q[i]) < 1e-8,
          "final block displacement matches (dof " + std::to_string(i) + ")");
  std::print("polytope_on_fem: all checks pass\n");
  return 0;
}
