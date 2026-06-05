#include <fstream>
#include <print>

#include "siconos/siconos.hpp"

namespace siconos::config {
using disk = model::lagrangian_ds;
using nslaw = model::newton_impact_friction;
using diskdisk_r = collision::diskdisk_r;
using diskfdisk_r = collision::diskfdisk_r;
using diskfsegment_r = collision::diskfsegment_r;
using segment_shape = collision::shape::segment;
using disk_shape = collision::shape::disk;
using translated_disk_shape = collision::translated<disk_shape>;

using fc2d = simul::nonsmooth_problem<FrictionContactProblem>;
using osnspb = simul::one_step_nonsmooth_problem<fc2d>;
using solver_options = simul::solver_options;
using interaction =
    simul::interaction<nslaw, diskdisk_r, diskfsegment_r, diskfdisk_r>;
using topo = simul::topology<disk, interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using pointd = collision::point<disk>;
using pointl = collision::point<collision::shape::segment>;
// using pointtds = collision::point<translated_disk_shape>;
using neighborhood = collision::neighborhood<pointd, pointl>;
using space_filter = collision::space_filter<topo, neighborhood>;
using interaction_manager = simul::interaction_manager<space_filter>;
using simulation = simul::time_stepping<td, osi, osnspb, topo>;

using io = io::io<osi>;

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 3>, iparam<"ncgroups", 1>>;
};
}  // namespace siconos::config

int main(int argc, char* argv[])
{
  using namespace siconos;
  using storage::pattern::wrap;
  using namespace storage;

  auto data = storage::make<
      config::env, config::simulation, config::interaction_manager,
      config::io, config::segment_shape, config::disk_shape,
      storage::with_properties<
          storage::wrapped<config::disk, some::unbounded_collection>,
          storage::wrapped<config::diskdisk_r, some::unbounded_collection>,
          storage::wrapped<config::diskfsegment_r,
                           some::unbounded_collection>,
          storage::wrapped<config::pointl, some::unbounded_collection>,
          storage::wrapped<config::pointd, some::unbounded_collection>,
          //          storage::wrapped<config::pointtds,
          //          some::unbounded_collection>,
          storage::wrapped<config::interaction, some::unbounded_collection>,
          storage::wrapped<config::segment_shape, some::unbounded_collection>,
          storage::wrapped<config::disk_shape, some::unbounded_collection>,
          storage::wrapped<config::translated_disk_shape,
                           some::unbounded_collection>,
          storage::attached<config::disk, storage::pattern::symbol<"shape">,
                            storage::some::item_ref<config::disk_shape>>,
          storage::time_invariant<
              storage::pattern::attr_t<config::disk, "fext">>,
          storage::diagonal<
              storage::pattern::attr_t<config::disk, "mass_matrix">>,
          storage::assembled_diagonal<
              storage::pattern::attr_t<typename config::osi::assembled_osi_t,
                                       "mass_matrix_assembled">>>>();

  // unsigned int nDof = 3;         // degrees of freedom for the disk
  double t0 = 0;               // initial computation time
  double tmax = 1;             // final computation time
  double h = 0.005;            // time step
  double position_init = 1.0;  // initial position for lowest bead.
  double velocity_init = 0.0;  // initial velocity for lowest bead.
  double theta = 0.5;          // theta for MoreauJeanOSI integrator
  double radius = 0.1;         // Disk radius
  double m = 1.;               // Disk mass
  double g = 9.81;             // Gravity

  unsigned int ndisks = 3;

  std::print("====> Model loading ...\n");

  // --------------------------
  // -- The dynamical_system --
  // --------------------------

  // one disk shape for all disks
  auto disk_shape = storage::add<config::disk_shape>(data);
  disk_shape.radius() = radius;

  auto disks = storage::add<config::disk>(data, ndisks);

  for (auto [i, disk] : disks | view::enumerate) {
    storage::prop<"id">(disk) = static_cast<uint>(i + 1);

    disk.q() = {0, position_init * static_cast<uint>(i + 1), 0};
    disk.velocity() = {0, velocity_init, 0};
    disk.mass_matrix().diagonal() << m, m, 2. / 5. * m * radius * radius;

    prop<"shape">(disk) = disk_shape;

    // -- Set external forces (weight) --
    disk.fext() = {0., -m * g, 0.};
    //  d2.fext() = {-m * g, 0., 0.};
  }

  // ------------------
  // -- The relation --
  // ------------------

  // -- nslaw --
  double e = 0.9;
  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  auto fc2d = storage::add<config::fc2d>(data);
  fc2d.create();
  fc2d.instance()->dimension = 2;
  fc2d.instance()->mu = 0;

  // ------------------
  // --- Simulation ---
  // ------------------
  auto simulation = storage::add<config::simulation>(data);

  // auto io = storage::add<config::io>(data);

  simulation.one_step_integrator().theta() = theta;
  simulation.one_step_integrator().constraint_activation_threshold() = 0.;
  simulation.time_discretization().t0() = t0;
  simulation.time_discretization().h() = h;

  simulation.time_discretization().tmax() = tmax;
  // -- set the formulation for the one step nonsmooth problem --
  auto osnspb = simulation.one_step_nonsmooth_problem();
  osnspb.problem() = fc2d;

  // -- set the options --
  auto so = storage::add<simul::solver_options>(data);
  so.create(SICONOS_FRICTION_2D_LEMKE);
  osnspb.options() = so;

  auto ngbh = storage::add<config::neighborhood>(data);

  ngbh.create(0.6);  // radius

  auto segment = storage::add<config::segment_shape>(data);
  auto diskdisk_r = storage::add<config::diskdisk_r>(data);
  auto ground_r = storage::add<config::diskfsegment_r>(data);
  ground_r.segment() = segment;
  segment.x1() = -10.;
  segment.y1() = 0.;
  segment.x2() = 10.;
  segment.y2() = 0.;
  segment.initialize();
  segment.maxpoints() = 10;

  auto spacef = storage::add<config::space_filter>(data);
  spacef.neighborhood() = ngbh;
  spacef.diskdisk_r() = diskdisk_r;
  spacef.nslaw() = nslaw;
  spacef.insert_diskfsegment_r(ground_r);
  spacef.make_points();
  ngbh.add_point_sets(0);
  // =========================== End of model definition
  // ===========================
  // ================================= Computation
  // =================================

  //  auto fd = io::open<ascii>("result.dat");

  //  auto disks = storage::handles<config::disk>(data, 0);
  auto disk1 = (disks | view::take(1)).front();
  auto disk2 = (disks | view::take(2)).back();

  // fix this for constant fext
  simulation.initialize();

  std::fstream result_file("result-many.dat");

  // https://stackoverflow.com/questions/72767354/how-to-flush-fmt-output-in-debug-mode
  // std::ofstream cout("result.dat");
  std::print(result_file,
             "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
             simulation.current_step() * simulation.time_step(),
             storage::attr<"q">(disk1, simulation.current_step())(1),
             storage::attr<"q">(disk2, simulation.current_step())(1),
             storage::attr<"velocity">(disk1, simulation.current_step())(1),
             storage::attr<"velocity">(disk2, simulation.current_step())(1),
             0., 0.);

  // auto& vds = storage::prop_values<config::interaction, "vd">(
  //     data, simulation.current_step());

  while (simulation.has_next_event()) {
    ngbh.update(0);
    ngbh.search();
    spacef.update_index_set0(simulation.current_step());

    auto ninvds = simulation.compute_one_step();

    // auto& positions = io.positions(0);

    double p0, lambda;
    if (ninvds > 0) {
      p0 = get_vector(simulation.one_step_integrator().p0_vector_assembled(),
                      0)(0);
      lambda = get_vector(
          simulation.one_step_integrator().lambda_vector_assembled(), 0)(0);
    }
    else {
      p0 = 0;
      lambda = 0;
    }

    std::print(result_file,
               "{:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e}\n",
               simulation.current_step() * simulation.time_step(),
               storage::attr<"q">(disk1, simulation.current_step())(1),
               storage::attr<"q">(disk2, simulation.current_step())(1),
               storage::attr<"velocity">(disk1, simulation.current_step())(1),
               storage::attr<"velocity">(disk2, simulation.current_step())(1),
               p0, lambda);
  }
  //  io::close(fd);
}
