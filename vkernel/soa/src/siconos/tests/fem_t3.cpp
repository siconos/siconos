

#include <FENode.hpp>
#include <FemTools.hpp>
#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <Mesh.hpp>
#include <MeshUtils.hpp>
#include <SiconosKernel.hpp>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "siconos/config/config.hpp"
#include "siconos/config/environment.hpp"
#include "siconos/model/fem.hpp"
#include "siconos/model/lagrangian_ds.hpp"
#include "siconos/model/lagrangian_r.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/interaction.hpp"
#include "siconos/simul/one_step_integrator.hpp"
#include "siconos/simul/one_step_nonsmooth_problem.hpp"
#include "siconos/simul/time_discretization.hpp"
#include "siconos/simul/time_stepping.hpp"
#include "siconos/simul/topology.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
struct fem_ds : model::elastic_lagrangian_ds {};
using ball = model::lagrangian_ds;
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<lcp> {};
using nslaw = model::newton_impact;
struct relation : model::lagrangian_r<nslaw::size> {};
struct rt_relation : model::rt_lagrangian_r {};
struct interaction : simul::interaction<nslaw, relation> {};
struct rt_ct_interaction : simul::rt_ct_interaction<nslaw, rt_relation> {};
struct rt_rt_interaction : simul::rt_rt_interaction<nslaw, rt_relation> {};
struct topo
    : simul::topology<ball, interaction, fem_ds, storage::pattern::empty_item,
                      rt_rt_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

template <typename T>
struct env : standard_environment<T> {
  using params = map<iparam<"dof", 1>>;
};
struct make : storage::make<
                  env, fem_ds, simulation,
                  storage::with_properties<
                      storage::wrapped<config::interaction,
                                       storage::some::unbounded_collection>,
                      storage::wrapped<config::relation,
                                       storage::some::unbounded_collection>,
                      storage::wrapped<config::rt_rt_interaction,
                                       storage::some::unbounded_collection>,
                      storage::wrapped<config::rt_ct_interaction,
                                       storage::some::unbounded_collection>,
                      storage::wrapped<config::rt_relation,
                                       storage::some::unbounded_collection>,
                      storage::unbounded<storage::attr_t<fem_ds, "q">>,
                      storage::unbounded<storage::attr_t<fem_ds, "velocity">>,
                      storage::unbounded<storage::attr_t<fem_ds, "fext">>,
                      storage::sparse<storage::attr_t<fem_ds, "mass_matrix">>,
                      storage::sparse<storage::attr_t<fem_ds, "k_matrix">>>> {
};

}  // namespace siconos::config

int main(int args, char* argv[])
{
  double Ly = 1.0;

  auto gmsh_filename = "square_200.msh";

  // Applied forces
  siconos::algebra::SiconosVector nodal_forces{2};
  nodal_forces << 0., -1e7;
  // Boundary Conditions
  std::vector<int> node_dof_index(2);
  node_dof_index[0] = 0;
  node_dof_index[1] = 1;

  siconos::mechanics::fem::Tags tags;
  tags[siconos::mechanics::fem::MeshTags::bulk_material] = 1;
  tags[siconos::mechanics::fem::MeshTags::boundary_conditions] = 2;
  tags[siconos::mechanics::fem::MeshTags::applied_forces] = 3;

  siconos::mechanics::fem::Material mat{7800, 210e9, 1. / 3};

  auto FEsolid = siconos::mechanics::fem::build_dynamicalsystem_from_gmsh(
      gmsh_filename, tags, mat, nodal_forces, node_dof_index);

  double t0 = 0;     // initial computation time
  double T = 1e-02;  // final computation time
  auto solid =
      std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(t0, T);
  // add the dynamical system in the non smooth dynamical system
  solid->insertDynamicalSystem(FEsolid);
  // Contact Conditions
  auto femodel = FEsolid->FEModel();
  double e = 0.0;
  //  auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
  siconos::algebra::SiconosVector initial_gap{1};
  initial_gap << Ly * 5e-4;
  std::vector<siconos::algebra::SiconosDenseMatrix> Hv;

  namespace storage = siconos::storage;
  namespace config = siconos::config;
  using storage::handle;

  auto data = siconos::config::make();
  handle fe_solid = storage::add<config::fem_ds>(data);
  // fe_solid.dof() = FEsolid->dimension();
  handle nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  handle lcp = storage::add<config::lcp>(data);
  lcp.create(SICONOS_LCP_LEMKE);

  // ------------------
  // --- Simulation ---
  // ------------------
  double h = 1e-05;    // time step
  double theta = 1.0;  // theta for MoreauJeanOSI integrator
  // -- one step integrator + time discretisation + one step nonsmooth
  // -- problem + simulation

  handle simul = storage::add<config::simulation>(data);
  simul.one_step_integrator().theta() = theta;
  simul.one_step_integrator().constraint_activation_threshold() = 0.;
  simul.time_discretization().t0() = t0;
  simul.time_discretization().h() = h;
  simul.time_discretization().tmax() = T;

  handle osnspb = simul.one_step_nonsmooth_problem();
  osnspb.problem() = lcp;

  handle so = storage::add<siconos::simul::solver_options>(data);
  so.create(SICONOS_LCP_LEMKE);
  osnspb.options() = so;

  for (auto node : femodel->nodes()) {
    if (fabs(node->y()) <= 1e-16 and fabs(node->x()) >= 1e-16) {
      std::cout << "contact node number : " << node->num() << " " << node->y()
                << "\n";
      auto idx_y = node->global_dof_index()[1];
      Hv.emplace_back(1, FEsolid->dimension());
      Hv.back().setZero();
      Hv.back()(0, idx_y) = 1.0;

      handle rt_rel = storage::add<config::rt_relation>(data);
      rt_rel.h_matrix().resize(1, FEsolid->dimension());
      rt_rel.h_matrix() = Hv.back();
      rt_rel.b() = initial_gap;

      handle inter = simul.topology().link(fe_solid);
      inter.nslaw() = nslaw;
      inter.relation() = rt_rel;
    }
  }

  fe_solid.mass_matrix() = FEsolid->mass();
  fe_solid.k_matrix() = FEsolid->stiffnessMatrix();
  fe_solid.fext() = FEsolid->fext();
  fe_solid.q() = *FEsolid->q();
  storage::attr<"q">(fe_solid, 1) = storage::attr<"q">(fe_solid, 0);

  fe_solid.velocity() = *FEsolid->velocity();
  storage::attr<"velocity">(fe_solid, 1) = storage::attr<"velocity">(fe_solid, 0);

  auto indices = FEsolid->boundaryConditions()->velocityIndices();
  auto& bc_vel = storage::prop<"bc_velocities_0">(fe_solid);
  bc_vel.resize(indices.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    bc_vel[i] = indices[i];
  }

  int N = ceil((T - t0) / h);  // Number of time steps

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  siconos::algebra::Index outputSize = 6;
  siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);

  simul.initialize();

  // Access initial state at step 0
  auto q = storage::attr<"q">(fe_solid, 0);
  auto v = storage::attr<"velocity">(fe_solid, 0);
  auto dim = q.size();  // Get DOF dimension from state vector

  dataPlot(0, 0) = simul.time_discretization().t0();
  dataPlot(0, 1) = q(dim - 1);
  dataPlot(0, 2) = v(dim - 1);
  dataPlot(0, 3) = 0.0;  // Initial impulse
  dataPlot(0, 4) = q(169);
  dataPlot(0, 5) = v(169);

  auto filename =
      siconos::mechanics::fem::prepareWriteDisplacementforPython("T3");
  auto mesh = femodel->mesh();
  siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q,
                                                      filename);

  auto ds_index = storage::prop<"index">(fe_solid);

  // --- Time loop ---
  std::cout << "====> Start computation ... \n";
  int k = 1;
  auto start = std::chrono::system_clock::now();
  while (simul.has_next_event()) {
    // compute_one_step returns the number of involved dynamical systems
    uint ninvds = simul.compute_one_step();

    auto step = simul.current_step();

    // Access current state from SOA storage
    auto q_current = storage::attr<"q">(fe_solid, step);
    auto v_current = storage::attr<"velocity">(fe_solid, step);

    // Get reaction impulse from assembled vectors only if there are active
    // interactions
    double p_val = 0.0;
    if (ninvds > 0) {
      handle osi = simul.one_step_integrator();
      p_val = get_vector(osi.p0_vector_assembled(), ds_index)(0);
    }

    dataPlot(k, 0) = step * simul.time_step();
    dataPlot(k, 1) = q_current(dim - 1);
    dataPlot(k, 2) = v_current(dim - 1);
    dataPlot(k, 3) = p_val;
    dataPlot(k, 4) = q_current(169);
    dataPlot(k, 5) = v_current(169);

    if (k % 1 == 0)
      siconos::mechanics::fem::writeDisplacementforPython(
          *mesh, *femodel, q_current, filename);

    k++;
    siconos::tools::progressBar((double)k / N);
  }

  auto end = std::chrono::system_clock::now();
  auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();
  std::cout << "\nEnd of computation - Number of iterations done: " << k - 1;
  std::cout << "\nComputation time : " << elapsed << " ms\n";

  // --- Output files ---
  std::cout << "====> Output file writing ...\n";
  dataPlot.conservativeResize(k, outputSize);
  siconos::algebra::io::write("fem_t3.dat", dataPlot,
                              siconos::algebra::io::ASCII_OUT,
                              siconos::algebra::io::WriteType::nodim);
  double eps = 1e-11;
  if (siconos::algebra::io::compareRefFile(dataPlot, "T3_square_200.ref",
                                           eps) >= eps)
    return 1;

  return 0;
}
