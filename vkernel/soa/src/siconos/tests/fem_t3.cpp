

#include <FENode.hpp>
#include <FemTools.hpp>
#include <FiniteElementModel.hpp>
#include <Material.hpp>
#include <Mesh.hpp>
#include <MeshUtils.hpp>

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "siconos/model/fem.hpp"
#include "siconos/siconos.hpp"
#include "siconos/storage/handle.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using fem = model::finite_element_linear_tids;
using fem_ds = model::rt_lagrangian_ds;
using material = model::material;
using ball = model::lagrangian_ds;
struct lcp : simul::nonsmooth_problem<LinearComplementarityProblem> {};
struct osnspb : simul::one_step_nonsmooth_problem<lcp> {};
using nslaw = model::newton_impact;
struct relation : model::lagrangian_r<nslaw::size> {};
struct rt_relation : model::rt_lagrangian_r {};
struct interaction : simul::interaction<nslaw, relation> {};
struct fem_interaction : simul::rt_interaction<nslaw, rt_relation> {};
struct topo : simul::topology<ball, interaction, fem_ds, fem_interaction> {};
struct osi : simul::one_step_integrator<topo>::moreau_jean {};
struct td : simul::time_discretization<> {};
struct simulation : simul::time_stepping<td, osi, osnspb> {};

using params = map<iparam<"dof", 3>>;

struct make
    : storage::make<
          standard_environment<params>, fem_ds, fem, material, simulation,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<ball, "fext">>,
              storage::diagonal<storage::attr_t<ball, "mass_matrix">>,
              storage::assembled_diagonal<storage::attr_t<
                  typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {
};

}  // namespace siconos::config

int main(int args, char* argv[])
{
  double Ly = 1.0;
  //  std::shared_ptr<Mesh> mesh = create2dMesh2x1();
  //  std::shared_ptr<Mesh> mesh = create2dMeshnxm(50, 15 , 3., Ly);
  // string gmsh_filename = "./mesh_data/triangle_felippa.msh";
  // string gmsh_filename = "./mesh_data/triangle_reference.msh";
  // string gmsh_filename = "./mesh_data/square_6.msh";
  auto gmsh_filename = "./mesh_data/square_200.msh";
  // string gmsh_filename = "./mesh_data/square_2720.msh";

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
  auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(e);
  siconos::algebra::SiconosVector initial_gap{1};
  initial_gap << Ly * 5e-4;
  std::vector<siconos::algebra::SiconosDenseMatrix> Hv;

  namespace storage = siconos::storage;
  namespace config = siconos::config;
  using storage::handle;

  auto data = siconos::config::make();
  handle fe_solid = storage::add<config::fem_ds>(data);
  handle soa_nslaw = storage::add<config::nslaw>(data);
  soa_nslaw.e() = e;

  // -- one step integrator + time discretisation + one step nonsmooth
  // -- problem + simulation
  handle simul = storage::add<config::simulation>(data);

  for (auto node : femodel->nodes()) {
    if (fabs(node->y()) <= 1e-16 and fabs(node->x()) >= 1e-16) {
      std::cout << "contact node number : " << node->num() << " " << node->y()
                << "\n";
      auto idx_y = node->global_dof_index()[1];
      Hv.emplace_back(1, FEsolid->dimension());
      Hv.back().setZero();
      Hv.back()(0, idx_y) = 1.0;

      handle rt_rel = storage::add<config::rt_relation>(data);
      storage::attr<"h_matrix">(rt_rel) = Hv.back();  // Copy H matrix from AOS
      storage::attr<"b">(rt_rel) = initial_gap;        // Copy gap/e vector from AOS

//      handle rt_inter = storage::add<config::fem_interaction>(data);
//      rt_inter.nslaw() = soa_nslaw;
//      rt_inter.relation() = rt_rel;

      auto topo = simul.topology();
      topo.link(fe_solid);
    }
  }



  const auto& massMatrix = FEsolid->mass();
  const auto& stiffnessMatrix = FEsolid->stiffnessMatrix();

  // Print dimensions
  std::cout << "Mass matrix size: " << massMatrix.rows() << " x "
            << massMatrix.cols() << std::endl;
  std::cout << "Stiffness matrix size: " << stiffnessMatrix.rows() << " x "
            << stiffnessMatrix.cols() << std::endl;

  fe_solid.mass_matrix() = massMatrix;
  fe_solid.k_matrix() = stiffnessMatrix;

  // ------------------
  // --- Simulation ---
  // ------------------
  double h = 1e-05;    // time step
  double theta = 1.0;  // theta for MoreauJeanOSI integrator

  int N = ceil((T - t0) / h);  // Number of time steps

  // --- Get the values to be plotted ---
  // -> saved in a matrix dataPlot
  siconos::algebra::Index outputSize = 6;
  siconos::algebra::SiconosDenseMatrix dataPlot(N + 1, outputSize);

  /// OLD AOS CODE
  // SOA/vkernel style computation using handles
  // SOA/vkernel style computation using handles
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

  auto filename = siconos::mechanics::fem::prepareWriteDisplacementforPython(
      "T3_square_200.ref");
  auto mesh = femodel->mesh();
  siconos::mechanics::fem::writeDisplacementforPython(*mesh, *femodel, q,
                                                      filename);
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
      p_val = get_vector(osi.lambda_vector_assembled(), 0)(0);
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
}
