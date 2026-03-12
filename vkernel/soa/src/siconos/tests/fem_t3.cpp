
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

#include "siconos/model/fem.hpp"
#include "siconos/siconos.hpp"
#include "siconos/utils/print.hpp"

namespace siconos::config {
using fem = model::finite_element_linear_tids;
using fem_ds = model::rt_lagrangian_ds;
using material = model::material;
using ball = model::lagrangian_ds;
using lcp = simul::nonsmooth_problem<LinearComplementarityProblem>;
using osnspb = simul::one_step_nonsmooth_problem<lcp>;
using nslaw = model::newton_impact;
using relation = model::lagrangian_r<nslaw::size>;
using interaction = simul::interaction<nslaw, relation>;
using fem_interaction = simul::rt_interaction<nslaw, relation>;
using topo =
    simul::topology<ball, interaction, fem_ds, fem_interaction>;
using osi = simul::one_step_integrator<topo>::moreau_jean;
using td = simul::time_discretization<>;
using simulation = simul::time_stepping<td, osi, osnspb>;

using params = map<iparam<"dof", 3>>;

struct make
    : storage::make<
          standard_environment<params>, fem_ds, fem, material, simulation,
          storage::with_properties<
              storage::time_invariant<storage::attr_t<ball, "fext">>,
              storage::diagonal<storage::attr_t<ball, "mass_matrix">>,
              storage::assembled_diagonal<
                storage::attr_t<typename osi::assembled_osi_t, "mass_matrix_assembled">>>> {};

}  // namespace siconos::config

int main(int args, char* argv[])
{
  using namespace siconos;
  using namespace mechanics::fem;

  using Matrix = siconos::algebra::SiconosMatrix;
  using Vector = siconos::algebra::SiconosVector;

  auto gmsh_filename = "./data/square_200.msh";

  auto mesh = createMeshFromGMSH2(gmsh_filename);

  writeMeshforPython(mesh);

  double Ly = 1.0;
  int bulk_material_tag = 1;
  int boundary_condition_tag = 2;
  int applied_force_tag = 3;
  double density = 7800.;

  auto data = config::make();

  auto mat1 = storage::add<config::material>(data);
  mat1.instance().reset(new Material(density, 210e9, 1 / 3.));

  std::map<unsigned int, std::shared_ptr<Material>> materials = {
      {bulk_material_tag, mat1.instance()}};

  auto fe_solid = storage::add<config::fem>(data);
  fe_solid.instance().reset(new FiniteElementLinearTIDS(
      mesh, materials));
  auto fe_model = fe_solid.instance()->FEModel();

  auto nodal_forces = std::make_shared<Vector>(2);
  nodal_forces->setZero();
  (*nodal_forces)(1) = -1e7;
  fe_solid.instance()->applyNodalForces(applied_force_tag, nodal_forces);

  auto node_dof_index = std::make_shared<std::vector<int>>(0);
  node_dof_index->push_back(0);
  node_dof_index->push_back(1);

  fe_solid.instance()->applyDirichletBoundaryConditions(
      boundary_condition_tag, node_dof_index);
  //  fe_solid.instance()->boundaryConditions()->display();

  double t0 = 0;       // initial computation time
  double T = 1e-02;    // final computation time
  double h = 1e-05;    // time step
  double theta = 1.0;  // theta for MoreauJeanOSI integrator

  /*------------------------------------------------- Contact Conditions  */
  double e = 0.0;

  auto nslaw = storage::add<config::nslaw>(data);
  nslaw.e() = e;

  // Create the vector with the size only
  auto initial_gap = std::make_shared<Vector>(1);

  // Set the value using the provided value
  (*initial_gap)(0) = Ly * 5e-4;

  auto inter = storage::add<config::fem_interaction>(data);

}
