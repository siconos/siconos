// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SolidLinearTIDSTest.hpp"

#include "DruckerPragerPlasticityNSL.hpp"
#include "FEM-MoreauJeanGOSI.hpp"
#include "FEM-Plasticity.hpp"
#include "FENode.hpp"
#include "FemTools.hpp"
#include "FiniteElementModel.hpp"
#include "Interaction.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SolidLinearTIDS.hpp"
#include "StressLinearTIR.hpp"
#include "TimeDiscretisation.hpp"
#include "TimeStepping.hpp"
#include "VonMisesPlasticityNSL.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(SolidLinearTIDSTest);

void SolidLinearTIDSTest::setUp() {
  mesh = siconos::mechanics::fem::createMeshFromGMSH2("square_6.msh");
  auto sd = siconos::algebra::io::readDenseMatrix("S.dat");
  double tol = 1e-14;
  S = sd.sparseView(tol);
  B = siconos::algebra::io::readSparseMatrix("B.dat");
}

void SolidLinearTIDSTest::tearDown() {}

// Mass, K, C
void SolidLinearTIDSTest::testBuildSolidLinearTIDS1() {
  std::cout << "--> Test: constructor 1." << std::endl;
  std::map<int, const siconos::mechanics::fem::Material> materials = {{1, mat}};
  siconos::mechanics::fem::SolidLinearTIDS solid{mesh, materials};

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", solid.dimension() == 10, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", solid.dimension() == 10, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", solid.stressDimension() == 12,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testBuildSolidLinearTIDS1 : ", solid.FEModel()->elements().size() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildSolidLinearTIDS1 : ", solid.stress_read()(0) == 0.0,
                               true);
  auto norm = (solid.elasticityMatrix() - S).norm();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testS : ", norm < 1.0e-16, true);
  norm = (solid.BMatrix() - B).norm();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testB: ", norm < 1.0e-16, true);
  std::cout << "✅ Basic constructor test ended with success." << std::endl;
}

// S = D^{-1} per element for plane strain and plane stress, and K = B^T S^{-1} B
void SolidLinearTIDSTest::testPlaneStressCompliance() {
  std::cout << "--> Test: plane strain / plane stress compliance." << std::endl;
  using siconos::mechanics::fem::AnalysisType2D;
  const double E = 210e9, nu = 0.3;
  for (auto type : {AnalysisType2D::plane_strain, AnalysisType2D::plane_stress}) {
    siconos::mechanics::fem::Material m{7800, E, nu, type, 2.0};
    std::map<int, const siconos::mechanics::fem::Material> materials = {{1, m}};
    siconos::mechanics::fem::SolidLinearTIDS solid{mesh, materials};

    Eigen::Matrix3d D = Eigen::Matrix3d::Zero();
    if (type == AnalysisType2D::plane_strain) {
      double c = E / ((1 + nu) * (1 - 2 * nu));
      D << c * (1 - nu), c * nu, 0, c * nu, c * (1 - nu), 0, 0, 0, c * (1 - 2 * nu) / 2;
    } else {
      double c = E / (1 - nu * nu);
      D << c, c * nu, 0, c * nu, c, 0, 0, 0, c * (1 - nu) / 2;
    }
    Eigen::MatrixXd S_dense(solid.elasticityMatrix());
    Eigen::Matrix3d Dinv = D.inverse();
    auto nel = static_cast<Eigen::Index>(solid.FEModel()->elements().size());
    double err = 0.0;
    for (Eigen::Index e = 0; e < nel; e++)
      err = std::max(err, (S_dense.block<3, 3>(3 * e, 3 * e) - Dinv).cwiseAbs().maxCoeff());
    CPPUNIT_ASSERT_MESSAGE("S blocks = D^-1", err <= 1e-14 * Dinv.cwiseAbs().maxCoeff());
    // no coupling between elements
    double offdiag = 0.0;
    for (Eigen::Index e = 0; e < nel; e++) S_dense.block<3, 3>(3 * e, 3 * e).setZero();
    offdiag = S_dense.cwiseAbs().maxCoeff();
    CPPUNIT_ASSERT_MESSAGE("S block diagonal", offdiag == 0.0);

    // K (thickness included) = B^T S^{-1} B (B scaled by sqrt(area * thickness))
    siconos::algebra::SiconosSparseMatrix K(solid.dimension(), solid.dimension());
    solid.FEModel()->computeStiffnessMatrix(K, materials);
    Eigen::MatrixXd Bd(solid.BMatrix());
    Eigen::MatrixXd Sinv = Eigen::MatrixXd(solid.elasticityMatrix()).inverse();
    Eigen::MatrixXd K_from_BS = Bd.transpose() * Sinv * Bd;
    double errK = (Eigen::MatrixXd(K) - K_from_BS).cwiseAbs().maxCoeff();
    CPPUNIT_ASSERT_MESSAGE("K = B^T S^-1 B", errK <= 1e-12 * K_from_BS.cwiseAbs().maxCoeff());
  }
  std::cout << "✅ plane strain / plane stress compliance test ended with success."
            << std::endl;
}

// lumped mass = row sums of the consistent mass, total mass = rho * thickness * area
void SolidLinearTIDSTest::testLumpedMass() {
  std::cout << "--> Test: lumped mass." << std::endl;
  using siconos::mechanics::fem::MassType;
  const double rho = 7800, thickness = 2.0;
  siconos::mechanics::fem::Material m{
      rho, 210e9, 0.3, siconos::mechanics::fem::AnalysisType2D::plane_stress, thickness};
  std::map<int, const siconos::mechanics::fem::Material> materials = {{1, m}};
  siconos::mechanics::fem::SolidLinearTIDS solid{mesh, materials, MassType::lumped};
  auto fem = solid.FEModel();
  auto n = solid.dimension();
  siconos::algebra::SiconosSparseMatrix Mc(n, n), Ml(n, n);
  fem->computeMassMatrix(Mc, materials);
  fem->computeMassMatrix(Ml, materials, MassType::lumped);
  Eigen::MatrixXd Mcd(Mc), Mld(Ml);
  Eigen::VectorXd rowsum = Mcd.rowwise().sum();
  double scale = rowsum.cwiseAbs().maxCoeff();
  CPPUNIT_ASSERT_MESSAGE("lumped diagonal = row sums",
                         (Mld.diagonal() - rowsum).cwiseAbs().maxCoeff() <= 1e-14 * scale);
  Mld.diagonal().setZero();
  CPPUNIT_ASSERT_MESSAGE("lumped mass diagonal", Mld.cwiseAbs().maxCoeff() == 0.0);
  // the DS uses the lumped mass
  Eigen::MatrixXd Mds(solid.mass());
  CPPUNIT_ASSERT_MESSAGE("DS lumped mass",
                         (Mds.diagonal() - rowsum).cwiseAbs().maxCoeff() <= 1e-14 * scale);
  // total mass (2 dofs per node): 2 rho thickness area, area from the B scaling sum
  double area = 0.0;
  for (auto& fe : fem->elements()) {
    const auto& nd = fe->nodes();
    area += 0.5 * std::abs((nd[1]->x() - nd[0]->x()) * (nd[2]->y() - nd[0]->y()) -
                           (nd[2]->x() - nd[0]->x()) * (nd[1]->y() - nd[0]->y()));
  }
  double total = rowsum.sum();
  CPPUNIT_ASSERT_MESSAGE("total mass",
                         std::abs(total - 2 * rho * thickness * area) <= 1e-12 * total);
  std::cout << "✅ lumped mass test ended with success." << std::endl;
}

// Plane stress uniaxial tension of square_6.msh (left edge: u_x = 0, right edge: v_x = v)
// with Von Mises plasticity solved by fem::nonsmooth_formulations::Plasticity. With T3
// elements the uniform uniaxial state is exact: sigma_xx = E eps while elastic, then
// sigma_xx = sigma_y in every element, sigma_yy = sigma_xy = 0.
void SolidLinearTIDSTest::testVonMisesTension() {
  std::cout << "--> Test: Von Mises uniaxial tension." << std::endl;
  using namespace siconos::mechanics::fem;
  const double E = 210e3, nu = 0.3, rho = 7.8e-9, sigma_y = 210., v = 100., h = 1e-6;
  const int n_steps = 20;  // yield at eps = sigma_y / E = 1e-3, i.e. step 10
  std::map<int, const Material> materials = {
      {1, Material(rho, E, nu, AnalysisType2D::plane_stress)}};
  auto solid = std::make_shared<SolidLinearTIDS>(mesh, materials, MassType::lumped);
  solid->applyDirichletBoundaryConditions(2, {0});
  solid->applyDirichletBoundaryConditions(3, {0}, v);
  auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., n_steps * h);
  nsds->insertDynamicalSystem(solid);
  auto nslaw = std::make_shared<siconos::modeling::VonMisesPlasticityNSL>(sigma_y);
  auto nel = solid->FEModel()->elements().size();
  std::vector<siconos::algebra::SiconosMatrix> C(nel);
  for (size_t e = 0; e < nel; e++) {
    C[e] = solid->elementStressSelection(e);
    auto inter = std::make_shared<siconos::modeling::Interaction>(
        nslaw, std::make_shared<StressLinearTIR>(C[e]));
    nsds->link(inter, solid);
  }
  auto osi = std::make_shared<integrators::MoreauJeanGOSI>(1.0);
  auto td = std::make_shared<siconos::simulation::TimeDiscretisation>(0., h);
  auto osnspb = std::make_shared<nonsmooth_formulations::Plasticity>();
  osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-10;
  osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 100000;
  auto s = std::make_shared<siconos::simulation::TimeStepping>(nsds, td, osi, osnspb);

  // The prescribed velocity jumps from 0 to v at the first step, and the yield starts at
  // step 10: the inertia gives transients (relative error ~1e-3) that decay in 3 steps.
  auto sigma = solid->stress_read();
  int k = 0;
  double max_err = 0.0, max_err_transient = 0.0, max_vm = 0.0;
  while (s->hasNextEvent()) {
    s->computeOneStep();
    k++;
    double sxx_expected = std::min(E * v * k * h, sigma_y);  // length 1
    double step_err = 0.0;
    for (size_t e = 0; e < nel; e++) {
      Eigen::Vector3d se = sigma.segment<3>(3 * e) / solid->elementStressScaling(e);
      step_err = std::max(step_err, std::abs(se(0) - sxx_expected) / sigma_y);
      step_err = std::max(step_err, std::max(std::abs(se(1)), std::abs(se(2))) / sigma_y);
      double vm = std::sqrt(se(0) * se(0) + se(1) * se(1) - se(0) * se(1) + 3 * se(2) * se(2));
      max_vm = std::max(max_vm, vm / sigma_y);
    }
    bool transient = (k <= 3) || (k >= 11 && k <= 13);
    (transient ? max_err_transient : max_err) =
        std::max(transient ? max_err_transient : max_err, step_err);
    s->nextStep();
  }
  std::cout << "max relative error on the stresses: " << max_err
            << " (transients: " << max_err_transient << "), max VM / sigma_y: " << max_vm
            << std::endl;
  CPPUNIT_ASSERT_MESSAGE("uniaxial Von Mises tension", max_err <= 1e-6);
  CPPUNIT_ASSERT_MESSAGE("uniaxial Von Mises tension, transients", max_err_transient <= 1e-3);
  CPPUNIT_ASSERT_MESSAGE("yield condition", max_vm <= 1.0 + 1e-10);
  std::cout << "✅ Von Mises uniaxial tension test ended with success." << std::endl;
}

// DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain: the cone ||sigma_T|| <= eta sigma_N + c'
// of the numerics (sigma_N = -(sxx + syy)/sqrt2, sigma_T = ((sxx - syy)/sqrt2, sqrt2 sxy)) is
// the Mohr-Coulomb criterion R <= -p sin(phi) + c cos(phi) of the plane
void SolidLinearTIDSTest::testDruckerPragerLaw() {
  std::cout << "--> Test: Drucker-Prager / Mohr-Coulomb law." << std::endl;
  const double c = 2e3, phi = 0.4363323129985824, psi = 0.17453292519943295;
  auto law =
      siconos::modeling::DruckerPragerPlasticityNSL::mohrCoulombPlaneStrain(c, phi, psi);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sin(phi), law->eta(), 1e-15);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sin(phi) - std::sin(psi), law->theta(), 1e-15);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0) * c * std::cos(phi), law->cohesion(), 1e-9);
  // random stresses: both criteria give the same yield function up to the factor sqrt2
  for (int k = 0; k < 100; k++) {
    double sxx = 5e3 * std::sin(1.0 + 3 * k), syy = 5e3 * std::cos(2.0 + 5 * k),
           sxy = 3e3 * std::sin(0.5 + 7 * k);
    double sN = -(sxx + syy) / std::sqrt(2.0);
    double sT = std::hypot((sxx - syy) / std::sqrt(2.0), std::sqrt(2.0) * sxy);
    double f_cone = sT - law->eta() * sN - law->cohesion();
    double p = 0.5 * (sxx + syy), R = std::hypot(0.5 * (sxx - syy), sxy);
    double f_mc = R + p * std::sin(phi) - c * std::cos(phi);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0) * f_mc, f_cone, 1e-9 * c);
  }
  std::cout << "✅ Drucker-Prager / Mohr-Coulomb law test ended with success." << std::endl;
}
