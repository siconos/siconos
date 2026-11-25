// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SolidLinearTIDSTest.hpp"

#include "FENode.hpp"
#include "FemTools.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"
#include "SolidLinearTIDS.hpp"
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
