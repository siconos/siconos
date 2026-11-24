// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "SolidLinearTIDSTest.hpp"

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
