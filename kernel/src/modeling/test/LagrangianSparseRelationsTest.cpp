/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "LagrangianSparseRelationsTest.hpp"

#include "BlockVector.hpp"
#include "Interaction.hpp"
#include "LagrangianSparseDS.hpp"
#include "LagrangianSparseRheonomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "StorageTools.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianSparseRelationsTest);

void LagrangianSparseRelationsTest::setUp() {}

void LagrangianSparseRelationsTest::tearDown() {}

// Complete test of LagrangianSparseRelations: see CamFollower example.

//  constructor:
void LagrangianSparseRelationsTest::testBuildLagrangianSparseRelationsheonomousR() {
  std::cout << "--> Start rheonomous (sparse storage) relation test" << std::endl;
  auto rel = std::make_shared<siconos::modeling::LagrangianSparseRheonomousR>();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "Type ", rel->getType() == siconos::modeling::RelationType::LagrangianSparse, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "Subtype", rel->getSubType() == siconos::modeling::RelationSubType::RheonomousR, true);
  siconos::algebra::Index ndof = 3;
  auto q0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  auto v0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  *q0 << 1., 2., 3.;
  v0->setZero();
  double tref = 0.4;
  siconos::algebra::SiconosVector href{1};
  href << (*q0)(1) + (*v0)(2) + tref;

  auto hfunc = [](const siconos::algebra::BlockVector& pos, double time,
                  Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result << pos(1) + pos(4) + time;
  };

  auto jachqfunc = [](const siconos::algebra::BlockVector& pos, double time,
                      siconos::algebra::SiconosSparseMatrix& result) {
    // Using coeffRef might be very unefficient, depending on the
    // original structure of the
    // matrix, but this is just for tests.

    result.coeffRef(0, 0) = pos(0);  // just to check that pos is callable ...
    result.coeffRef(0, 2) = 2. * time;
    result.coeffRef(0, 4) = pos(4);
    result.makeCompressed();  // Required
    std::cout << "Compute jacobian h ...\n";
  };

  auto jacref = std::make_shared<siconos::algebra::SiconosSparseMatrix>(1, 6);
  jacref->insert(0, 0) = (*q0)(0);
  jacref->insert(0, 2) = 2. * tref;
  jacref->insert(0, 5) = (*v0)(2);
  jacref->makeCompressed();

  auto hdotfunc = [](const siconos::algebra::BlockVector& pos, double time,
                     Eigen::Ref<siconos::algebra::MapVectorType> result) {
    result << pos(1) + pos(4) + time;
  };

  siconos::algebra::SiconosVector result_y{1};
  result_y << 0.;

  // Create a nsds
  auto lds1 = std::make_shared<siconos::modeling::LagrangianSparseDS>(
      *q0, *v0, siconos::algebra::alias_t);
  auto lds2 = std::make_shared<siconos::modeling::LagrangianSparseDS>(
      *q0, *v0, siconos::algebra::alias_t);
  auto nslaw = std::make_shared<siconos::modeling::NewtonImpactNSL>(0.);
  try {
    auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, rel);
    auto nsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(0., 1.);

    nsds->insertDynamicalSystem(lds1);
    nsds->insertDynamicalSystem(lds2);

    // First test - user-defined function for every operators
    rel->setComputehFunction(hfunc);
    rel->setComputeJacobianhOver_qFunction(jachqfunc);
    rel->setComputehdotFunction(hdotfunc);

    nsds->link(inter, lds1, lds2);

    // rel->initialize(*inter);

    siconos::algebra::BlockVector block{2, ndof};
    block.setVectorPtr(0, q0);
    block.setVectorPtr(1, v0);

    rel->computeh(block, tref, result_y);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - rheonomous (sparse storage) - Compute h",
                                 result_y.isApprox(href), true);

    rel->computehdot(block, tref);

    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - rheonomous (sparse storage) - Compute hdot",
                                 result_y.isApprox(href), true);

    rel->computeJacobianhOver_q(block, tref);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("test - rheonomous (sparse storage) - Compute jacobian h_q",
                                 rel->jacobianhOver_q().isApprox(*jacref), true);
    rel->computeJacobianhOver_q(block, 0.);  // reset

    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - rheonomous (sparse storage) - Compute jacobian h_q (2)",
        rel->jacobianhOver_q().isApprox(*jacref), false);

    // Set consant value - Copy
    rel->setConstantJacobianhOver_q(*jacref, siconos::algebra::copy_t);
    rel->initialize(*inter);
    (*q0)(1) *= 4;
    rel->computeJacobianhOver_q(block, 1.3);  // this should not have any effect
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - rheonomous (sparse storage) - constant jacobian h_q (copy)",
        rel->jacobianhOver_q().isApprox(*jacref), true);

    // back to setcompute, just to check
    auto jachqfunc2 = [](const siconos::algebra::BlockVector& pos, double time,
                         siconos::algebra::SiconosSparseMatrix& result) {
      // Using coeffRef might be very unefficient, depending on the
      // original structure of the
      // matrix, but this is just for tests.

      result.coeffRef(0, 0) = 2 * pos(0);  // just to check that pos is callable ...
      result.coeffRef(0, 2) = 4. * time;
      result.coeffRef(0, 4) = 2 * pos(4);
      result.makeCompressed();  // Required
      std::cout << "Compute jacobian h v2 ...\n";
    };

    rel->setComputeJacobianhOver_qFunction(jachqfunc2);
    rel->initialize(*inter);
    rel->computeJacobianhOver_q(block, tref);
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - rheonomous (sparse storage) -  Compute jacobian h_q (3)",
        rel->jacobianhOver_q().isApprox(2. * *jacref), true);

    // Set consant value - alias

    auto jacmap = std::make_shared<Eigen::Map<siconos::algebra::SiconosSparseMatrix>>(
        jacref->rows(), jacref->cols(), jacref->nonZeros(), jacref->outerIndexPtr(),
        jacref->innerIndexPtr(), jacref->valuePtr());

    rel->setConstantJacobianhOver_q(*jacmap, siconos::algebra::alias_t);
    rel->initialize(*inter);
    rel->computeJacobianhOver_q(block, 1.3);  // this should not have any effect
    CPPUNIT_ASSERT_EQUAL_MESSAGE(
        "test - rheonomous (sparse storage) - constant jacobian h_q (alias)",
        rel->jacobianhOver_q().isApprox(*jacref), true);

  } catch (...) {
    siconos::exception::process();
  }

  // For a complete test/example with initialize(inter), compute call and so on
  // see CamFollower-Rheonomous
  std::cout << "✅  rheonomous (sparse storage) relation test ok" << std::endl;
}
