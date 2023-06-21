/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "CableDSTest.hpp"

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosVectorIterator.hpp"
#include "SiconosAlgebraTypes.hpp"  // for UblasType

//#include "ioMatrix.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)		\
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(CableDSTest);

void CableDSTest::setUp()
{
  ndof = 10;
  q0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  v0 = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  mass = std::make_shared<siconos::algebra::SimpleMatrix>(ndof, ndof, siconos::algebra::UblasType::SPARSE);
  // TODO : read a specific reference problem somewhere ...
}

void CableDSTest::tearDown() {}

void CableDSTest::testNoFext()
{

  auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testNoFext: ", cable->fExt() == nullptr, true);
  // ...
}

void CableDSTest::testConstantFext()
{
  auto externalForces = std::make_shared<siconos::algebra::SiconosVector>(ndof);
  // fill external forces as you want ...
  for (auto& val : *externalForces) {
    val = 12;
  }

  auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass);
  cable->setFExtPtr(externalForces);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 1: ", cable->fExt() != nullptr, true);
  auto fext = cable->fExt();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 2: ", fext != nullptr, true);
  for (auto val : *fext) {
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testCFext 3: ", val == 12, true);
  }
  // ...
}

void CableDSTest::testVariableFext()

{
  // A lambda function example, used to compute external forces
  auto myforces = [](double time, std::shared_ptr<siconos::algebra::SiconosVector> result) {
    assert(result);
    for (auto& v : *result)
      v = cos(time);
  };

  auto cable = std::make_shared<siconos::mechanics::fem::CableDS>(q0, v0, mass, myforces);
  auto fext = cable->fExt();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 1: ", fext != nullptr, true);

  cable->computeFExt(3.);
  for (auto val : *fext)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 2: ", val == cos(3.), true);

  auto positions = cable->q();
  auto velocities = cable->velocity();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 4: ", positions->size() == cable->dimension(),
                               true);
  cable->computeForces(5., positions, velocities);
  for (auto val : *fext)
    CPPUNIT_ASSERT_EQUAL_MESSAGE(" testVFext 5: ", val == cos(5.), true);

  // Example on how to proceed to compare matrices:
  // matrix is the new matrix to be checked. SomeReferenceFile.ref contains a former matrix
  // saved in a file and used as reference.
  // auto matrix = mass;
  // double eps = 1e-12;
  // auto error = ioMatrix::compareRefFile(*matrix, "SomeReferenceFile.ref", eps);
  // Check error ...
}
