/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "LagrangianDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianDSTest);


void LagrangianDSTest::setUp()
{
  q0.reset(new SiconosVector(3));
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0.reset(new SiconosVector(3));
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  u0.reset(new SiconosVector(2));
  (*u0)(0) = 4;
  (*u0)(1) = 5;

  mass.reset(new SimpleMatrix(3, 3));
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;

}

void LagrangianDSTest::tearDown()
{}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS4()
{
  std::cout << "--> Test: constructor 4." <<std::endl;

  SP::LagrangianDS ds(new LagrangianDS(q0, velocity0, mass));
  double time = 1.5;
  ds->initialize(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4A : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4B : ", ds->number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4D : ", ds->ndof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS4L : ", *(ds->mass()) == *(mass), true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildLagrangianDS4L : ", ds->computeKineticEnergy() == 87.0, true);
  std::cout << "--> Constructor 4 test ended with success." <<std::endl;
}

// constructor from data
void LagrangianDSTest::testBuildLagrangianDS5()
{
  std::cout << "--> Test: constructor 5." <<std::endl;
  std::string plugin = "TestPlugin:computeMass";
  SP::DynamicalSystem ds(new LagrangianDS(q0, velocity0, plugin));
  double time = 1.5;
  ds->initialize(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5A : ", Type::value(*ds) == Type::LagrangianDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5B : ", ds->number() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5D : ", std11::static_pointer_cast<LagrangianDS>(ds)->ndof() == 3, true);

  std11::static_pointer_cast<LagrangianDS>(ds)->computeMass();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianDS5L : ", *(std11::static_pointer_cast<LagrangianDS>(ds)->mass()) == (*mass), true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testBuildLagrangianDS5L: ", std11::static_pointer_cast<LagrangianDS>(ds)->computeKineticEnergy() == 87.0, true);
  std::cout << "--> Constructor 5 test ended with success." <<std::endl;
}

void LagrangianDSTest::testcomputeDS()
{
  // std::cout << "-->Test: computeDS." <<std::endl;
  // DynamicalSystem * ds(new LagrangianDS(tmpxml2));
  // SP::LagrangianDS copy =  std11::static_pointer_cast<LagrangianDS>(ds);
  // double time = 1.5;
  // ds->initialize("EventDriven", time);
  // ds->computeRhs(time);
  // std::cout << "-->Test: computeDS." <<std::endl;
  // ds->computeJacobianRhsx(time);
  // std::cout << "-->Test: computeDS." <<std::endl;
  // SimpleMatrix M(3, 3);
  // M(0, 0) = 1;
  // M(1, 1) = 2;
  // M(2, 2) = 3;
  // SP::SiconosMatrix jx = ds->jacobianRhsx();
  // SP::SiconosVector vf = ds->rhs();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->vector(0)) == *velocity0, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(M, *(vf->vector(1))) == (copy->getFExt() - copy->getFInt() - copy->getFGyr()) , true);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 0))) == (copy->getJacobianFL(0)) , true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(M, *(jx->block(1, 1))) == (copy->getJacobianFL(1)) , true);
  // std::cout << "--> computeDS test ended with success." <<std::endl;

}
void LagrangianDSTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of LagrangianDS tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
