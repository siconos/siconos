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
#include "LagrangianLinearTIDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIDSTest);


void LagrangianLinearTIDSTest::setUp()
{
  q0.reset(new SiconosVector(3));
  (*q0)(0) = 1;
  (*q0)(1) = 2;
  (*q0)(2) = 3;
  velocity0.reset(new SiconosVector(3));
  (*velocity0)(0) = 4;
  (*velocity0)(1) = 5;
  (*velocity0)(2) = 6;

  mass.reset(new SimpleMatrix(3, 3));
  (*mass)(0, 0) = 1;
  (*mass)(1, 1) = 2;
  (*mass)(2, 2) = 3;

  K.reset(new SimpleMatrix("K.dat", true));
  C.reset(new SimpleMatrix("C.dat", true));




void LagrangianLinearTIDSTest::tearDown()
{}


// Mass, K, C
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  SP::LagrangianLinearTIDS ds(new LagrangianLinearTIDS(8, *q0, *velocity0, *mass, *K, *C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2A : ", Type::value(*ds) == Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2B : ", ds->number() == 8, true);
  std::cout << "--> Test: constructor 2." <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getStepsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2F : ", ds->getVelocity0() == *velocity0, true);
  std::cout << "--> Test: constructor 2." <<std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2G : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2H : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getMass() == *mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2J : ", ds->getK() == *K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2K : ", ds->getC() == *C, true);
  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");

  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getFExt() == time* *x01, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}

// only  Mass
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;
  SP::LagrangianLinearTIDS ds(new LagrangianLinearTIDS(8, *q0, *velocity0, *mass));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2A : ", Type::value(*ds) == Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2B : ", ds->number() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getStepsInMemory() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2D : ", ds->getNdof() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2E : ", ds->getQ0() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2F : ", ds->getVelocity0() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2G : ", ds->getQ() == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2H : ", ds->getVelocity() == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getMass() == *mass, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2J : ", !ds->K(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2K : ", !ds->C(), false);

  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");

  double time = 1.5;
  ds->initialize("TimeStepping", time);

  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2I : ", ds->getFExt() == time* *x01, true);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

void LagrangianLinearTIDSTest::testcomputeDS()
{
  std::cout << "-->Test: computeDS." <<std::endl;
  SP::DynamicalSystem ds(new LagrangianLinearTIDS(tmpxml1));
  SP::LagrangianLinearTIDS copy = std11::static_pointer_cast<LagrangianLinearTIDS>(ds);
  double time = 1.5;
  ds->initialize("EventDriven", time);
  SP::SiconosMatrix jx = ds->jacobianRhsx();
  SP::SiconosVector vf = ds->rhs();

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSI : ", *(vf->vector(0)) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSJ : ", prod(*mass, *(vf->vectorPtr(1))) == (copy->getFExt() - prod(*K, *(copy->q())) - prod(*C, *(copy->getVelocity()))) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(*mass, *(jx->block(1, 0))) == (-1.0 * *K) , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testComputeDSL : ", prod(*mass, *(jx->block(1, 1))) == (-1.0 * *C) , true);
  std::cout << "--> computeDS test ended with success." <<std::endl;


}
void LagrangianLinearTIDSTest::End()
{
  std::cout << "==============================================" <<std::endl;
  std::cout << " ===== End of LagrangianLinearTIDS tests =====" <<std::endl;
  std::cout << "==============================================" <<std::endl;
}
