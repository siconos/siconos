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
#include "BlockMatrix.hpp"

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

  rhsK.reset(new SimpleMatrix(*K));
  rhsC.reset(new SimpleMatrix(*C));
  minus_inv_M.reset(new SimpleMatrix(3,3));
  (*minus_inv_M)(0,0) = -1.; (*minus_inv_M)(1,1) = -0.5; (*minus_inv_M)(2,2) = -1./3;
  prod(*minus_inv_M, *K, *rhsK, true);
  prod(*minus_inv_M, *C, *rhsC, true);
}

void LagrangianLinearTIDSTest::tearDown()
{}


// Mass, K, C
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::LagrangianLinearTIDS ds(new LagrangianLinearTIDS(q0, velocity0, mass, K, C));
  SiconosVector zero(3);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", Type::value(*ds) == Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->q()) == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->velocity()) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->acceleration() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->mass()) == *(mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->K()) == *(K), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->C()) == *(C), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->p(0) == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->p(2) == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->forces() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->fInt() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->fExt() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->fGyr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->inverseMass() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianFIntq() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianFIntqDot() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianFGyrq() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianFGyrqDot() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianqForces() == K, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->jacobianqDotForces() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->computeKineticEnergy() == 87.0, true);

  double time = 1.;
  ds->initRhs(time);
  SiconosVector x0(*q0, *velocity0);
  SiconosVector acc0(3);
  prod(*K, *q0, acc0, true);
  prod(*C, *velocity0, acc0, false);
  prod(*minus_inv_M, acc0, acc0, true);
  SiconosVector rhs0(*velocity0, acc0);
  
  SP::SiconosMatrix m0(new SimpleMatrix(3,3, Siconos::ZERO));
  SimpleMatrix i0(3,3); // new SimpleMatrix(3,3));//, Siconos::IDENTITY));
  i0(0,0) = i0(1,1) = i0(2,2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->x0()) == x0, true);
  
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0,0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0,1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1,0)) == *rhsK, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1,1)) == *rhsC, true);
  SiconosVector f0(3);
  prod(*K, *q0, f0, true);
  prod(*C, *velocity0, f0, false);
  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");

  time = 1.5;
  ds->computeForces(time, q0, velocity0);
  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  add(f0, -time * *x01, f0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->fExt()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->forces()) == -1.* f0, true);
  ds->computeRhs(time);
  prod(*minus_inv_M, f0, f0, true);
  SiconosVector rhs1(*velocity0, f0); 
  ds->computeJacobianRhsx(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->rhs()) == rhs1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0,0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(0,1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1,0)) == *rhsK, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS1 : ", *(ds->jacobianRhsx()->block(1,1)) == *rhsC, true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}


// Initial conditions and mass
void LagrangianLinearTIDSTest::testBuildLagrangianLinearTIDS2()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::LagrangianLinearTIDS ds(new LagrangianLinearTIDS(q0, velocity0, mass));
  SiconosVector zero(3);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", Type::value(*ds) == Type::LagrangianLinearTIDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->dimension() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->q0() == q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->velocity0() == velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->q()) == *q0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->velocity()) == *velocity0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->acceleration() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->mass()) == *(mass), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->K() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->C() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->p(1)) == zero, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->p(0) == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->p(2) == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->forces() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->fInt() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->fExt() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->fGyr() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->inverseMass() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianFIntq() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianFIntqDot() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianFGyrq() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianFGyrqDot() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianqForces() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->jacobianqDotForces() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->computeKineticEnergy() == 87.0, true);

  double time = 1.;
  ds->initRhs(time);
  SiconosVector x0(*q0, *velocity0);
  SiconosVector rhs0(*velocity0, zero);
  SP::SiconosMatrix m0(new SimpleMatrix(3,3, Siconos::ZERO));
  SimpleMatrix i0(3,3); // new SimpleMatrix(3,3));//, Siconos::IDENTITY));
  i0(0,0) = i0(1,1) = i0(2,2) = 1.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", ds->n() == 2 * 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->x0()) == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->rhs()) == rhs0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(0,0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(0,1)) == i0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(1,0)) == *m0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->jacobianRhsx()->block(1,1)) == *m0, true);
  SiconosVector f0(3);
  ds->setComputeFExtFunction("TestPlugin.so", "computeFExt");
  time = 1.5;
  ds->computeForces(time, q0, velocity0);
  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  add(f0, -time * *x01, f0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->fExt()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIDS2 : ", *(ds->forces()) == -1.* f0, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}


