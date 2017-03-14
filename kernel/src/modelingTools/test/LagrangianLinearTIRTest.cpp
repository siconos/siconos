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
#include "LagrangianLinearTIRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianLinearTIRTest);


void LagrangianLinearTIRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SiconosVector(1));
  (*e)(0) = 0.1;
}

void LagrangianLinearTIRTest::tearDown()
{}

// data constructor (1)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR1d : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, F, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR2g : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}

// data constructor (5)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR3g : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// data constructor (4)
void LagrangianLinearTIRTest::testBuildLagrangianLinearTIR4()
{
  std::cout << "--> Test: constructor 4." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianLinearTIR4d : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 4 test ended with success." <<std::endl;
}
// set C as a matrix and then plug it


// setCPtr
void LagrangianLinearTIRTest::testSetCPtr()
{
  std::cout << "--> Test: setCPtr." <<std::endl;
  SP::SimpleMatrix tmp(new SimpleMatrix(*C));
  tmp->zero();
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(tmp));
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C() == C, true);
  std::cout << "--> setCPtr test ended with success." <<std::endl;
}

// setFPtr
void LagrangianLinearTIRTest::testSetFPtr()
{
  std::cout << "--> Test: setFPtr." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F() == F, true);
  std::cout << "--> setFPtr test ended with success." <<std::endl;
}

// set E

// setEPtr
void LagrangianLinearTIRTest::testSetEPtr()
{
  std::cout << "--> Test: setEPtr." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e() == e, true);
  std::cout << "--> setEPtr test ended with success." <<std::endl;
}



void LagrangianLinearTIRTest::testGetJacPtr()
{
  std::cout << "--> Test: jac." <<std::endl;
  SP::LagrangianLinearTIR folr(new LagrangianLinearTIR(C));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJachq: ", folr->jachq() == C, true);

  std::cout << "--> testGetJacPtr test ended with success." <<std::endl;
}
