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
#include "LagrangianCompliantLinearTIRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(LagrangianCompliantLinearTIRTest);


void LagrangianCompliantLinearTIRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  D.reset(new SimpleMatrix("matD.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SiconosVector(1));
  (*e)(0) = 0.1;
}

void LagrangianCompliantLinearTIRTest::tearDown()
{}

// data constructor (1)
void LagrangianCompliantLinearTIRTest::testBuildLagrangianCompliantLinearTIR1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR1a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR1c : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR1d : ", folr->getSubType() == RELATION::CompliantLinearTIR, true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

// data constructor (5)
void LagrangianCompliantLinearTIRTest::testBuildLagrangianCompliantLinearTIR2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D, F, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR2f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR2g : ", folr->getSubType() == RELATION::CompliantLinearTIR, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}

// data constructor (5)
void LagrangianCompliantLinearTIRTest::testBuildLagrangianCompliantLinearTIR3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D, e));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR3a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR3d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR3f : ", folr->getType() == RELATION::Lagrangian, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildLagrangianCompliantLinearTIR3g : ", folr->getSubType() == RELATION::CompliantLinearTIR, true);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// setCPtr
void LagrangianCompliantLinearTIRTest::testSetCPtr()
{
  std::cout << "--> Test: setCPtr." <<std::endl;
  SP::SimpleMatrix tmpC(new SimpleMatrix(*C));
  tmpC->zero();
  SP::SimpleMatrix tmpD(new SimpleMatrix(*D));
  tmpD->zero();
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(tmpC, tmpD));
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C() == C, true);
  std::cout << "--> setCPtr test ended with success." <<std::endl;
}

// set D

// setDPtr
void LagrangianCompliantLinearTIRTest::testSetDPtr()
{
  std::cout << "--> Test: setDPtr." <<std::endl;
  SP::SimpleMatrix tmp(new SimpleMatrix(*D));
  tmp->zero();
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, tmp));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->D() == D, true);
  std::cout << "--> setDPtr test ended with success." <<std::endl;
}

// set F

// setFPtr
void LagrangianCompliantLinearTIRTest::testSetFPtr()
{
  std::cout << "--> Test: setFPtr." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D));
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F() == F, true);
  std::cout << "--> setFPtr test ended with success." <<std::endl;
}

// set E

// setEPtr
void LagrangianCompliantLinearTIRTest::testSetEPtr()
{
  std::cout << "--> Test: setEPtr." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D));
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e() == e, true);
  std::cout << "--> setEPtr test ended with success." <<std::endl;
}



void LagrangianCompliantLinearTIRTest::testGetJacPtr()
{
  std::cout << "--> Test: jac." <<std::endl;
  SP::LagrangianCompliantLinearTIR folr(new LagrangianCompliantLinearTIR(C, D));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJachq: ", folr->jachq() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJachlambda: ", folr->jachlambda() == D, true);

  std::cout << "--> setBPtr test ended with success." <<std::endl;
}

void LagrangianCompliantLinearTIRTest::End()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== End of LagrangianCompliantLinearTIR Tests ===== " <<std::endl;
  std::cout << "=========================================== " <<std::endl;
}
