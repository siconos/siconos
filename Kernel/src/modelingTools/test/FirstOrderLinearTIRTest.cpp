/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "FirstOrderLinearTIRTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearTIRTest);


void FirstOrderLinearTIRTest::setUp()
{
  C.reset(new SimpleMatrix("matC.dat", true));
  D.reset(new SimpleMatrix("matD.dat", true));
  B.reset(new SimpleMatrix("matB.dat", true));
  F.reset(new SimpleMatrix("matF.dat", true));
  e.reset(new SiconosVector(1));
  (*e)(0) = 0.1;
}

void FirstOrderLinearTIRTest::tearDown()
{}

// data constructor (1)
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, B));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1b : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1c : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR1d : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 1 test ended with success." <<std::endl;
}

// data constructor (2)
void FirstOrderLinearTIRTest::testBuildFirstOrderLinearTIR2()
{
  std::cout << "--> Test: constructor 2." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, D, F, e, B));
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2a : ", folr->C() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2b : ", folr->D() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2c : ", folr->F() == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2d : ", folr->e() == e, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2e : ", folr->B() == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2f : ", folr->getType() == RELATION::FirstOrder, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearTIR2g : ", folr->getSubType() == RELATION::LinearTIR, true);
  std::cout << "--> Constructor 2 test ended with success." <<std::endl;
}


// set C as a matrix and then plug it


// setCPtr
void FirstOrderLinearTIRTest::testSetCPtr()
{
  std::cout << "--> Test: setCPtr." <<std::endl;
  SP::SiconosMatrix tmp(new SimpleMatrix(*C));
  tmp->zero();
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(tmp, B));
  folr->setCPtr(C);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetCPtr : ", folr->C() == C, true);
  std::cout << "--> setCPtr test ended with success." <<std::endl;
}

// set D

// setDPtr
void FirstOrderLinearTIRTest::testSetDPtr()
{
  std::cout << "--> Test: setDPtr." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, B));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetDPtr: ", folr->D() == D, true);
  std::cout << "--> setDPtr test ended with success." <<std::endl;
}

// set F

// setFPtr
void FirstOrderLinearTIRTest::testSetFPtr()
{
  std::cout << "--> Test: setFPtr." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, B));
  folr->setFPtr(F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetFPtr: ", folr->F() == F, true);
  std::cout << "--> setFPtr test ended with success." <<std::endl;
}

// set E


// setEPtr
void FirstOrderLinearTIRTest::testSetEPtr()
{
  std::cout << "--> Test: setEPtr." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, B));
  folr->setEPtr(e);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetEPtr: ", folr->e() == e, true);
  std::cout << "--> setEPtr test ended with success." <<std::endl;
}

// set B


// setBPtr
void FirstOrderLinearTIRTest::testSetBPtr()
{
  std::cout << "--> Test: setBPtr." <<std::endl;
  SP::SiconosMatrix tmp(new SimpleMatrix(*B));
  tmp->zero();
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, tmp));
  folr->setBPtr(B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr: ", folr->B() == B, true);
  std::cout << "--> setBPtr test ended with success." <<std::endl;
}



void FirstOrderLinearTIRTest::testGetJacPtr()
{
  std::cout << "--> Test: jac." <<std::endl;
  SP::FirstOrderLinearTIR folr(new FirstOrderLinearTIR(C, B));
  folr->setDPtr(D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJach: ", folr->jachx() == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJach: ", folr->jachlambda() == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetJach: ", folr->jacglambda() == B, true);

  std::cout << "--> setBPtr test ended with success." <<std::endl;
}

void FirstOrderLinearTIRTest::End()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== End of FirstOrderLinearTIR Tests ===== " <<std::endl;
  std::cout << "=========================================== " <<std::endl;
}
