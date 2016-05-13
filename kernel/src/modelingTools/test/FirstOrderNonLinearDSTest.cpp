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
#include "FirstOrderNonLinearDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderNonLinearDSTest);


void FirstOrderNonLinearDSTest::setUp()
{
  x0.reset(new SiconosVector(3));
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  J0.reset(new SimpleMatrix("matJ0.dat", true));
  M.reset(new SimpleMatrix("matM.dat", true));

}

void FirstOrderNonLinearDSTest::tearDown()
{}

// constructor from data
void FirstOrderNonLinearDSTest::testBuildFirstOrderNonLinearDS3()
{
  std::cout << "--> Test: constructor 3." <<std::endl;
  SP::FirstOrderNonLinearDS ds(new FirstOrderNonLinearDS(*x0, "TestPlugin:computeF", "TestPlugin:computeJacobianfx"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3A : ", Type::value(*ds) == Type::FirstOrderNonLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3C : ", ds->getN() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3D : ", ds->getX0() == *x0, true);
  double time = 1.5;
  ds->initialize("TimeStepping", time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS3E : ", ds->getRhs() == time* *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2H : ", ds->getF() == time* *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", ds->getJacobianfx() == *J0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS2I : ", !ds->M(), false);
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// setX0
void FirstOrderNonLinearDSTest::testSetX0()
{
  std::cout << "--> Test: setX0." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX0(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0 : ", ds1->getX0() == *x0, true);
  std::cout << "--> setX0 test ended with success." <<std::endl;
}

// setX0Ptr
void FirstOrderNonLinearDSTest::testSetX0Ptr()
{
  std::cout << "--> Test: setX0Ptr." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX0Ptr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX0Ptr : ", ds1->getX0() == *x0, true);
  std::cout << "--> setX0Ptr test ended with success." <<std::endl;
}

// setX
void FirstOrderNonLinearDSTest::testSetx()
{
  std::cout << "--> Test: setX." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetX : ", ds1->getx() == *x0, true);
  std::cout << "--> setX test ended with success." <<std::endl;
}

// setXPtr
void FirstOrderNonLinearDSTest::testSetXPtr()
{
  std::cout << "--> Test: setXPtr." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setXPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetXPtr : ", ds1->getx() == *x0, true);
  std::cout << "--> setXPtr test ended with success." <<std::endl;
}

// setR
void FirstOrderNonLinearDSTest::testSetR()
{
  std::cout << "--> Test: setR." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setR(*x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetR : ", ds1->getR() == *x0, true);
  std::cout << "--> setR test ended with success." <<std::endl;
}

// setRPtr
void FirstOrderNonLinearDSTest::testSetRPtr()
{
  std::cout << "--> Test: setRPtr." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setRPtr(x0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetRPtr : ", ds1->getR() == *x0, true);
  std::cout << "--> setRPtr test ended with success." <<std::endl;
}

// set JacobianX
void FirstOrderNonLinearDSTest::testSetJacobianfx()
{
  std::cout << "--> Test: setJacobianfx." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setJacobianfx(*J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianX : ", ds1->getJacobianfx() == *J0, true);
  std::cout << "--> setJacobianfx test ended with success." <<std::endl;
}

// setJacobianXPtr
void FirstOrderNonLinearDSTest::testSetJacobianfxPtr()
{
  std::cout << "--> Test: setJacobianfxPtr." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setJacobianfxPtr(J0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetJacobianfxPtr : ", ds1->getJacobianfx() == *J0, true);
  std::cout << "--> setJacobianfxPtr test ended with success." <<std::endl;
}

// init
void FirstOrderNonLinearDSTest::testInitMemory()
{
  std::cout << "--> Test: initMemory." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->initMemory(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem1 : ", ds1->xMemory()->getMemorySize() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem3 : ", ds1->rMemory()->getMemorySize() == 2, true);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem4 : ", ds1->xMemory()->getNbVectorsInMemory() == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testInitMem6 : ", ds1->rMemory()->getNbVectorsInMemory() == 0, true);
  std::cout << "--> initMemory test ended with success." <<std::endl;
}


// swap
void FirstOrderNonLinearDSTest::testSwap()
{
  std::cout << "--> Test: swap." <<std::endl;
  SP::FirstOrderNonLinearDS ds1(new FirstOrderNonLinearDS(tmpxml2));
  ds1->setX(*x0);
  ds1->setR(*x0);
  ds1->initMemory(1);
  ds1->swapInMemory();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap1 : ", *((ds1->xMemory()->getVectorMemory())[0]) == *x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSwap3 : ", *((ds1->rMemory()->getVectorMemory())[0]) == *x0, true);
  std::cout << "--> swap test ended with success." <<std::endl;
}


// plugins: plugins loading is already in testBuildFirstOrderNonLinearDS2

void FirstOrderNonLinearDSTest::End()
{
  std::cout << "===============================================" <<std::endl;
  std::cout << " ===== End of FirstOrderNonLinearDS tests =====" <<std::endl;
  std::cout << "===============================================" <<std::endl;
}
