
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
#include "FirstOrderLinearDSTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(FirstOrderLinearDSTest);


void FirstOrderLinearDSTest::setUp()
{
  x0.reset(new SiconosVector(3));
  (*x0)(0) = 1;
  (*x0)(1) = 2;
  (*x0)(2) = 3;

  b0.reset(new SiconosVector(3));
  (*b0)(0) = 4;
  (*b0)(1) = 5;
  (*b0)(2) = 6;

  A0.reset(new SimpleMatrix("matA0.dat", true));
}
void FirstOrderLinearDSTest::tearDown()
{}

// constructor from initial state only
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS0()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(x0));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->n() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->M() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->invM() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->b() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->A() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->f() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->jacobianfx() == NULL, true);
  double time = 1.5;
  SiconosVector zero(3);
  ds->initRhs(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", *(ds->rhs()) == zero , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS0 : ", ds->jacobianRhsx() == NULL, true);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeM(time);

  ds->setComputeMFunction("TestPlugin", "computeM");
  ds->computeM(time);
  SimpleMatrix Mref(3,3);
  Mref(0,0) = 1. * time; Mref(1,1) = 2. * time; Mref(2,2) = 3. * time;
  std::cout << "MLMLMQLSQML " << std::numeric_limits<double>::epsilon() << std::endl ;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", *(ds->M()) == Mref, true);
  std::cout << "--> Constructor 0 test ended with success." <<std::endl;
}

// constructor from initial state and plugins
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(x0, "TestPlugin:computeA", "TestPlugin:computeb"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->n() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->x0() == x0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->M() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->invM() == NULL, true);

  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;
  double time = 1.5;
  ds->computeA(time);
  ds->computeb(time);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->b()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->A()) == 2 * *A0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", ds->f() == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->jacobianfx()) == *(ds->A()), true);
  ds->setComputeMFunction("TestPlugin", "computeM");
  ds->computeM(time);
  SimpleMatrix Mref(3,3);
  Mref(0,0) = 1. * time; Mref(1,1) = 2. * time; Mref(2,2) = 3. * time; 
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderNonLinearDS1 : ", *(ds->M()) == Mref, true);
  ds->initRhs(time);
  SimpleMatrix invM(3,3);
  invM(0,0) = 1. / time; invM(1,1) = 1./ (2. * time); invM(2,2) = 1./(3. * time);
  SiconosVector tmp(3);
  tmp = (time* *x01 + 2 * prod(*A0, *x0));
  prod(invM, tmp, tmp);
  prod(invM, 2* *A0, Mref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->rhs()) == tmp , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->jacobianRhsx()) == Mref, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->b()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1 : ", *(ds->A()) == 2 * *A0, true);
  
  std::cout << "--> Constructor 3 test ended with success." <<std::endl;
}

// setAPtr
void FirstOrderLinearDSTest::testSetAPtr()
{
  std::cout << "--> Test: setAPtr." <<std::endl;
  SP::FirstOrderLinearDS ds1(new FirstOrderLinearDS(x0));
  ds1->setAPtr(A0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetAPtr : ", ds1->A() == A0, true);
  std::cout << "--> setAPtr test ended with success." <<std::endl;
}

// setBPtr
void FirstOrderLinearDSTest::testSetBPtr()
{
  std::cout << "--> Test: setBPtr." <<std::endl;
  SP::FirstOrderLinearDS ds1(new FirstOrderLinearDS(x0));
  ds1->setbPtr(b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ds1->b() == b0, true);
  std::cout << "--> setBPtr test ended with success." <<std::endl;
}
