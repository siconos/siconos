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


// constructor from data
void FirstOrderLinearDSTest::testBuildFirstOrderLinearDS1()
{
  std::cout << "--> Test: constructor 1." <<std::endl;
  SP::FirstOrderLinearDS ds(new FirstOrderLinearDS(x0, "TestPlugin:computeA", "TestPlugin:computeb"));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1A : ", Type::value(*ds) == Type::FirstOrderLinearDS, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1B : ", ds->n() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1C : ", ds->x0() == x0, true);

  double time = 1.5;
  ds->initialize(time);
  ds->computeA(time);
  ds->computeb(time);
  ds->computeRhs(time);
  SP::SiconosVector x01(new SiconosVector(3));
  (*x01)(0) = 0;
  (*x01)(1) = 1;
  (*x01)(2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1D : ", *(ds->b()) == time* *x01, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1E : ", *(ds->A()) == 2 * *A0, true);
  //  ds->rhs()->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS1F : ", *(ds->rhs()) == (time* *x01 + 2 * prod(*A0, *x0)) , true);

  ds->computeRhs(time);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildFirstOrderLinearDS2M : ", *(ds->rhs()) == (2 * prod(*A0 , *x0) + * (ds->b())), true);
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
  ds1->setb(b0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBPtr : ", ds1->b() == b0, true);
  std::cout << "--> setBPtr test ended with success." <<std::endl;
}

// plugins: plugins loading is already in testBuildFirstOrderLinearDS2

void FirstOrderLinearDSTest::End()
{
  std::cout << "============================================" <<std::endl;
  std::cout << " ===== End of FirstOrderLinearDS tests =====" <<std::endl;
  std::cout << "============================================" <<std::endl;
}
