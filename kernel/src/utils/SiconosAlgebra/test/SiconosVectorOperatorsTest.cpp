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
#include "SiconosConfig.h"

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosAlgebra.hpp"
#include "SiconosVectorOperatorsTest.hpp"
#include "TypeName.hpp"
#include "SiconosVectorOperators.hpp"

using namespace boost::numeric::ublas;

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosVectorOperatorsTest);

void SiconosVectorOperatorsTest::setUp()
{
}

void SiconosVectorOperatorsTest::tearDown()
{}


void SiconosVectorOperatorsTest::testSetValue()
{
  std::cout  << "--> Test: SetValue." << std::endl;
  SiconosVector v(1);
  apply_visitor<SetValue>(storage(v), 0, 1.0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("SetValue: v(0) == 1.0", (apply_visitor<GetValue, double>(storage(v), 0)) == 1.0, true);
}

void SiconosVectorOperatorsTest::testGetValue()
{
  std::cout  << "--> Test: GetValue." << std::endl;
  SiconosVector v(3);
  v(0) = 1.0; v(1) = 2.0; v(2) = 3.0;
  std::cout << Type::name(storage(v)) << std::endl;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("GetValue: v(0) == 1.0", (apply_visitor<GetValue, double>(storage(v), 0)) == 1.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("GetValue: v(1) == 2.0", (apply_visitor<GetValue, double>(storage(v), 1)) == 2.0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("GetValue: v(2) == 3.0", (apply_visitor<GetValue, double>(storage(v), 2)) == 3.0, true);
}

void SiconosVectorOperatorsTest::testCopy()
{
  std::cout  << "--> Test: GetCopy." <<std::endl;

  SiconosVector v(3);
  v(0) = 1.0; v(1) = 2.0; v(2) = 3.0;

  SiconosVector w;

  apply_visitor<Copy>(storage(v), storage(w));

  CPPUNIT_ASSERT_EQUAL_MESSAGE("GetValue: w(0) == v(0)", (apply_visitor<GetValue, double>(storage(w), 0) == apply_visitor<GetValue, double>(storage(v), 0)), true);
}
