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
#ifndef __SiconosVectorTest__
#define __SiconosVectorTest__

#include <cppunit/extensions/HelperMacros.h>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include <cmath>
#include <vector>

class SiconosVectorOperatorsTest : public CppUnit::TestFixture
{


private:
  CPPUNIT_TEST_SUITE(SiconosVectorOperatorsTest);

  CPPUNIT_TEST(testSetValue);
  CPPUNIT_TEST(testGetValue);
  CPPUNIT_TEST(testCopy);

  CPPUNIT_TEST_SUITE_END();

  void testSetValue();
  void testGetValue();
  void testCopy();

public:
  void setUp();
  void tearDown();

};

#endif




