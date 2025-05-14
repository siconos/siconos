/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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
#ifndef __SiconosMatrixTest__
#define __SiconosMatrixTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class SiconosMatrixTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SiconosMatrixTest);

  // test suite
  CPPUNIT_TEST_SUITE(SiconosMatrixTest);

  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testSymm);
  CPPUNIT_TEST(testfillTriplet);
  CPPUNIT_TEST(testSetBlock);
  CPPUNIT_TEST(testProd);
  CPPUNIT_TEST(testProd2);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testNormInf();
  void testSymm();
  void testfillTriplet();
  void testSetBlock();
  void testProd();
  void testProd2();
  void End();

  unsigned int size;
  std::shared_ptr<siconos::algebra::SiconosDenseMatrix> SicM;
  siconos::algebra::SiconosSparseMatrix Asparse;
  std::shared_ptr<const siconos::algebra::SiconosDenseMatrix> A;
  std::string fic1, fic2;
  double tol;

 public:
  void setUp();
  void tearDown();
};

#endif
