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
#ifndef __SimpleMatrixTest__
#define __SimpleMatrixTest__

#include <cppunit/extensions/HelperMacros.h>

#include "BlockMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

class SimpleMatrixTest : public CppUnit::TestFixture {
 private:
  ACCEPT_SERIALIZATION(SimpleMatrixTest);

  // test suite
  CPPUNIT_TEST_SUITE(SimpleMatrixTest);

  CPPUNIT_TEST(testNormInf);
  CPPUNIT_TEST(testSetBlock);
  CPPUNIT_TEST(testSetBlock2);
  CPPUNIT_TEST(testOperators6Bis);
  CPPUNIT_TEST(testOperators8_5);
  CPPUNIT_TEST(testOperators8_6);
  CPPUNIT_TEST(testOperators9);
  CPPUNIT_TEST(testProd);
  CPPUNIT_TEST(testProdBis);
  CPPUNIT_TEST(testProd4);
  CPPUNIT_TEST(testProd5);
  CPPUNIT_TEST(testProd6);
  // CPPUNIT_TEST(testGemv);
  // CPPUNIT_TEST(testGemm);
  // CPPUNIT_TEST(testPLUFactorizationInPlace);
  // CPPUNIT_TEST(testFactorize);
  CPPUNIT_TEST(End);
  CPPUNIT_TEST_SUITE_END();

  void testNormInf();
  void testSetBlock();
  void testSetBlock2();
  void testOperators6Bis();
  void testOperators8_5();
  void testOperators8_6();
  void testOperators9();
  void testProd();
  void testProdBis();
  void testProdTer();
  void testProd4();
  void testProd5();
  void testProd6();
  // void testGemm();
  // void testGemv();
  // void testFromAndFillCSC();
  // void testPLUFactorizationInPlace();
  // void testFactorize();
  // void testSolve();
  void End();

  unsigned int size, size2;
  std::shared_ptr<siconos::algebra::SiconosMatrix> SicM, m1, m2, m3, m4, m5, m6, m7, m8, C, Cb,
      Cb2;
  std::shared_ptr<const siconos::algebra::SiconosMatrix> A, B, Ab, Bb;
  std::shared_ptr<siconos::algebra::SiconosMatrix> SimM;
  std::string fic1, fic2;
  std::shared_ptr<siconos::algebra::SiconosVector> vect1, vect2, vect3, vect_sparse_1;
  siconos::algebra::SiconosMatrix D{2, 2};
  double tol;

 public:
  void setUp();
  void tearDown();
};

#endif
