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

#include "BlockMatrixTest.hpp"

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <string>

#include "SiconosAlgebraTypes.hpp"
#include "SiconosConfig.h"
#include "SiconosMatrixOp.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(BlockMatrixTest);

using BlockMatrix = siconos::algebra::BlockMatrix;

void BlockMatrixTest::setUp() {
  tol = 1e-12;

  B = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2, 1);
  C = std::make_shared<siconos::algebra::SimpleMatrix>(2, 4, 2);
  D = std::make_shared<siconos::algebra::SimpleMatrix>(2, 1, 3);
  E = std::make_shared<siconos::algebra::SimpleMatrix>(3, 2, 4);
  F = std::make_shared<siconos::algebra::SimpleMatrix>(3, 4, 5);
  G = std::make_shared<siconos::algebra::SimpleMatrix>(3, 1, 6);

  m.resize(6);
  m[0] = B;
  m[1] = C;
  m[2] = D;
  m[3] = E;
  m[4] = F;
  m[5] = G;

  tRow.resize(2);
  tRow[0] = 2;
  tRow[1] = 5;
  tCol.resize(3);
  tCol[0] = 2;
  tCol[1] = 6;
  tCol[2] = 7;

  mapRef = std::make_shared<BlocksMatrix>(2, 1);
  (*mapRef)(0, 0) = B;
  (*mapRef)(1, 0) = E;
}

void BlockMatrixTest::tearDown() { m.clear(); }

//______________________________________________________________________________

void BlockMatrixTest::testConstructor0()  // constructor with a vector of
                                          // std::shared_ptr<SiconosMatrix>
{
  std::cout << "====================================" << std::endl;
  std::cout << "=== Block Matrix tests start ...=== " << std::endl;
  std::cout << "====================================" << std::endl;
  std::cout << "--> Test: constructor 0." << std::endl;
  auto test = std::make_shared<BlockMatrix>(m, 2, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(0) == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(0, 2) == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 1) == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->block(1, 2) == G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", *test->tabCol() == tCol, true);
  std::cout << "--> Constructor 0 test ended with success." << std::endl;
}

void BlockMatrixTest::testConstructor1()  // Copy constructor, from a BlockMatrix
{
  std::cout << "--> Test: constructor 1." << std::endl;
  auto ref = std::make_shared<BlockMatrix>(m, 2, 3);
  auto test = std::make_shared<BlockMatrix>(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size(0) == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 0)) == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 1)) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(0, 2)) == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 0)) == *E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 1)) == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(test->block(1, 2)) == *G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test->tabCol() == tCol, true);
  std::cout << "--> Constructor 1(copy) test ended with success." << std::endl;
}

void BlockMatrixTest::testConstructor2()  // Copy constructor, from a SiconosMatrix(Block)
{
  std::cout << "--> Test: constructor 2." << std::endl;
  auto ref = std::make_shared<BlockMatrix>(m, 2, 3);
  auto test = std::make_shared<BlockMatrix>(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size(0) == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size(1) == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 0)) == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 1)) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(0, 2)) == *D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 0)) == *E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 1)) == *F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *(test->block(1, 2)) == *G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test->tabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test->tabCol() == tCol, true);

  std::cout << "--> Constructor 2(copy) test ended with success." << std::endl;
}
void BlockMatrixTest::testConstructor3()  // Copy constructor, from a SiconosMatrix(Simple)
{
  std::cout << "--> Test: constructor 3." << std::endl;
  auto ref = std::make_shared<siconos::algebra::SimpleMatrix>(5, 7, 2.3);
  auto test = std::make_shared<BlockMatrix>(*ref);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size(0) == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size(1) == 7, true);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", fabs((*test)(i, j) - 2.3) < tol,
                                   true);
  auto tab = test->tabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 5, true);

  tab = test->tabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 7, true);

  std::cout << "--> Constructor 3(copy) test ended with success." << std::endl;
}

void BlockMatrixTest::testConstructor4()  // Constructor from 4 std::shared_ptr<SiconosMatrix>
{
  std::cout << "--> Test: constructor 4." << std::endl;
  auto test = std::make_shared<BlockMatrix>(B, C, E, F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size(0) == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size(1) == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->block(1, 1) == F, true);
  auto tab = test->tabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[1] == 5, true);
  tab = test->tabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tab->size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", (*tab)[1] == 6, true);
  std::cout << "--> Constructor 4 test ended with success." << std::endl;
}

// void BlockMatrixTest::testResize()
// {
//   std::cout << "--> Test: resize." <<std::endl;
//   auto test= std::make_shared<BlockMatrix>(m,2,3);
//   test->resize(3,4);
//   (*test)(2,0) = B;
//   (*test)(2,1) = C;
//   (*test)(2,2) = D;
//   (*test)(2,3) = C;
//   (*test)(0,3) = C;
//   (*test)(1,3) = F;

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->isBlock() == true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size(0) == 7, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size(1) == 11, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(0,3) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,0) == E, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,1) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,2) == G, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(1,3) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->block(2,3) == C, true);

//   std::cout << "--> resize test ended with success." <<std::endl;
// }

void BlockMatrixTest::testNormInf() {
  std::cout << "--> Test: normInf." << std::endl;
  auto test = std::make_shared<BlockMatrix>(m, 2, 3);
  test->zero();
  double n = 12;
  (*test)(4, 3) = n;
  (*test)(2, 1) = n - 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", test->normInf() == n, true);
  std::cout << "--> normInf test ended with success." << std::endl;
}

void BlockMatrixTest::testZero() {
  std::cout << "--> Test: zero." << std::endl;
  auto A = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  A->eye();
  auto H = std::make_shared<siconos::algebra::SimpleMatrix>(2, 4);
  H->eye();
  auto I = std::make_shared<siconos::algebra::SimpleMatrix>(5, 2);
  I->eye();
  auto J = std::make_shared<siconos::algebra::SimpleMatrix>(5, 4);
  J->eye();

  std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>> v(4);
  v[0] = A;
  v[1] = H;
  v[2] = I;
  v[3] = J;
  auto test = std::make_shared<BlockMatrix>(v, 2, 2);
  test->zero();
  unsigned int n1 = test->size(0);
  unsigned int n2 = test->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*test)(i, j) == 0, true);
  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*A)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*J)(i, j) == 0, true);
  }

  std::cout << "--> zero test ended with success." << std::endl;
}

void BlockMatrixTest::testEye() {
  std::cout << "--> Test: eye." << std::endl;
  auto A = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  auto H = std::make_shared<siconos::algebra::SimpleMatrix>(2, 4);
  auto I = std::make_shared<siconos::algebra::SimpleMatrix>(5, 2);
  auto J = std::make_shared<siconos::algebra::SimpleMatrix>(5, 4);

  std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>> v(4);
  v[0] = A;
  v[1] = H;
  v[2] = I;
  v[3] = J;
  auto test = std::make_shared<BlockMatrix>(v, 2, 2);
  test->eye();
  unsigned int n1 = test->size(0);
  unsigned int n2 = test->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 0, true);

  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 0, true);
    }
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i) {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j) {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 0, true);
    }
  }

  std::cout << "--> eye test ended with success." << std::endl;
}
// Add tests with getDense ...

void BlockMatrixTest::testGetSetRowCol() {
  std::cout << "--> Test: get, set Row and Col." << std::endl;
  auto tmp = std::make_shared<siconos::algebra::SiconosVector>(6);
  auto tmp1 = std::make_shared<siconos::algebra::SiconosVector>(6);
  (*tmp1)(0) = 1;
  (*tmp1)(2) = 2;
  auto A = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  auto H = std::make_shared<siconos::algebra::SimpleMatrix>(2, 4);
  auto I = std::make_shared<siconos::algebra::SimpleMatrix>(5, 2);
  auto J = std::make_shared<siconos::algebra::SimpleMatrix>(5, 4);
  std::vector<std::shared_ptr<siconos::algebra::SiconosMatrix>> v(4);
  v[0] = A;
  v[1] = H;
  v[2] = I;
  v[3] = J;
  auto test = std::make_shared<BlockMatrix>(v, 2, 2);

  test->setRow(1, *tmp1);
  test->getRow(1, *tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp == *tmp1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(1, 0) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*H)(1, 0) == 2, true);

  auto tmp2 = std::make_shared<siconos::algebra::SiconosVector>(7);
  auto tmp3 = std::make_shared<siconos::algebra::SiconosVector>(7);
  (*tmp3)(0) = 1;
  (*tmp3)(2) = 2;
  test->setCol(1, *tmp3);
  test->getCol(1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp2 == *tmp3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(0, 1) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*I)(0, 1) == 2, true);

  std::cout << "--> get, set Row and Col tests ended with success." << std::endl;
}

void BlockMatrixTest::testAssignment() {
  std::cout << "--> Test: assignment." << std::endl;
  auto Btmp = std::make_shared<siconos::algebra::SimpleMatrix>(2, 2);
  auto Ctmp = std::make_shared<siconos::algebra::SimpleMatrix>(2, 5);
  auto Dtmp = std::make_shared<siconos::algebra::SimpleMatrix>(3, 2);
  auto Etmp = std::make_shared<siconos::algebra::SimpleMatrix>(3, 5);

  auto test = std::make_shared<BlockMatrix>(Btmp, Ctmp, Dtmp, Etmp);
  // Block = Siconos(Simple)
  unsigned int size0 = test->size(0), size1 = test->size(1);
  auto ref = std::make_shared<siconos::algebra::SimpleMatrix>(size0, size1);
  for (unsigned int i = 0; i < size0; ++i)
    for (unsigned int j = 0; j < size1; ++j) (*ref)(i, j) = i + j;
  *test = *ref;
  for (unsigned int i = 0; i < size0; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testAssignment: ", fabs((*test)(i, j) - (*ref)(i, j)) < tol, true);

  // Block = Block
  test->zero();
  auto ref2 = std::make_shared<BlockMatrix>(m, 2, 3);
  *test = *ref2;
  for (unsigned int i = 0; i < size0; ++i)
    for (unsigned int j = 0; j < size1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testAssignment: ", fabs((*test)(i, j) - (*ref2)(i, j)) < tol, true);
  std::cout << "-->  test assignment ended with success." << std::endl;
}

void BlockMatrixTest::testOperators1() {
  std::cout << "--> Test: operators1." << std::endl;
  double tol = 1e-10;
  auto Ab = std::make_shared<BlockMatrix>(m, 2, 3);
  auto Cb = std::make_shared<BlockMatrix>(*Ab);
  auto A = std::make_shared<siconos::algebra::SimpleMatrix>(5, 7);

  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j) (*A)(i, j) = i + j;

  double a = 2.3;
  int a1 = 2;

  // Block *= scal or /= scal
  *Cb *= a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators1: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);

  *Cb *= a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators1: ", fabs((*Cb)(i, j) - a1 * a * (*Ab)(i, j)) < tol, true);

  *Cb /= a;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators1: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol, true);
  *Cb /= a1;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j)) < tol,
                                   true);

  // Block +=  Simple
  *Cb += *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);

  // Block -=  Block
  *Cb -= *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*A)(i, j)) < tol,
                                   true);

  // Block += Block
  *Cb += *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);

  // Block -= Simple
  *Cb -= *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators1: ", fabs((*Cb)(i, j) - (*Ab)(i, j)) < tol,
                                   true);

  std::cout << "-->  test operators1 ended with success." << std::endl;
}

void BlockMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of BlockMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
