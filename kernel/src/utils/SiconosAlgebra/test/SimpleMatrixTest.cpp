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
#include "SimpleMatrixTest.hpp"

#include "NumericsToolsNamespace.h"  // for NM_csc, NM_free ...
#include "SiconosAlgebraTypes.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"        // for prod
#include "SiconosMatrixVectorOp.hpp"  // for prod
#include "Tools.hpp"
// namespace ublas = boost::numeric::ublas;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

CPPUNIT_TEST_SUITE_REGISTRATION(SimpleMatrixTest);

#define DEBUG_MESSAGES
#include "siconos_debug.h"

// Note FP: tests are (rather) complete for Dense objects but many are missing for other cases
// (Triang, Symm etc ...).

using SiconosMatrix = siconos::algebra::SiconosMatrix;

void SimpleMatrixTest::setUp() {
  tol = 1e-9;

  fic1 = "mat1.dat";  // 2 X 2
  fic2 = "mat2.dat";  // 2 X 3
  SicM = std::make_shared<SiconosMatrix>(siconos::algebra::readMatrixFromFile(fic1));
  SimM = std::make_shared<SiconosMatrix>(siconos::algebra::readMatrixFromFile(fic2));

  std::vector<double> v3(2, 0);
  std::vector<double> v4(2, 0);
  std::vector<double> v5(3, 0);
  v4[0] = 6;
  v4[1] = 9;
  v5[0] = 8;
  v5[1] = 9;
  v5[2] = 10;
  // Dense

  for (unsigned i = 0; i < D.rows(); ++i)
    for (unsigned j = 0; j < D.cols(); ++j) D(i, j) = 3 * i + j;

  // BlockMat
  size = 10;
  size2 = 10;

  C = std::make_shared<SiconosMatrix>(size, size);
  A = std::make_shared<SiconosMatrix>(siconos::algebra::readMatrixFromFile("A.dat"));
  B = std::make_shared<SiconosMatrix>(siconos::algebra::readMatrixFromFile("B.dat"));

  m1 = std::make_shared<SiconosMatrix>(size - 2, size - 2);
  m2 = std::make_shared<SiconosMatrix>(size - 2, 2);
  m3 = std::make_shared<SiconosMatrix>(2, size - 2);
  m4 = std::make_shared<SiconosMatrix>(2, 2);
  m5 = std::make_shared<SiconosMatrix>(size2 - 2, size2 - 2);
  m6 = std::make_shared<SiconosMatrix>(size2 - 2, 2);
  m7 = std::make_shared<SiconosMatrix>(2, size2 - 2);
  m8 = std::make_shared<SiconosMatrix>(2, 2);
  Ab = (std::make_shared<siconos::algebra::BlockMatrix>(m1, m2, m3, m4))->toSiconosMatrix();
  Cb = std::make_shared<siconos::algebra::BlockMatrix>(m5, m6, m7, m8)->toSiconosMatrix();
}

void SimpleMatrixTest::tearDown() {}

//______________________________________________________________________________

void SimpleMatrixTest::testNormInf() {
  std::cout << "--> Test: normInf." << std::endl;
  double n = siconos::algebra::normInf(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7, true);
  std::cout << "--> normInf test ended with success." << std::endl;
}

void SimpleMatrixTest::testSetBlock() {
  std::cout << "--> Test: testSetBlock." << std::endl;

  // Copy of a sub-block of a Simple into a Simple
  SiconosMatrix result{10, 8};
  result.setRandom();
  SiconosMatrix mIn{2, 3};
  mIn << 1, 2, 3, 4, 5, 6;
  siconos::algebra::setBlock(1, 2, mIn, result);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", (result.block(1, 2, 2, 3) - mIn).norm() < tol,
                               true);
  std::cout << "-->  setBlock test ended with success." << std::endl;
}

void SimpleMatrixTest::testProdBis() {
  std::cout << "--> Test: ublas::prod. mat-vect (bis)" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->setZero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  siconos::algebra::prod(*A, *x, *y);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  siconos::algebra::matrixBlockVector_prod(*A, *xB, *y);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i)-sum) < tol, true);
  }

  // Block = Simple * Simple
  siconos::algebra::matrixVector_prod_toBlock(*A, *x, *yB);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*yB)(i)-sum) < tol, true);
  }
  std::cout << "-->  test ublas::prodBis ended with success." << std::endl;
}
void SimpleMatrixTest::testProd4()  // y += A*x
{
  std::cout << "--> Test: ublas::prod. mat-vect (4)" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->setZero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  y->setZero();
  siconos::algebra::prod(*A, *x, *y, false);
  siconos::algebra::prod(*A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  y->setZero();
  siconos::algebra::matrixBlockVector_prod(*A, *xB, *y, false);
  siconos::algebra::matrixBlockVector_prod(*A, *xB, *y, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i)-sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->setZero();
  siconos::algebra::matrixVector_prod_toBlock(*A, *x, *yB, false);
  siconos::algebra::matrixVector_prod_toBlock(*A, *x, *yB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*yB)(i)-sum) < tol, true);
  }

  std::cout << "-->  test ublas::prod4 ended with success." << std::endl;
}

void SimpleMatrixTest::testProd5()  // y += a*A*x
{
  std::cout << "--> Test: ublas::prod. mat-vect (5)" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->setZero();

  // Matrix - vector ublas::product
  double a = 3.0;
  // Simple = Simple * Simple
  y->setZero();
  siconos::algebra::prod(a, *A, *x, *y, false);
  siconos::algebra::prod(a, *A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * a * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i)-sum) < tol, true);
  }
  std::cout << "-->  test ublas::prod5 ended with success." << std::endl;
}

void SimpleMatrixTest::testProd6()  // y += trans(A)*x
{
  std::cout << "--> Test: ublas::prod. mat-vect (6)" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->setZero();

  auto tmp = std::make_shared<SiconosMatrix>(*A);
  tmp->transposeInPlace();
  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  y->setZero();
  siconos::algebra::transposeMatrixVector_prod(*x, *A, *y);
  siconos::algebra::transposeMatrixVector_prod(*x, *A, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*y)(i)-sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->setZero();
  siconos::algebra::transposeMatrixVector_prod_toBlock(*x, *A, *yB);
  siconos::algebra::transposeMatrixVector_prod_toBlock(*x, *A, *yB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i)-sum) < tol, true);
  }

  std::cout << "-->  test ublas::prod6 ended with success." << std::endl;
}

void SimpleMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
