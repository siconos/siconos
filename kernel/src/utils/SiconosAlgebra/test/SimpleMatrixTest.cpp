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

  // vect1 = std::make_shared<siconos::algebra::SiconosVector>(v3);
  // vect2 = std::make_shared<siconos::algebra::SiconosVector>(
  //     v4);  // vect2 != vect1, but vect2 == SimM second column
  // vect3 = std::make_shared<siconos::algebra::SiconosVector>(
  //     v5);  // vect3 != vect1, but vect3 == SimM second row

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
  double n = SicM->normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7, true);
  std::cout << "--> normInf test ended with success." << std::endl;
}

void SimpleMatrixTest::testSetBlock() {
  std::cout << "--> Test: testSetBlock." << std::endl;

  // Copy of a sub-block of a Simple into a Simple
  auto MIn = std::make_shared<SiconosMatrix>(10, 10);
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 10; ++j) (*MIn)(i, j) = i + j;

  auto MOut = std::make_shared<SiconosMatrix>(5, 5);

  std::vector<std::size_t> subDim(2);
  std::vector<std::size_t> subPos(4);
  subDim[0] = 2;
  subDim[1] = 3;
  subPos[0] = 1;
  subPos[1] = 2;
  subPos[2] = 1;
  subPos[3] = 2;

  siconos::algebra::setBlock(*MIn, MOut, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3]; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*MOut)(i, j) - (*MIn)(i, j)) < tol,
                                   true);

  // Copy of a sub-block of a Simple into a Block
  Cb->setZero();
  siconos::algebra::setBlock(*MIn, Cb, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3]; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*Cb)(i, j) - (*MIn)(i, j)) < tol,
                                   true);

  // Copy of a sub-block of a Block into a Simple

  MOut = std::make_shared<SiconosMatrix>(5, 5);
  siconos::algebra::setBlock(*Ab, MOut, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3]; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*MOut)(i, j) - (*Ab)(i, j)) < tol,
                                   true);

  std::cout << "-->  setBlock test ended with success." << std::endl;
}

void SimpleMatrixTest::testSetBlock2() {
  std::cout << "--> Test: testSetBlock2." << std::endl;
  // Copy of a Simple into a sub-block of Simple
  auto MOut = std::make_shared<SiconosMatrix>(10, 10);

  auto MIn = std::make_shared<SiconosMatrix>(5, 5);
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 5; ++j) (*MIn)(i, j) = i + j;

  MOut->setBlock(2, 3, *MIn);

  for (unsigned int i = 2; i < 7; ++i)
    for (unsigned int j = 3; j < 8; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testSetBlock2: ", fabs((*MOut)(i, j) - (*MIn)(i - 2, j - 3)) < tol, true);

  // Copy of a Block into a sub-block of Simple

  // auto MInB = std::make_shared<siconos::algebra::BlockMatrix>(m4, m4, m4, m4);
  // MOut->setBlock(2, 3, *MInB);

  // for (unsigned int i = 2; i < 6; ++i)
  //   for (unsigned int j = 3; j < 7; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testSetBlock2: ", fabs((*MOut)(i, j) - (*MInB)(i - 2, j - 3)) < tol, true);

  std::cout << "-->  setBlock2 test ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators9() {
  std::cout << "--> Test: operator9." << std::endl;

  // C = a*A or A/a

  double a = 2.2;
  int a1 = 3;

  // Simple = a * Simple or Simple/a
  *C = a * *A;
  for (unsigned int i = 0; i < C->rows(); ++i)
    for (unsigned int j = i; j < C->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol,
                                   true);
  *C = a1 * *A;
  for (unsigned int i = 0; i < C->rows(); ++i)
    for (unsigned int j = i; j < C->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a1 * (*A)(i, j)) < tol, true);

  // *C = *A / a;
  // for (unsigned int i = 0; i < C->rows(); ++i)
  //   for (unsigned int j = i ; j < C->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a) <
  //     tol, true);
  // *C = *A / a1;
  // for (unsigned int i = 0; i < C->rows(); ++i)
  //   for (unsigned int j = i ; j < C->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a1) <
  //     tol, true);

  // Simple = a * Block

  *C = a * *Ab;
  for (unsigned int i = 0; i < C->rows(); ++i)
    for (unsigned int j = i; j < C->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);
  ;
  *C = a1 * *Ab;
  for (unsigned int i = 0; i < C->rows(); ++i)
    for (unsigned int j = i; j < C->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  // *C = *Ab / a;
  // for (unsigned int i = 0; i < C->rows(); ++i)
  //   for (unsigned int j = i ; j < C->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a) <
  //     tol, true);
  // *C = *Ab / a1;
  // for (unsigned int i = 0; i < C->rows(); ++i)
  //   for (unsigned int j = i ; j < C->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a1) <
  //     tol, true);

  // Block = a * Block
  *Cb = a * *Ab;
  for (unsigned int i = 0; i < Cb->rows(); ++i)
    for (unsigned int j = i; j < Cb->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *Cb = a1 * *Ab;
  for (unsigned int i = 0; i < Cb->rows(); ++i)
    for (unsigned int j = i; j < Cb->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  // *Cb = *Ab / a;
  // for (unsigned int i = 0; i < Cb->rows(); ++i)
  //   for (unsigned int j = i ; j < Cb->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a) <
  //     tol, true);
  // *Cb = *Ab / a1;
  // for (unsigned int i = 0; i < Cb->rows(); ++i)
  //   for (unsigned int j = i ; j < Cb->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a1)
  //     < tol, true);

  // Block = a * Simple
  *Cb = a * *A;
  for (unsigned int i = 0; i < Cb->rows(); ++i)
    for (unsigned int j = i; j < Cb->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);
  *Cb = a1 * *A;
  for (unsigned int i = 0; i < Cb->rows(); ++i)
    for (unsigned int j = i; j < Cb->cols(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a1 * (*A)(i, j)) < tol, true);

  // *Cb = *A / a;
  // for (unsigned int i = 0; i < Cb->rows(); ++i)
  //   for (unsigned int j = i ; j < Cb->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a) <
  //     tol, true);
  // *Cb = *A / a1;
  // for (unsigned int i = 0; i < Cb->rows(); ++i)
  //   for (unsigned int j = i ; j < Cb->cols(); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a1) <
  //     tol, true);
  std::cout << "-->  test operators9 ended with success." << std::endl;
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

  // Block = Simple * Block
  yB->setZero();
  //  siconos::algebra::prod(*A ,*xB,*yB,false);
  //  siconos::algebra::prod(*A ,*xB,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->cols(); ++j)
  //      sum += 2*(*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
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
  // Simple = Simple * Block
  y->setZero();
  //  siconos::algebra::prod(a,*A ,*xB,*y,false);
  //  siconos::algebra::prod(a,*A ,*xB,*y,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->cols(); ++j)
  //      sum += a*2*(*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i) - sum)< tol, true);
  //  }

  // Block = Simple * Simple
  yB->setZero();
  //  siconos::algebra::prod(a,*A ,*x,*yB,false);
  //  siconos::algebra::prod(a,*A ,*x,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->cols(); ++j)
  //      sum += a*2*(*A)(i,j)*(*x)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  //
  // Block = Simple * Block
  yB->setZero();
  //  siconos::algebra::prod(a,*A ,*xB,*yB,false);
  //  siconos::algebra::prod(a,*A ,*xB,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->cols(); ++j)
  //      sum += a*2*(*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
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
  // Simple = Simple * Block
  y->setZero();
  //  siconos::algebra::prod(*xB,*A,*y);
  //  siconos::algebra::prod(*xB,*A,*y,false);
  //
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< size; ++j)
  //      sum += 2*(*tmp)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*y)(i) - sum)< tol, true);
  //  }

  // Block = Simple * Simple
  yB->setZero();
  siconos::algebra::transposeMatrixVector_prod_toBlock(*x, *A, *yB);
  siconos::algebra::transposeMatrixVector_prod_toBlock(*x, *A, *yB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i)-sum) < tol, true);
  }

  // Block = Simple * Block
  yB->setZero();
  //  siconos::algebra::prod(*xB,*A ,*yB);
  //  siconos::algebra::prod(*xB,*A ,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->cols(); ++j)
  //      sum += 2*(*tmp)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  std::cout << "-->  test ublas::prod6 ended with success." << std::endl;
}

void SimpleMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
