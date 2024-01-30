/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

// #include <boost/numeric/bindings/ublas/matrix.hpp>
// #include <boost/numeric/ublas/banded.hpp>
// #include <boost/numeric/ublas/matrix_sparse.hpp>
// #include <boost/numeric/ublas/symmetric.hpp>
// #include <boost/numeric/ublas/triangular.hpp>

#include "NumericsToolsNamespace.h"  // for NM_csc, NM_free ...
#include "SiconosAlgebraTypes.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosMatrixOp.hpp"        // for prod
#include "SiconosMatrixVectorOp.hpp"  // for prod
#include "SimpleMatrix.hpp"
#include "Tools.hpp"
// namespace ublas = boost::numeric::ublas;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

CPPUNIT_TEST_SUITE_REGISTRATION(SimpleMatrixTest);

#define DEBUG_MESSAGES
#include "siconos_debug.h"

// Note FP: tests are (rather) complete for Dense objects but many are missing for other cases
// (Triang, Symm etc ...).

using SimpleMatrix = siconos::algebra::SimpleMatrix;

void SimpleMatrixTest::setUp() {
  tol = 1e-9;

  fic1 = "mat1.dat";  // 2 X 2
  fic2 = "mat2.dat";  // 2 X 3
  SicM = std::make_shared<SimpleMatrix>(siconos::algebra::readMatrixFromFile(fic1));
  SimM = std::make_shared<SimpleMatrix>(siconos::algebra::readMatrixFromFile(fic2));

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
  D = std::make_shared<siconos::algebra::SiconosMatrix>(2, 2);
  for (unsigned i = 0; i < D->size(0); ++i)
    for (unsigned j = 0; j < D->size(1); ++j) (*D)(i, j) = 3 * i + j;

  // // Triang
  // T = std::make_shared<siconos::algebra::TriangMat>(3, 3);
  // for (unsigned i = 0; i < T->size1(); ++i)
  //   for (unsigned j = i; j < T->size2(); ++j) (*T)(i, j) = 3 * i + j;
  // T2 = std::make_shared<siconos::algebra::TriangMat>(4, 4);
  // for (unsigned i = 0; i < T2->size1(); ++i)
  //   for (unsigned j = i; j < T2->size2(); ++j) (*T2)(i, j) = 3 * i + j;

  // // Sym
  // S = std::make_shared<siconos::algebra::SymMat>(3, 3);
  // for (unsigned i = 0; i < S->size1(); ++i)
  //   for (unsigned j = i; j < S->size2(); ++j) (*S)(i, j) = 3 * i + j;
  // S2 = std::make_shared<siconos::algebra::SymMat>(4, 4);
  // for (unsigned i = 0; i < S2->size1(); ++i)
  //   for (unsigned j = i; j < S2->size2(); ++j) (*S2)(i, j) = 3 * i + j;

  // // Sparse
  // SP = std::make_shared<siconos::algebra::SparseMat>(4, 4);
  // for (unsigned i = 0; i < SP->size1(); ++i)
  //   for (unsigned j = 0; j < SP->size2(); ++j) (*SP)(i, j) = 3 * i + j;

  // SP2 = std::make_shared<siconos::algebra::SparseMat>(4, 4);
  // for (unsigned i = 0; i < SP2->size1(); ++i)
  //   for (unsigned j = 0; j < SP->size2() - 1; ++j)
  //     if (i != j) (*SP2)(i, j) = 3 * i + j;

  // SP3 = std::make_shared<siconos::algebra::SparseMat>(3, 3);

  // (*SP3)(0, 0) = 1.0;
  // (*SP3)(0, 2) = 2.0;
  // (*SP3)(1, 2) = 3.0;
  // (*SP3)(2, 0) = 4.0;
  // (*SP3)(2, 1) = 5.0;
  // (*SP3)(2, 2) = 5.0;

  // SP4 = std::make_shared<siconos::algebra::SparseMat>(3, 3);
  // (*SP4)(0, 0) = 1.0;
  // (*SP4)(0, 2) = 4.0;
  // (*SP4)(1, 1) = 1.0;
  // (*SP4)(1, 2) = 5.0;
  // (*SP4)(2, 0) = 4.0;
  // (*SP4)(2, 1) = 5.0;
  // (*SP4)(2, 2) = 77.0;

  // Sparse Coordinate
  // SP_coor = std::make_shared<siconos::algebra::SparseCoordinateMat>(4, 4);
  // for (unsigned i = 0; i < SP->size1(); ++i)
  //   for (unsigned j = 0; j < SP->size2(); ++j) (*SP_coor)(i, j) = 3 * i + j;

  // // Banded
  // Band = std::make_shared<siconos::algebra::BandedMat>(4, 4, 1, 1);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     (*Band)(i, j) = 3 * i + j;
  // Band2 = std::make_shared<siconos::algebra::BandedMat>(4, 3, 1, 1);
  // for (signed i = 0; i < signed(Band2->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band2->size2())); ++j)
  //     (*Band2)(i, j) = 3 * i + j;

  // // Zero
  // Z = std::make_shared<siconos::algebra::ZeroMat>(3, 3);
  // Z2 = std::make_shared<siconos::algebra::ZeroMat>(4, 4);
  // // Identity
  // I = std::make_shared<siconos::algebra::IdentityMat>(3, 3);
  // I2 = std::make_shared<siconos::algebra::IdentityMat>(4, 4);

  // BlockMat
  size = 10;
  size2 = 10;

  C = std::make_shared<SimpleMatrix>(size, size);
  A = std::make_shared<SimpleMatrix>(siconos::algebra::readMatrixFromFile("A.dat"));
  B = std::make_shared<SimpleMatrix>(siconos::algebra::readMatrixFromFile("B.dat"));

  m1 = std::make_shared<SimpleMatrix>(size - 2, size - 2);
  m2 = std::make_shared<SimpleMatrix>(size - 2, 2);
  m3 = std::make_shared<SimpleMatrix>(2, size - 2);
  m4 = std::make_shared<SimpleMatrix>(2, 2);
  m5 = std::make_shared<SimpleMatrix>(size2 - 2, size2 - 2);
  m6 = std::make_shared<SimpleMatrix>(size2 - 2, 2);
  m7 = std::make_shared<SimpleMatrix>(2, size2 - 2);
  m8 = std::make_shared<SimpleMatrix>(2, 2);
  Ab = (std::make_shared<siconos::algebra::BlockMatrix>(m1, m2, m3, m4))->toSiconosMatrix();
  Bb = (std::make_shared<siconos::algebra::BlockMatrix>(3 * *Ab)->toSiconosMatrix());
  Cb = std::make_shared<siconos::algebra::BlockMatrix>(m5, m6, m7, m8)->toSiconosMatrix();
}

void SimpleMatrixTest::tearDown() {}

//______________________________________________________________________________

void SimpleMatrixTest::testConstructor0()  // constructor with TYP and dim
{
  std::cout << "====================================" << std::endl;
  std::cout << "=== Simple Matrix tests start ...=== " << std::endl;
  std::cout << "====================================" << std::endl;
  std::cout << "--> Test: constructor 0." << std::endl;
  auto test = std::make_shared<SimpleMatrix>(2, 3);

  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testConstructor0 : ", test->num() == siconos::algebra::UblasType::DENSE, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(0) == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size(1) == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->normInf() < tol, true);
  std::cout << "--> Constructor 0 test ended with success." << std::endl;
}

void SimpleMatrixTest::testConstructor1()  // Copy constructor, from a SimpleMatrix
{
  std::cout << "--> Test: constructor 1." << std::endl;
  auto test = std::make_shared<SimpleMatrix>(*SimM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test == *SimM, true);
  std::cout << "--> Constructor 1 (copy) test ended with success." << std::endl;
}

void SimpleMatrixTest::testConstructor2()  // Copy constructor, from a SiconosMatrix
{
  std::cout << "--> Test: constructor 2." << std::endl;
  auto test = std::make_shared<SimpleMatrix>(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test == *SicM, true);
  std::cout << "--> Constructor 2 (copy) test ended with success." << std::endl;
}

void SimpleMatrixTest::testConstructor3()  // Copy constructor, from a BlockMatrix
{
  std::cout << "--> Test: constructor 3." << std::endl;
  auto test = std::make_shared<SimpleMatrix>(*Ab);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", *test == *Ab, true);
  std::cout << "--> Constructor 3 (copy) test ended with success." << std::endl;
}

void SimpleMatrixTest::testConstructor4() {
  std::cout << "--> Test: constructor 4." << std::endl;
  auto test = std::make_shared<SimpleMatrix>(*D);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
      // "testConstructor4 : ", test->num() == siconos::algebra::UblasType::DENSE, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", (*test - *D).normInf() == 0,   // SAM : TODO
  //                              true);
  std::cout << "--> Constructor 4 test ended with success." << std::endl;
}

// void SimpleMatrixTest::testConstructor5() {
//   std::cout << "--> Test: constructor 5." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*T);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor5 : ", test->num() == siconos::algebra::UblasType::TRIANGULAR, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", norm_inf(test->getTriang() - *T) == 0,
//                                true);
//   std::cout << "--> Constructor 5 test ended with success." << std::endl;
// }

// void SimpleMatrixTest::testConstructor6() {
//   std::cout << "--> Test: constructor 6." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*S);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor6 : ", test->num() == siconos::algebra::UblasType::SYMMETRIC, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(test->getSym() - *S) == 0,
//                                true);
//   std::cout << "--> Constructor 6 test ended with success." << std::endl;
// }

// void SimpleMatrixTest::testConstructor7() {
//   std::cout << "--> Test: constructor 7." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*SP);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor7 : ", test->num() == siconos::algebra::UblasType::SPARSE, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", norm_inf(test->getSparse() - *SP) == 0,
//                                true);
//   std::cout << "--> Constructor 7 test ended with success." << std::endl;
// }

// void SimpleMatrixTest::testConstructor8() {
//   std::cout << "--> Test: constructor 8." << std::endl;
//   std::cout << "--> Constructor 8 test ended with success." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*Band);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor8 : ", test->num() == siconos::algebra::UblasType::BANDED, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor8 : ", norm_inf(test->getBanded() - *Band) == 0,
//                                true);
// }

// void SimpleMatrixTest::testConstructor9()  // constructor with TYP and dim and input value
// {
//   std::cout << "--> Test: constructor 9." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(2, 3, 4.5);

//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor9 : ", test->num() == siconos::algebra::UblasType::DENSE, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test->size(0) == 2, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test->size(1) == 3, true);
//   for (unsigned int i = 0; i < 2; ++i)
//     for (unsigned int j = 0; j < 3; ++j) {
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", (*test)(i, j) == 4.5, true);
//     }
//   std::cout << "--> Constructor 9 test ended with success." << std::endl;
// }

// void SimpleMatrixTest::testConstructor10() {
//   std::cout << "--> Test: constructor 10." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(fic1);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor10 : ", *test == *SicM, true);
//   std::cout << "--> Constructor 10 test ended with success." << std::endl;
// }

// void SimpleMatrixTest::testConstructor11() {
//   std::cout << "--> Test: constructor 11." << std::endl;
//   std::cout << "--> Constructor 11 test ended with success." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*Z);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor11 : ", test->num() == siconos::algebra::UblasType::ZERO, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor11 : ", test->normInf() == 0, true);
// }

// void SimpleMatrixTest::testConstructor12() {
//   std::cout << "--> Test: constructor 12." << std::endl;
//   std::cout << "--> Constructor 12 test ended with success." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(*I);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor12 : ", test->num() == siconos::algebra::UblasType::IDENTITY, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor12 : ", test->normInf() == 1, true);
// }

// void SimpleMatrixTest::testConstructor13() {
//   std::cout << "--> Test: constructor 13." << std::endl;
//   auto test = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor3 : ", test->num() == siconos::algebra::UblasType::SPARSE, true);
//   std::cout << "--> Constructor 13 test ended with success." << std::endl;
// }
// void SimpleMatrixTest::testConstructor14() {
//   std::cout << "--> Test: constructor 14." << std::endl;
//   auto test =
//       std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE_COORDINATE);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testConstructor14 : ", test->num() == siconos::algebra::UblasType::SPARSE_COORDINATE,
//       true);
//   std::cout << "--> Constructor 14 test ended with success." << std::endl;
// }
// Add tests with getDense ...

// Add tests with getDense ...

void SimpleMatrixTest::testZero() {
  std::cout << "--> Test: zero." << std::endl;
  auto tmp = std::make_shared<SimpleMatrix>(*SimM);
  tmp->zero();
  unsigned int n1 = tmp->size(0);
  unsigned int n2 = tmp->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*tmp)(i, j) == 0, true);
  std::cout << "--> zero test ended with success." << std::endl;
}

void SimpleMatrixTest::testEye() {
  std::cout << "--> Test: eye." << std::endl;
  auto tmp = std::make_shared<SimpleMatrix>(*SimM);
  tmp->eye();
  unsigned int n1 = tmp->size(0);
  unsigned int n2 = tmp->size(1);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i != j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 0, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 1, true);
  std::cout << "--> eye test ended with success." << std::endl;
}

void SimpleMatrixTest::testResize() {
  std::cout << "--> Test: resize." << std::endl;
  auto tmp = std::make_shared<SimpleMatrix>(*SicM);
  tmp->conservativeResize(3, 4);
  unsigned int n1 = SicM->size(0);
  unsigned int n2 = SicM->size(1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 4, true);

  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i, j) - (*SicM)(i, j)) < tol,
                                   true);
  //   for(unsigned int i = n1; i<3; ++i)
  //     for(unsigned int j=0;j<4;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i,j)) < tol, true);
  //   for(unsigned int j = n2; j<4; ++j)
  //     for(unsigned int i=0;i<3;++i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", fabs((*tmp)(i,j)) < tol, true)
  ;
  // Check the effect of bool = false (ie preserve == false in boost resize)
  //   tmp->resize(6,8, false);
  //   tmp->display();
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 6, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 8, true);
  //   for(unsigned int i = 0; i<6; ++i)
  //     for(unsigned int j=0;j<8;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  //   // Reduction ...
  //   tmp->resize(1,2, false);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(0) == 1, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size(1) == 2, true);
  //   for(unsigned int i = 0; i<1; ++i)
  //     for(unsigned int j=0;j<2;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  std::cout << "--> resize test ended with success." << std::endl;
}

void SimpleMatrixTest::testNormInf() {
  std::cout << "--> Test: normInf." << std::endl;
  double n = SicM->normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7, true);
  std::cout << "--> normInf test ended with success." << std::endl;
}

void SimpleMatrixTest::testSetBlock() {
  std::cout << "--> Test: testSetBlock." << std::endl;

  // Copy of a sub-block of a Simple into a Simple
  auto MIn = std::make_shared<SimpleMatrix>(10, 10);
  for (unsigned int i = 0; i < 10; ++i)
    for (unsigned int j = 0; j < 10; ++j) (*MIn)(i, j) = i + j;

  auto MOut = std::make_shared<SimpleMatrix>(5, 5);

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
  Cb->zero();
  siconos::algebra::setBlock(*MIn, Cb, subDim, subPos);

  for (unsigned int i = subPos[2]; i < subPos[2] + subDim[0]; ++i)
    for (unsigned int j = subPos[3]; j < subPos[3] + subDim[1]; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", fabs((*Cb)(i, j) - (*MIn)(i, j)) < tol,
                                   true);

  // Copy of a sub-block of a Block into a Simple

  MOut = std::make_shared<SimpleMatrix>(5, 5);
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
  auto MOut = std::make_shared<SimpleMatrix>(10, 10);

  auto MIn = std::make_shared<SimpleMatrix>(5, 5);
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

void SimpleMatrixTest::testGetSetRowCol() {
  std::cout << "--> Test: get, set Row and Col." << std::endl;

  auto vIn = std::make_shared<siconos::algebra::SiconosVector>(10);
  vIn->setConstant(1.2);
  auto vBIn = std::make_shared<siconos::algebra::BlockVector>();
  auto v1 = std::make_shared<siconos::algebra::SiconosVector>(3, 2);
  auto v2 = std::make_shared<siconos::algebra::SiconosVector>(5, 3);
  auto v3 = std::make_shared<siconos::algebra::SiconosVector>(2, 4);
  vBIn->insertPtr(v1);
  vBIn->insertPtr(v2);
  vBIn->insertPtr(v3);

  // Set row with a SiconosVector
  C->setRow(4, *vIn);
  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - 1.2) < tol, true);

  // Set col with a SiconosVector
  C->setCol(4, *vIn);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - 1.2) < tol, true);

  //  C->setCol(4, *vBIn);
  //  for (unsigned int i = 0; i< C->size(1); ++i)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4,i)- (*vBIn)(i)) < tol,
  //    true);

  *C = *A;  // reset C
  vIn->zero();
  vBIn->zero();
  // get row and copy it into a SiconosVector
  *vIn = C->row(4);
  for (unsigned int i = 0; i < C->size(1); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(4, i) - (*vIn)(i)) < tol,
                                 true);

  // get col and copy it into a SiconosVector
  *vIn = C->col(4);
  for (unsigned int i = 0; i < C->size(0); ++i)
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", fabs((*C)(i, 4) - (*vIn)(i)) < tol,
                                 true);

  std::cout << "--> get, set Row and Col tests ended with success." << std::endl;
}

void SimpleMatrixTest::testTrans() {
  std::cout << "--> Test: trans." << std::endl;

  // Transpose in place ...
  auto ref = std::make_shared<SimpleMatrix>(*D);
  auto tRef = std::make_shared<SimpleMatrix>(*ref);

  tRef->transposeInPlace();
  for (unsigned int i = 0; i < ref->size(0); ++i)
    for (unsigned int j = 0; j < ref->size(1); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j), true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i), true);

  // Transpose of another matrix ...
  // Dense
  tRef->zero();

  *tRef = ref->transpose();
  for (unsigned int i = 0; i < ref->size(0); ++i)
    for (unsigned int j = 0; j < ref->size(1); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j), true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i), true);

  // // Sym
  // ref = std::make_shared<SimpleMatrix>(*S);
  // tRef = std::make_shared<SimpleMatrix>(*ref);
  // std::cout << siconos::tools::enum_to_string(ref->num()) << " "
  //           << siconos::tools::enum_to_string(tRef->num()) << "\n";
  // tRef->trans(*ref);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef) == (*ref), true);
  // // Sparse
  // ref = std::make_shared<SimpleMatrix>(*SP);
  // tRef = std::make_shared<SimpleMatrix>(*ref);
  // tRef->trans(*ref);
  //   for(unsigned int i = 0; i<ref->size(0); ++i)
  //     {
  //       for(unsigned int j = 0 ; j< ref->size(1); ++j)
  //  if(i==j)
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(i,j) , true);
  //  else
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(j,i) , true);
  //     }
  // Banded
  //   ref = std::make_shared<SimpleMatrix>(*Band);
  //   tRef = std::make_shared<SimpleMatrix>(*ref);
  //   *tRef = trans(*ref);
  //   for(unsigned int i = 0; i<ref->size(0); ++i)
  //     for(unsigned int j = 0 ; j< ref->size(1); ++j)
  //       if(i==j)
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(i,j) , true);
  //       else
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(j,i) , true);
  std::cout << "-->  test trans ended with success." << std::endl;
}

void SimpleMatrixTest::testAssignment0() {
  std::cout << "--> Test: assignment0." << std::endl;

  // Simple = Simple

  auto ref = std::make_shared<SimpleMatrix>(*D);
  auto tRef = std::make_shared<SimpleMatrix>(*SicM);
  // Dense = any type:
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*T);
  // auto tRef3 = std::make_shared<SimpleMatrix>(3, 3);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*S);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref), true);

  // auto tRef4 = std::make_shared<SimpleMatrix>(4, 4);
  // ref = std::make_shared<SimpleMatrix>(*SP);
  // *tRef4 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef4) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Band);
  // *tRef4 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef4) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef3) == (*ref), true);

  // // Triang = Triang, Zero or Identity
  // ref = std::make_shared<SimpleMatrix>(*T);
  // tRef = std::make_shared<SimpleMatrix>(*T);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // // Sym = Sym, Zero or Id
  // ref = std::make_shared<SimpleMatrix>(*S);
  // tRef = std::make_shared<SimpleMatrix>(*S);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // // Sparse => Sparse or Zero
  // ref = std::make_shared<SimpleMatrix>(*SP);
  // tRef = std::make_shared<SimpleMatrix>(*SP);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Z2);
  // *tRef = *ref;

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // // // Sparse coordinate => Sparse
  // ref = std::make_shared<SimpleMatrix>(*SP_coor);
  // // ref->display();
  // tRef = std::make_shared<SimpleMatrix>(*SP);
  // tRef->zero();
  // *tRef = *ref;
  // // tRef->displayExpert();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  // // Banded = Banded, Id or Zero
  // ref = std::make_shared<SimpleMatrix>(*Band);
  // tRef = std::make_shared<SimpleMatrix>(*Band);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*Z2);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*I2);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment0: ", (*tRef) == (*ref), true);

  std::cout << "-->  test assignment0 ended with success." << std::endl;
}

void SimpleMatrixTest::testAssignment1() {
  std::cout << "--> Test: assignment1." << std::endl;

  // Simple = Siconos(Block)

  *C = *Ab;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment1: ", (*C) == (*Ab), true);
  std::cout << "-->  test assignment1 ended with success." << std::endl;
}

void SimpleMatrixTest::testAssignment2() {
  std::cout << "--> Test: assignment2." << std::endl;

  // Simple = Siconos(Simple)

  auto ref = std::make_shared<SimpleMatrix>(*D);
  auto tRef = std::make_shared<SimpleMatrix>(*SicM);
  // Dense = any type:
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*T);
  // auto tRef3 = std::make_shared<SimpleMatrix>(3, 3);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*S);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref), true);

  // auto tRef4 = std::make_shared<SimpleMatrix>(4, 4);
  // ref = std::make_shared<SimpleMatrix>(*SP);
  // *tRef4 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef4) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Band);
  // *tRef4 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef4) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef3 = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef3) == (*ref), true);
  // // Triang = Triang, Zero or Identity
  // ref = std::make_shared<SimpleMatrix>(*T);
  // tRef = std::make_shared<SimpleMatrix>(*T);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  // // Sym = Sym, Zero or Id
  // ref = std::make_shared<SimpleMatrix>(*S);
  // tRef = std::make_shared<SimpleMatrix>(*S);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*Z);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*I);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  // // Sparse = Sparse or Zero
  // ref = std::make_shared<SimpleMatrix>(*SP);
  // tRef = std::make_shared<SimpleMatrix>(*SP);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Z2);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  // // Banded = Banded, Id or Zero
  // ref = std::make_shared<SimpleMatrix>(*Band);
  // tRef = std::make_shared<SimpleMatrix>(*Band);
  // tRef->zero();
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);

  // ref = std::make_shared<SimpleMatrix>(*Z2);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  // ref = std::make_shared<SimpleMatrix>(*I2);
  // *tRef = *ref;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment2: ", (*tRef) == (*ref), true);
  std::cout << "-->  test assignment2 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators1() {
  std::cout << "--> Test: operators1." << std::endl;
  //+=, -=, *=, /=

  auto tmp = std::make_shared<SimpleMatrix>(*D);
  // Dense *=, /=
  double a = 2.2;
  int a1 = 2;
  *tmp *= a;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators: ", fabs((*tmp)(i, j) - a * (*D)(i, j)) < tol, true);

  *tmp *= a1;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators: ", fabs((*tmp)(i, j) - a * a1 * (*D)(i, j)) < tol, true);

  *tmp /= a;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators: ", fabs((*tmp)(i, j) - a1 * (*D)(i, j)) < tol, true);

  *tmp /= a1;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - (*D)(i, j)) < tol,
                                   true);

  // Dense +=, -= Dense

  *tmp += *SicM;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators: ", fabs((*tmp)(i, j) - (*SicM)(i, j) - (*D)(i, j)) < tol, true);

  *tmp -= *SicM;
  for (unsigned int i = 0; i < tmp->size(0); ++i)
    for (unsigned int j = 0; j < tmp->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*tmp)(i, j) - (*D)(i, j)) < tol,
                                   true);

  // Dense +=, -= Block
  C->zero();
  *C += *Ab;
  *C += *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*C)(i, j) - 2 * (*Ab)(i, j)) < tol,
                                   true);
  *C -= *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*C)(i, j) - (*Ab)(i, j)) < tol,
                                   true);

  std::cout << "-->  test operators1 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators2() {
  std::cout << "--> Test: operators2." << std::endl;
  // // +=, -=, *=, /= triangular
  // auto tmp = std::make_shared<SimpleMatrix>(*T);
  // auto tmp2 = std::make_shared<SimpleMatrix>(*T);
  // *tmp += *tmp2;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*T)(i, j), true);

  // int mult = 2;
  // double mult0 = 2.2;
  // *tmp *= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol, true);

  // *tmp *= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*T)(i, j)) < tol, true);

  // *tmp /= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol, true);

  // *tmp /= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*T)(i, j), true);

  // *tmp -= *tmp2;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getTriang() - *T) == 0, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", tmp->num() == siconos::algebra::UblasType::TRIANGULAR, true);

  std::cout << "-->  test operators2 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators3() {
  std::cout << "--> Test: operators3." << std::endl;
  // +=, -=, *=, /= Symmetric
  // auto tmp = std::make_shared<SimpleMatrix>(*S);
  // auto tmp2 = std::make_shared<SimpleMatrix>(*S);
  // *tmp += *tmp2;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*S)(i, j), true);

  // int mult = 2;
  // double mult0 = 2.2;
  // *tmp *= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol, true);

  // *tmp *= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*S)(i, j)) < tol, true);

  // *tmp /= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol, true);

  // *tmp /= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*S)(i, j), true);

  // *tmp -= *tmp2;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSym() - *S) == 0, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", tmp->num() == siconos::algebra::UblasType::SYMMETRIC, true);

  std::cout << "-->  test operators3 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators4() {
  std::cout << "--> Test: operators4." << std::endl;
  // +=, -=, *=, /= sparse
  // auto tmp = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp2 = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp3 = std::make_shared<SimpleMatrix>(*T2);

  // auto tmp4 = std::make_shared<SimpleMatrix>(*Band);
  // auto tmp5 = std::make_shared<SimpleMatrix>(*S2);

  // *tmp += *tmp2;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * (*SP)(i, j)) < tol,
  //                                  true);

  // int mult = 2;
  // double mult0 = 2.2;
  // *tmp *= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol, true);

  // *tmp *= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*SP)(i, j)) < tol, true);

  // *tmp /= mult;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol, true);

  // *tmp /= mult0;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2 * (*SP)(i, j)) < tol,
  //                                  true);

  // *tmp -= *tmp2;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSparse() - *SP) == 0, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", tmp->num() == siconos::algebra::UblasType::SPARSE, true);

  // // += -= a triangular
  // *tmp += *tmp3;
  // for (unsigned int i = 0; i < tmp->size(0); ++i) {
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol,
  //                                  true);
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp3)(i, j)) < tol, true);
  // }

  // *tmp -= *tmp3;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol,
  //                                  true);

  // // += -= a banded
  // *tmp -= *tmp;
  // *tmp += *tmp4;
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*Band)(i, j)) < tol,
  //                                  true);

  // *tmp -= *tmp4;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) < tol, true);

  // // += -= a sym

  // *tmp += *tmp5;
  // for (unsigned int i = 0; i < tmp->size(0); ++i) {
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp5)(j, i)) < tol, true);
  //   for (unsigned int j = i; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", ((*tmp)(i, j) - (*SP)(i, j) - (*tmp5)(i, j)) < tol, true);
  // }

  // *tmp -= *tmp5;
  // for (unsigned int i = 0; i < tmp->size(0); ++i)
  //   for (unsigned int j = 0; j < tmp->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP)(i, j)) < tol,
  //                                  true);

  std::cout << "-->  test operators4 ended with success." << std::endl;
}
void SimpleMatrixTest::testOperators4bis() {
  std::cout << "--> Test: operators4bis." << std::endl;
  // +=, -=, *=, /= sparse
//   auto tmp = std::make_shared<SimpleMatrix>(*SP_coor);
//   auto tmp2 = std::make_shared<SimpleMatrix>(*SP_coor);
//   auto tmp3 = std::make_shared<SimpleMatrix>(*T2);

//   auto tmp4 = std::make_shared<SimpleMatrix>(*Band);
//   auto tmp5 = std::make_shared<SimpleMatrix>(*S2);

//   auto tmp6 = std::make_shared<SimpleMatrix>(*SP);
//   *tmp += *tmp2;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * (*SP_coor)(i, j)) < tol, true);

//   int mult = 2;
//   double mult0 = 2.2;
//   *tmp *= mult0;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP_coor)(i, j)) < tol, true);

//   *tmp *= mult;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*SP_coor)(i, j)) < tol,
//           true);

//   *tmp /= mult;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP_coor)(i, j)) < tol, true);

//   *tmp /= mult0;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2 * (*SP_coor)(i, j)) < tol, true);

//   *tmp -= *tmp2;
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators: ", norm_inf(tmp->getSparseCoordinate() - *SP_coor) == 0, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators: ", tmp->num() == siconos::algebra::UblasType::SPARSE_COORDINATE, true);

//   // += -= a triangular
//   *tmp += *tmp3;
//   for (unsigned int i = 0; i < tmp->size(0); ++i) {
//     for (unsigned int j = 0; j < i; ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j)) < tol,
//                                    true);
//     for (unsigned int j = i; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j) - (*tmp3)(i, j)) < tol, true);
//   }

//   *tmp -= *tmp3;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j)) < tol,
//                                    true);

//   // += -= a banded
//   *tmp -= *tmp;
//   *tmp += *tmp4;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*Band)(i, j)) < tol,
//                                    true);

//   *tmp -= *tmp4;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) < tol, true);

//   // += -= a sym

//   *tmp += *tmp5;
//   for (unsigned int i = 0; i < tmp->size(0); ++i) {
//     for (unsigned int j = 0; j < i; ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j) - (*tmp5)(j, i)) < tol, true);
//     for (unsigned int j = i; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j) - (*tmp5)(i, j)) < tol, true);
//   }

//   *tmp -= *tmp5;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j)) < tol,
//                                    true);

//   // += -= a sparse

//   *tmp += *tmp6;
//   for (unsigned int i = 0; i < tmp->size(0); ++i) {
//     for (unsigned int j = 0; j < i; ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j) - (*tmp6)(j, i)) < tol, true);
//     for (unsigned int j = i; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j) - (*tmp6)(i, j)) < tol, true);
//   }

//   *tmp -= *tmp6;
//   for (unsigned int i = 0; i < tmp->size(0); ++i)
//     for (unsigned int j = 0; j < tmp->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - (*SP_coor)(i, j)) < tol,
//                                    true);

//   std::cout << "-->  test operators4bis ended with success." << std::endl;
// }
// void SimpleMatrixTest::testOperators5() {
//   std::cout << "--> Test: operators5." << std::endl;
//   // +=, -=, *=, /= banded
//   auto tmp = std::make_shared<SimpleMatrix>(*Band);
//   auto tmp2 = std::make_shared<SimpleMatrix>(*Band);
//   *tmp += *tmp2;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*Band)(i, j),
//                                    true);

//   int mult = 2;
//   double mult0 = 2.2;
//   *tmp *= mult0;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*Band)(i, j)) < tol, true);

//   *tmp *= mult;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*Band)(i, j)) < tol, true);

//   *tmp /= mult;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*Band)(i, j)) < tol, true);

//   *tmp /= mult0;
//   for (signed i = 0; i < signed(Band->size1()); ++i)
//     for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*Band)(i, j), true);

//   *tmp -= *tmp2;
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getBanded() - *Band) == 0,
//                                true);

//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators: ", tmp->num() == siconos::algebra::UblasType::BANDED, true);

  std::cout << "-->  test operators5 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators6() {
  std::cout << "--> Test: operator6." << std::endl;

  // // ============= C = A + B =============

  // // Dense = Dense + Dense
  // C->zero();
  // *C = *A + *B;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  // C->zero();
  // // Dense = Dense + Block
  // *C = *A + *Ab;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  // C->zero();
  // // Dense = Block + Dense
  // *C = *Ab + *A;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  // C->zero();

  // // Dense = Block + Block
  // *C = *Ab + *Ab;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*Ab)(i, j) - (*Ab)(i, j)) < tol, true);

  // C->zero();

  // // Block = Dense + Dense
  // Cb->zero();
  // *Cb = *A + *B;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  // // Block = Dense + Block

  // Cb->zero();
  // *Cb = *A + *Ab;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  // // Block = Block + Dense

  // Cb->zero();
  // *Cb = *Ab + *A;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);

  // // Block = Block + Block

  // Cb->zero();
  // Ab->display();
  // std::cout << "Bb:" << std::endl;
  // Bb->display();
  // *Cb = *Ab + *Bb;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*Bb)(i, j)) < tol, true);
  // Cb->zero();

  // // ============= C = A - B =============

  // // Dense = Dense - Dense
  // C->zero();
  // *C = *A - *B;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  // C->zero();
  // // Dense = Dense - Block
  // *C = *A - *Ab;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  // C->zero();
  // // Dense = Block - Dense
  // *C = *Ab - *A;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) + (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  // C->zero();

  // // Dense = Block - Block
  // *C = *Ab - *Bb;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = 0; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*C)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);

  // C->zero();

  // // Block = Dense - Dense
  // Cb->zero();
  // *Cb = *A - *B;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  // // Block = Dense - Block

  // Cb->zero();
  // *Cb = *A - *Ab;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  // // Block = Block - Dense

  // Cb->zero();
  // *Cb = *Ab - *A;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*A)(i, j)) < tol, true);

  // // Block = Block - Block

  // Cb->zero();
  // *Cb = *Ab - *Bb;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = 0; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);
  // Cb->zero();
  std::cout << "-->  test operators6 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators6Bis() {
  std::cout << "--> Test: operator6Bis." << std::endl;

  // ============= C = A + B =============

  // Dense = Dense + Dense
  C->zero();
  siconos::algebra::add(*A, *B, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense + Block
  siconos::algebra::add(*A, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block + Dense
  siconos::algebra::add(*Ab, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block + Block
  siconos::algebra::add(*Ab, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Block = Dense + Dense
  Cb->zero();
  siconos::algebra::add(*A, *B, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*B)(i, j)) < tol, true);

  // Block = Dense + Block

  Cb->zero();
  siconos::algebra::add(*A, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  // Block = Block + Dense

  Cb->zero();
  siconos::algebra::add(*Ab, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*A)(i, j)) < tol, true);

  // Block = Block + Block

  Cb->zero();
  siconos::algebra::add(*Ab, *Bb, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) - (*Bb)(i, j)) < tol, true);
  Cb->zero();

  // ============= C = A - B =============

  // Dense = Dense - Dense
  C->zero();
  siconos::algebra::sub(*A, *B, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  C->zero();
  // Dense = Dense - Block
  siconos::algebra::sub(*A, *Ab, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  C->zero();
  // Dense = Block - Dense
  siconos::algebra::sub(*Ab, *A, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) + (*A)(i, j) - (*Ab)(i, j)) < tol, true);

  C->zero();

  // Dense = Block - Block
  siconos::algebra::sub(*Ab, *Bb, *C);
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = 0; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);

  C->zero();

  // Block = Dense - Dense
  Cb->zero();
  siconos::algebra::sub(*A, *B, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*B)(i, j)) < tol, true);

  // Block = Dense - Block

  Cb->zero();
  siconos::algebra::sub(*A, *Ab, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) + (*Ab)(i, j)) < tol, true);

  // Block = Block - Dense

  Cb->zero();
  siconos::algebra::sub(*Ab, *A, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*A)(i, j)) < tol, true);

  // Block = Block - Block

  Cb->zero();
  siconos::algebra::sub(*Ab, *Bb, *Cb);
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = 0; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators6Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) + (*Bb)(i, j)) < tol, true);
  Cb->zero();
  std::cout << "-->  test operators6Bis ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators6Ter() {
  std::cout << "--> Test: operator6Ter." << std::endl;

  // +, - for non-dense matrices.

  // // Triang +,-,* Triang
  // auto tmp = std::make_shared<SimpleMatrix>(*T);
  // auto tmp2 = std::make_shared<SimpleMatrix>(*T);
  // auto res = std::make_shared<SimpleMatrix>(3, 3, siconos::algebra::UblasType::TRIANGULAR);
  // *res = *tmp + *tmp2;
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6Ter: ", (*res)(i, j) == ((*T)(i, j) + (*T)(i, j)), true);

  // *res = *tmp - *tmp2;
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  // // siconos::algebra::prod(*tmp, *tmp2, *res);
  // // siconos::algebra::prod(*T,*T, *tmp, true);
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(res->getTriang() - *tmp) ==
  // // 0, true);

  // // Sym +,-,* Sym
  // tmp = std::make_shared<SimpleMatrix>(*S);
  // tmp2 = std::make_shared<SimpleMatrix>(*S);
  // res = std::make_shared<SimpleMatrix>(3, 3, siconos::algebra::UblasType::SYMMETRIC);
  // *res = *tmp + *tmp2;
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6Ter: ", (*res)(i, j) == ((*S)(i, j) + (*S)(i, j)), true);

  // *res = *tmp - *tmp2;
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  // // siconos::algebra::prod(*tmp , *tmp2, *res, true);

  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(res->getSym() -
  // // ublas::prod(*S, *S))
  // // == 0, true);

  // // Sparse +,-,* Sparse
  // tmp = std::make_shared<SimpleMatrix>(*SP);
  // tmp2 = std::make_shared<SimpleMatrix>(*SP);
  // res = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // *res = *tmp + *tmp2;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == (2.0 * (*tmp)), true);

  // // *res = siconos::algebra::prod(*tmp , *tmp2);
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(*res->sparse() -
  // // ublas::prod(*SP, *SP)) < tol, true);

  // *res = *tmp - *tmp2;
  // tmp->zero();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == *tmp, true);

  // // SparseCoordinate +,-,* SparseCoordinate

  // tmp = std::make_shared<SimpleMatrix>(*SP_coor);
  // tmp2 = std::make_shared<SimpleMatrix>(*SP_coor);
  // res = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE_COORDINATE);
  // *res = *tmp + *tmp2;

  // // res->displayExpert();
  // // tmp->displayExpert();
  // // tmp2->displayExpert();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == (2.0 * (*tmp)), true);

  // // *res = siconos::algebra::prod(*tmp , *tmp2);
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", norm_inf(*res->sparseCoordinate() -
  // // ublas::prod(*SP_coor, *SP_coor)) < tol, true);

  // *res = *tmp - *tmp2;
  // tmp->zero();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res) == *tmp, true);

  // // Banded +,- Banded
  // tmp = std::make_shared<SimpleMatrix>(*Band);
  // tmp2 = std::make_shared<SimpleMatrix>(*Band);
  // res = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::BANDED);
  // *res = *tmp + *tmp2;
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators6Ter: ", (*res)(i, j) == ((*Band)(i, j) + (*Band)(i, j)), true);
  // *res = *tmp - *tmp2;
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators6Ter: ", (*res)(i, j) == 0, true);

  std::cout << "-->  test operators6Ter6 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators7() {
  std::cout << "--> Test: operator7." << std::endl;
  auto tmp1 = std::make_shared<SimpleMatrix>(*D);
  tmp1->resize(4, 4);
  // auto tmp2 = std::make_shared<SimpleMatrix>(*T2);
  // auto tmp3 = std::make_shared<SimpleMatrix>(*S2);
  // auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp5 = std::make_shared<SimpleMatrix>(*Band);
  // auto tmp6 = std::make_shared<SimpleMatrix>(*Z2);
  // auto tmp7 = std::make_shared<SimpleMatrix>(*I2);

  auto res = std::make_shared<SimpleMatrix>(4, 4);

  // dense + ...
  // ... triang
  // add(*tmp1, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // add(*tmp1, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }
  // // ... Sparse
  // add(*tmp1, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = 0; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp4)(i, j)) < tol, true);
  // // ... Banded
  // add(*tmp1, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*Band)(i, j)) < tol, true);
  // // Zero
  // add(*tmp1, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp1, true);

  // // Id
  // add(*tmp1, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = 0; j < res->size(1); ++j) {
  //     if (i == j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - 1) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  //   }
  // }

  // dense - ...
  // ... triangular
  // sub(*tmp1, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // sub(*tmp1, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp3)(j, i)) < tol, true);
  // }
  // // ... Sparse
  // sub(*tmp1, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = 0; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*tmp4)(i, j)) < tol, true);
  // // ... Banded
  // sub(*tmp1, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + (*Band)(i, j)) < tol, true);

  // // Zero
  // sub(*tmp1, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp1, true);

  // // Id
  // sub(*tmp1, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = 0; j < res->size(1); ++j) {
  //     if (i == j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) + 1) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol, true);
  //   }
  // }
  // // triang + ...
  // // ... dense
  // add(*tmp2, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // add(*tmp2, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol,
  //                                  true);
  // }
  // // ... Sparse
  // add(*tmp2, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol,
  //                                  true);
  // }

  // // ... Banded
  // add(*tmp2, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*Band)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*Band)(i, j)) < tol, true);

  // // ... Zero
  // add(*tmp2, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp2, true);

  // // ... Identity
  // add(*tmp2, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp7)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j)) < tol, true);
  // }

  // // triang - ...
  // // ... dense
  // sub(*tmp2, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp1)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // sub(*tmp2, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol,
  //                                  true);
  // }
  // // ... Sparse
  // sub(*tmp2, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp4)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol,
  //                                  true);
  // }

  // // ... Banded
  // sub(*tmp2, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*Band)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*Band)(i, j)) < tol, true);

  // // ... Zero
  // sub(*tmp2, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp2, true);

  // // Identity
  // sub(*tmp2, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) + (*tmp7)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j)) < tol, true);
  // }

  // // sym + ...
  // // ... dense
  // add(*tmp3, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }
  // // ... triang
  // add(*tmp3, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol,
  //                                  true);
  // }
  // // ... Sparse
  // add(*tmp3, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }

  // // ... Banded
  // add(*tmp3, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*Band)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) - (*Band)(i, j)) < tol, true);

  // // ... Zero
  // add(*tmp3, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp3, true);

  // // ... identity
  // add(*tmp3, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }

  // // sym - ...
  // // ... dense
  // sub(*tmp3, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp1)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp1)(i, j)) < tol, true);
  // }
  // // ... triang
  // sub(*tmp3, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol,
  //                                  true);
  // }
  // // ... Sparse
  // sub(*tmp3, *tmp4, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp4)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp4)(i, j)) < tol, true);
  // }

  // // ... Banded
  // sub(*tmp3, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*Band)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*Band)(i, j)) < tol, true);

  // // ... Zero
  // sub(*tmp3, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp3, true);
  // // Identity
  // sub(*tmp3, *tmp7, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) + (*tmp7)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) + (*tmp7)(i, j)) < tol, true);
  // }

  // // sparse + ...
  // // ... dense
  // add(*tmp4, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = 0; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp4)(i, j)) < tol, true);
  // // ... triang
  // add(*tmp4, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // add(*tmp4, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }
  // // ... Banded
  // add(*tmp4, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*Band)(i, j)) < tol, true);

  // // ... zero
  // add(*tmp4, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp4, true);

  // // sparse - ...
  // // ... dense
  // sub(*tmp4, *tmp1, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i)
  //   for (unsigned int j = 0; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp1)(i, j)) < tol, true);
  // // ... triangular
  // sub(*tmp4, *tmp2, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp2)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol,
  //                                  true);
  // }
  // // ... Sym
  // sub(*tmp4, *tmp3, *res);
  // for (unsigned int i = 0; i < res->size(0); ++i) {
  //   for (unsigned int j = i; j < res->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp3)(i, j)) < tol, true);
  //   for (unsigned int j = 0; j < i; ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*tmp3)(j, i)) < tol, true);
  // }

  // // ... Banded
  // sub(*tmp4, *tmp5, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) + (*Band)(i, j)) < tol, true);

  // // ... zero
  // sub(*tmp4, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp4, true);

  // // Banded + ...
  // // ... dense
  // add(*tmp5, *tmp1, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp1)(i, j)) < tol,
  //                                  true);
  // }
  // // ... triang
  // add(*tmp5, *tmp2, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j)) < tol, true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j) - (*tmp5)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp2)(i, j)) < tol, true);
  // }

  // // ...sym
  // add(*tmp5, *tmp3, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j) - (*tmp5)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp3)(j, i)) < tol, true);
  // }

  // //... sparse
  // add(*tmp5, *tmp4, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size1())); j < signed(Band->size1()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp4)(i, j)) < tol,
  //                                  true);
  // }

  // // ... zero
  // add(*tmp5, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp5, true);
  // // ... identity
  // add(*tmp5, *tmp7, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) - (*tmp7)(i, j)) < tol,
  //                                  true);
  // }

  // // Banded - ...
  // // ... dense

  // sub(*tmp5, *tmp1, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp1)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp1)(i, j)) < tol,
  //                                  true);
  // }

  // // ... triang
  // sub(*tmp5, *tmp2, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp2)(i, j)) < tol, true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp2)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp2)(i, j)) < tol, true);
  // }

  // // ...sym
  // sub(*tmp5, *tmp3, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp3)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol, true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp3)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp3)(j, i)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     if (j >= i)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp3)(i, j)) < tol, true);
  //     else
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //           "testOperators: ", fabs((*res)(i, j) + (*tmp3)(j, i)) < tol, true);
  // }

  // //... sparse
  // sub(*tmp5, *tmp4, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j) - (*tmp5)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size1())); j < signed(Band->size1()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp4)(i, j)) < tol,
  //                                  true);
  // }

  // // ... zero
  // sub(*tmp5, *tmp6, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp5, true);
  // // ... identity
  // sub(*tmp5, *tmp7, *res);
  // for (signed i = 0; i < signed(Band->size1()); ++i) {
  //   for (signed j = 0; j < std::max(i - 1, 0); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp7)(i, j)) < tol,
  //                                  true);
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(Band->size2())); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*res)(i, j) - (*tmp5)(i, j) + (*tmp7)(i, j)) < tol, true);
  //   for (signed j = std::min(i + 2, signed(Band->size2())); j < signed(Band->size2()); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*res)(i, j) + (*tmp7)(i, j)) < tol,
  //                                  true);
  // }

  std::cout << "-->  test operators7 ended with success." << std::endl;
}

// void SimpleMatrixTest::testOperators8()
//  {
//    std::cout << "--> Test: operator8." <<std::endl;

//   // // Simple = Simple * Simple
//   // *C = ublas::prod(*A, *B);
//   // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8: ", norm_inf(*C->dense() -
//   ublas::prod(*A->dense(), *B->dense())) < tol, true);

//   // Block = Simple * Simple
//   // *Cb = ublas::prod(*A, *B);
//   // siconos::algebra::DenseMat Dtmp = ublas::prod(*A->dense(), *B->dense());
//   // auto tmp = std::make_shared<SimpleMatrix>(Dtmp);
//   // for (unsigned int i = 0; i < C->size(0); ++i)
//   //   for (unsigned int j = i ; j < C->size(1); ++j)
//   //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) <
//   tol, true);

//   // Others ...

//   auto tmp1 = std::make_shared<SimpleMatrix>(4, 4, 2.3);
//   auto tmp2 = std::make_shared<SimpleMatrix>(*T2);
//   auto tmp3 = std::make_shared<SimpleMatrix>(*S2);
//   auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
//   auto tmp5 = std::make_shared<SimpleMatrix>(*Band);

//   auto res = std::make_shared<SimpleMatrix>(4, 4, 0);

//   // Dense * ...
//   // triang
//   *res = ublas::prod(*tmp1, *tmp2);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp1->getDense(), tmp2->getTriang())) < tol, true);
//   // Sym
//   *res = ublas::prod(*tmp1, *tmp3);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp1->getDense(), tmp3->getSym())) < tol, true);
//   // Sparse
//   *res = ublas::prod(*tmp1, *tmp4);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp1->getDense(), tmp4->getSparse())) < tol, true);
//   // Banded
//   *res = ublas::prod(*tmp1, *tmp5);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp1->getDense(), tmp5->getBanded())) < tol, true);
//   // triang * ...
//   // dense
//   *res = ublas::prod(*tmp2, *tmp1);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp2->getTriang(), tmp1->getDense())) < tol, true);
//   // Sym
//   *res = ublas::prod(*tmp2, *tmp3);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp2->getTriang(), tmp3->getSym())) < tol, true);
//   // Sparse
//   *res = ublas::prod(*tmp2, *tmp4);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp2->getTriang(), tmp4->getSparse())) < tol, true);
//   // Banded
//   *res = ublas::prod(*tmp2, *tmp5);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp2->getTriang(), tmp5->getBanded())) < tol, true);
//   // sym * ...
//   // dense
//   *res = ublas::prod(*tmp3, *tmp1);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp3->getSym(), tmp1->getDense())) < tol, true);
//   // triang
//   *res = ublas::prod(*tmp3, *tmp2);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp3->getSym(), tmp2->getTriang())) < tol, true);
//   // Sparse
//   *res = ublas::prod(*tmp3, *tmp4);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp3->getSym(), tmp4->getSparse())) < tol, true);
//   // Banded
//   *res = ublas::prod(*tmp3, *tmp5);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp3->getSym(), tmp5->getBanded())) < tol, true);
//   // Sparse * ...
//   // dense
//   *res = ublas::prod(*tmp4, *tmp1);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp4->getSparse(), tmp1->getDense())) < tol, true);
//   // triang
//   *res = ublas::prod(*tmp4, *tmp2);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp4->getSparse(), tmp2->getTriang())) < tol, true);
//   // Sym
//   *res = ublas::prod(*tmp4, *tmp3);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp4->getSparse(), tmp3->getSym())) < tol, true);
//   // Banded
//   *res = ublas::prod(*tmp4, *tmp5);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp4->getSparse(), tmp5->getBanded())) < tol, true);
//   // Banded * ...
//   // dense
//   *res = ublas::prod(*tmp5, *tmp1);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp5->getBanded(), tmp1->getDense())) < tol, true);
//   // triang
//   *res = ublas::prod(*tmp5, *tmp2);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp5->getBanded(), tmp2->getTriang())) < tol, true);
//   // Sparse
//   *res = ublas::prod(*tmp5, *tmp4);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp5->getBanded(), tmp4->getSparse())) < tol, true);
//   // Sym
//   *res = ublas::prod(*tmp5, *tmp3);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(*res->dense() -
//   ublas::prod(tmp5->getBanded(), tmp3->getSym())) < tol, true);

//   std::cout << "-->  test operators8 ended with success." <<std::endl;
// }

void SimpleMatrixTest::testOperators8Bis() {
  std::cout << "--> Test: operator8Bis." << std::endl;
  // Simple = Simple * Simple
  // siconos::algebra::prod(*A, *B, *C);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators8Bis: ",
  //     norm_inf(*C->dense() - ublas::prod(*A->dense(), *B->dense())) < tol, true);

  // // Block = Simple * Simple
  // siconos::algebra::prod(*A, *B, *Cb);
  // siconos::algebra::DenseMat Dtmp = ublas::prod(*A->dense(), *B->dense());
  // auto tmp = std::make_shared<SimpleMatrix>(Dtmp);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) < tol,
  //                                  true);

  // Others ...

  // Others ...
//   auto tmp1 = std::make_shared<SimpleMatrix>(4, 4, 2.4);
//   auto tmp2 = std::make_shared<SimpleMatrix>(*T2);
//   auto tmp3 = std::make_shared<SimpleMatrix>(*S2);
//   auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
//   auto tmp5 = std::make_shared<SimpleMatrix>(*Band);

//   auto res = std::make_shared<SimpleMatrix>(4, 4);

//   // Dense * ...
//   // triang
//   siconos::algebra::prod(*tmp1, *tmp2, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp1->getDense(), tmp2->getTriang())) < tol, true);
//   // Sym
//   siconos::algebra::prod(*tmp1, *tmp3, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp1->getDense(), tmp3->getSym())) < tol, true);
//   // Sparse
//   siconos::algebra::prod(*tmp1, *tmp4, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp1->getDense(), tmp4->getSparse())) < tol, true);
//   // Banded
//   siconos::algebra::prod(*tmp1, *tmp5, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp1->getDense(), tmp5->getBanded())) < tol, true);
//   // triang * ...
//   // dense
//   siconos::algebra::prod(*tmp2, *tmp1, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp2->getTriang(), tmp1->getDense())) < tol, true);
//   // Sym
//   siconos::algebra::prod(*tmp2, *tmp3, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp2->getTriang(), tmp3->getSym())) < tol, true);
//   // Sparse
//   siconos::algebra::prod(*tmp2, *tmp4, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp2->getTriang(), tmp4->getSparse())) < tol, true);
//   // Banded
//   siconos::algebra::prod(*tmp2, *tmp5, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp2->getTriang(), tmp5->getBanded())) < tol, true);
//   // sym * ...
//   // dense
//   siconos::algebra::prod(*tmp3, *tmp1, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp3->getSym(), tmp1->getDense())) < tol, true);
//   // triang
//   siconos::algebra::prod(*tmp3, *tmp2, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp3->getSym(), tmp2->getTriang())) < tol, true);
//   // Sparse
//   siconos::algebra::prod(*tmp3, *tmp4, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp3->getSym(), tmp4->getSparse())) < tol, true);
//   // Banded
//   siconos::algebra::prod(*tmp3, *tmp5, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp3->getSym(), tmp5->getBanded())) < tol, true);
//   // Sparse * ...
//   // dense
//   siconos::algebra::prod(*tmp4, *tmp1, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp4->getSparse(), tmp1->getDense())) < tol, true);
//   // triang
//   siconos::algebra::prod(*tmp4, *tmp2, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp4->getSparse(), tmp2->getTriang())) < tol, true);
//   // Sym
//   siconos::algebra::prod(*tmp4, *tmp3, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp4->getSparse(), tmp3->getSym())) < tol, true);
//   // Banded
//   siconos::algebra::prod(*tmp4, *tmp5, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp4->getSparse(), tmp5->getBanded())) < tol, true);
//   // Banded * ...
//   // dense
//   siconos::algebra::prod(*tmp5, *tmp1, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp5->getBanded(), tmp1->getDense())) < tol, true);
//   // triang
//   siconos::algebra::prod(*tmp5, *tmp2, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp5->getBanded(), tmp2->getTriang())) < tol, true);
//   // Sparse
//   siconos::algebra::prod(*tmp5, *tmp4, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp5->getBanded(), tmp4->getSparse())) < tol, true);
//   // Sym
//   siconos::algebra::prod(*tmp5, *tmp3, *res);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Bis: ",
//       norm_inf(*res->dense() - ublas::prod(tmp5->getBanded(), tmp3->getSym())) < tol, true);

//   std::cout << "-->  test operators8Bis ended with success." << std::endl;
// }

// void SimpleMatrixTest::testOperators8Ter() {
//   std::cout << "--> Test: operator8Ter." << std::endl;
//   // Simple = Simple * Simple
//   siconos::algebra::axpy_prod(*A, *B, *C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Ter: ",
//       norm_inf(*C->dense() - ublas::prod(*A->dense(), *B->dense())) < tol, true);

//   // Simple += Simple * Simple
//   auto backUp = std::make_shared<SimpleMatrix>(*C);

//   siconos::algebra::axpy_prod(*A, *B, *C, false);

//   CPPUNIT_ASSERT_EQUAL_MESSAGE(
//       "testOperators8Ter: ",
//       norm_inf(*C->dense() - ublas::prod(*A->dense(), *B->dense()) - *backUp->dense()) < tol,
//       true);
//   // Block = Simple * Simple
//   siconos::algebra::axpy_prod(*A, *B, *Cb, true);
//   siconos::algebra::DenseMat Dtmp = ublas::prod(*A->dense(), *B->dense());
//   auto tmp = std::make_shared<SimpleMatrix>(Dtmp);
//   for (unsigned int i = 0; i < Cb->size(0); ++i)
//     for (unsigned int j = i; j < Cb->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j)) < tol,
//                                    true);

//   *backUp = *Cb;
//   // Block += Simple * Simple
//   siconos::algebra::axpy_prod(*A, *B, *Cb, false);
//   Dtmp = ublas::prod(*A->dense(), *B->dense());
//   *tmp = Dtmp;
//   for (unsigned int i = 0; i < Cb->size(0); ++i)
//     for (unsigned int j = i; j < Cb->size(1); ++j)
//       CPPUNIT_ASSERT_EQUAL_MESSAGE(
//           "testOperators: ", fabs((*Cb)(i, j) - (*tmp)(i, j) - (*backUp)(i, j)) < tol, true);

  std::cout << "-->  test operators8Ter ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators8_4()  // C += A*B
{
  std::cout << "--> Test: operator8_4." << std::endl;
  // // Simple = Simple * Simple
  // C->zero();
  // siconos::algebra::prod(*A, *B, *C, false);
  // siconos::algebra::prod(*A, *B, *C, false);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators8_4: ",
  //     norm_inf(*C->dense() - 2 * ublas::prod(*A->dense(), *B->dense())) < tol, true);

  // // Block = Simple * Simple
  // Cb->zero();
  // siconos::algebra::prod(*A, *B, *Cb, false);
  // siconos::algebra::prod(*A, *B, *Cb, false);
  // siconos::algebra::DenseMat Dtmp = ublas::prod(*A->dense(), *B->dense());
  // auto tmp = std::make_shared<SimpleMatrix>(Dtmp);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators: ", fabs((*Cb)(i, j) - 2 * (*tmp)(i, j)) < tol, true);
  std::cout << "-->  test operators8_4 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators8_5() {
  // == Test siconos::algebra::subprod ==

  std::cout << "--> Test: operator8_5." << std::endl;
  std::vector<std::size_t> coord(8);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(2);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto x3 = std::make_shared<siconos::algebra::SiconosVector>(5);
  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::BlockVector>();
  auto v = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->insertPtr(x1);
  x->insertPtr(x2);
  x->insertPtr(x3);
  for (unsigned int i = 0; i < size; ++i) {
    (*x)(i) = (double)i + 3;
    (*v)(i) = (double)i + 3;
  }

  // v == x but x is a 3-blocks vector.

  // Simple = Simple * Simple, all dense
  // siconos::algebra::subprod but with full matrix/vectors
  coord[0] = 0;
  coord[1] = size;
  coord[2] = 0;
  coord[3] = size;
  coord[4] = 0;
  coord[5] = size;
  coord[6] = 0;
  coord[7] = size;
  siconos::algebra::subprod(*A, *v, *y, coord, true);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators8_5: ",
  //     (*y - (siconos::algebra::prod(*A, *v))).normInf() < tol, true); // TODO

  // Simple = Simple * Block, all dense
  // siconos::algebra::subprod but with full matrix/vectors
  //  siconos::algebra::subprod(*A,*x,*y, coord, true);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", norm_inf(*y->dense()-
  //  ublas::prod(*A->dense(),*v->dense()))<tol, true);

  coord[0] = 0;
  coord[1] = 2;
  coord[2] = 1;
  coord[3] = 3;
  coord[4] = 3;
  coord[5] = 5;
  coord[6] = 2;
  coord[7] = 4;
  y->zero();
  // Simple = Simple * Simple, all dense
  siconos::algebra::subprod(*A, *v, *y, coord, true);
  double res = (*A)(0, 1) * (*v)(3) + (*A)(0, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*v)(3) + (*A)(1, 2) * (*v)(4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i) {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  }
  y->zero();
  // Simple = Simple * Block, all dense
  //  siconos::algebra::subprod(*A,*x,*y, coord, true);
  //  res = (*A)(0,1)*(*x)(3) + (*A)(0,2)*(*x)(4);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res-(*y)(2))<tol, true);
  //  res = (*A)(1,1)*(*x)(3) + (*A)(1,2)*(*x)(4);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res-(*y)(3))<tol, true);
  //  for (unsigned int i=0; i<size; ++i)
  //  {
  //    if (i!=2 && i!=3)
  //      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i))<tol, true);
  //  }
  //   // Others ...
  // Triang

  // auto A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::TRIANGULAR);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) (*A2)(i, j) = 3 * i + j;

  // siconos::algebra::subprod(*A2, *v, *y, coord, true);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  // }
  // // Sym
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::SYMMETRIC);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) (*A2)(i, j) = 3 * i + j;

  // siconos::algebra::subprod(*A2, *v, *y, coord, true);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  // }

  // // Sparse
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::SPARSE);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) A2->setValue(i, j, 3 * i + j);

  // siconos::algebra::subprod(*A2, *v, *y, coord, true);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  // }

  // // Banded
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::BANDED);
  // for (signed i = 0; i < signed(A2->size(0)); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(A2->size(1))); ++j)
  //     (*A2)(i, j) = 3 * i + j;
  // siconos::algebra::subprod(*A2, *v, *y, coord, true);
  // res = (*A2)(0, 1) * (*v)(3);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_5: ", fabs((*y)(i)) < tol, true);
  // }

  std::cout << "-->  test operators8_5 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators8_6() {
  // == Test siconos::algebra::subprod, with += ==

  std::cout << "--> Test: operator8_6." << std::endl;
  std::vector<std::size_t> coord(8);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(2);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(3);
  auto x3 = std::make_shared<siconos::algebra::SiconosVector>(5);
  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::BlockVector>();
  auto v = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->insertPtr(x1);
  x->insertPtr(x2);
  x->insertPtr(x3);
  for (unsigned int i = 0; i < size; ++i) {
    (*x)(i) = (double)i + 3;
    (*v)(i) = (double)i + 3;
  }

  // v == x but x is a 3-blocks vector.

  *y = *v;

  // Simple = Simple * Simple, all dense
  // siconos::algebra::subprod but with full matrix/vectors
  coord[0] = 0;
  coord[1] = size;
  coord[2] = 0;
  coord[3] = size;
  coord[4] = 0;
  coord[5] = size;
  coord[6] = 0;
  coord[7] = size;
  siconos::algebra::subprod(*A, *v, *y, coord, false);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators8_6: ",
  //     (*y - siconos::algebra::prod(*A, *v) - *v).normInf() < tol, true); // TODO

  // Simple = Simple * Block, all dense
  // siconos::algebra::subprod but with full matrix/vectors
  *y = *v;
  //  siconos::algebra::subprod(*A,*x,*y, coord, false);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", norm_inf(*y->dense()-
  //  ublas::prod(*A->dense(),*v->dense())- *v->dense())<tol, true);

  coord[0] = 0;
  coord[1] = 2;
  coord[2] = 1;
  coord[3] = 3;
  coord[4] = 3;
  coord[5] = 5;
  coord[6] = 2;
  coord[7] = 4;

  // Simple = Simple * Simple, all dense
  *y = *v;
  siconos::algebra::subprod(*A, *v, *y, coord, false);
  double res = (*A)(0, 1) * (*v)(3) + (*A)(0, 2) * (*v)(4) + (*v)(2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  res = (*A)(1, 1) * (*v)(3) + (*A)(1, 2) * (*v)(4) + (*v)(3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  for (unsigned int i = 0; i < size; ++i) {
    if (i != 2 && i != 3)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  }
  *y = *v;
  // Simple = Simple * Block, all dense
  //  siconos::algebra::subprod(*A,*x,*y, coord, false);
  //  res = (*A)(0,1)*(*x)(3) + (*A)(0,2)*(*x)(4) + (*v)(2);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res-(*y)(2))<tol, true);
  //  res = (*A)(1,1)*(*x)(3) + (*A)(1,2)*(*x)(4) + (*v)(3);
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res-(*y)(3))<tol, true);
  //  for (unsigned int i=0; i<size; ++i)
  //  {
  //    if (i!=2 && i!=3)
  //      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i)-(*v)(i))<tol, true);
  //  }

  //   // Others ...
  // Triang

  // auto A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::TRIANGULAR);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) (*A2)(i, j) = 3 * i + j;

  // *y = *v;
  // siconos::algebra::subprod(*A2, *v, *y, coord, false);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  // }

  // // Sym
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::SYMMETRIC);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) (*A2)(i, j) = 3 * i + j;

  // *y = *v;
  // siconos::algebra::subprod(*A2, *v, *y, coord, false);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  // }

  // // Sparse
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::SPARSE);
  // for (unsigned i = 0; i < A2->size(0); ++i)
  //   for (unsigned j = i; j < A2->size(1); ++j) A2->setValue(i, j, 3 * i + j);

  // *y = *v;
  // siconos::algebra::subprod(*A2, *v, *y, coord, false);
  // res = (*A2)(0, 1) * (*v)(3) + (*A2)(0, 2) * (*v)(4) + (*v)(2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  // }

  // // Banded
  // A2 = std::make_shared<SimpleMatrix>(10, 10, siconos::algebra::UblasType::BANDED);
  // for (signed i = 0; i < signed(A2->size(0)); ++i)
  //   for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(A2->size(1))); ++j)
  //     (*A2)(i, j) = 3 * i + j;
  // *y = *v;
  // siconos::algebra::subprod(*A2, *v, *y, coord, false);
  // res = (*A2)(0, 1) * (*v)(3) + (*v)(2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(2)) < tol, true);
  // res = (*A2)(1, 1) * (*v)(3) + (*A2)(1, 2) * (*v)(4) + (*v)(3);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs(res - (*y)(3)) < tol, true);
  // for (unsigned int i = 0; i < size; ++i) {
  //   if (i != 2 && i != 3)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators8_6: ", fabs((*y)(i) - (*v)(i)) < tol, true);
  // }

  std::cout << "-->  test operators8_6 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators9() {
  std::cout << "--> Test: operator9." << std::endl;

  // C = a*A or A/a

  double a = 2.2;
  int a1 = 3;

  // Simple = a * Simple or Simple/a
  *C = a * *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol,
                                   true);
  *C = a1 * *A;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a1 * (*A)(i, j)) < tol, true);

  // *C = *A / a;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i ; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a) <
  //     tol, true);
  // *C = *A / a1;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i ; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*A)(i, j) / a1) <
  //     tol, true);

  // Simple = a * Block

  *C = a * *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);
  ;
  *C = a1 * *Ab;
  for (unsigned int i = 0; i < C->size(0); ++i)
    for (unsigned int j = i; j < C->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*C)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  // *C = *Ab / a;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i ; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a) <
  //     tol, true);
  // *C = *Ab / a1;
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i ; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*C)(i, j) - (*Ab)(i, j) / a1) <
  //     tol, true);

  // Block = a * Block
  *Cb = a * *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);
  *Cb = a1 * *Ab;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a1 * (*Ab)(i, j)) < tol, true);

  // *Cb = *Ab / a;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i ; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a) <
  //     tol, true);
  // *Cb = *Ab / a1;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i ; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a1)
  //     < tol, true);

  // Block = a * Simple
  *Cb = a * *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);
  *Cb = a1 * *A;
  for (unsigned int i = 0; i < Cb->size(0); ++i)
    for (unsigned int j = i; j < Cb->size(1); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE(
          "testOperators9: ", fabs((*Cb)(i, j) - a1 * (*A)(i, j)) < tol, true);

  // *Cb = *A / a;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i ; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a) <
  //     tol, true);
  // *Cb = *A / a1;
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i ; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators9: ", fabs((*Cb)(i, j) - (*A)(i, j) / a1) <
  //     tol, true);
  std::cout << "-->  test operators9 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators9Bis() {
  std::cout << "--> Test: operator9Bis." << std::endl;

  // C = a*A or A/a

  double a = 2.2;

  // Simple = a * Simple or Simple/a
  // siconos::algebra::scal(a, *A, *C);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*C)(i, j) - a * (*A)(i, j)) < tol, true);

  // siconos::algebra::scal(1.0 / a, *A, *C);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*C)(i, j) - (*A)(i, j) / a) < tol, true);
  // // Simple = a * Block

  // siconos::algebra::scal(a, *Ab, *C);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*C)(i, j) - a * (*Ab)(i, j)) < tol, true);

  // siconos::algebra::scal(1.0 / a, *Ab, *C);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*C)(i, j) - (*Ab)(i, j) / a) < tol, true);

  // // Block = a * Block
  // siconos::algebra::scal(a, *Ab, *Cb);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*Cb)(i, j) - a * (*Ab)(i, j)) < tol, true);

  // siconos::algebra::scal(1.0 / a, *Ab, *Cb);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*Cb)(i, j) - (*Ab)(i, j) / a) < tol, true);

  // // Block = a * Simple
  // siconos::algebra::scal(a, *A, *Cb);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*Cb)(i, j) - a * (*A)(i, j)) < tol, true);

  // siconos::algebra::scal(1.0 / a, *A, *Cb);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Bis: ", fabs((*Cb)(i, j) - (*A)(i, j) / a) < tol, true);
  std::cout << "-->  test operators9Bis ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators9Ter() {
  std::cout << "--> Test: operator9Ter." << std::endl;

  // C += a*A or A/a

  double a = 2.2;
  C->zero();
  // // Simple = a * Simple or Simple/a
  // siconos::algebra::scal(a, *A, *C, false);
  // siconos::algebra::scal(a, *A, *C, false);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Ter: ", fabs((*C)(i, j) - 2 * a * (*A)(i, j)) < tol, true);

  // // Simple = a * Block
  // C->zero();
  // siconos::algebra::scal(a, *Ab, *C, false);
  // siconos::algebra::scal(a, *Ab, *C, false);
  // for (unsigned int i = 0; i < C->size(0); ++i)
  //   for (unsigned int j = i; j < C->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Ter: ", fabs((*C)(i, j) - 2 * a * (*Ab)(i, j)) < tol, true);

  // // Block = a * Block
  // Cb->zero();
  // siconos::algebra::scal(a, *Ab, *Cb, false);
  // siconos::algebra::scal(a, *Ab, *Cb, false);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Ter: ", fabs((*Cb)(i, j) - 2 * a * (*Ab)(i, j)) < tol, true);

  // // Block = a * Simple
  // Cb->zero();
  // siconos::algebra::scal(a, *A, *Cb, false);
  // siconos::algebra::scal(a, *A, *Cb, false);
  // for (unsigned int i = 0; i < Cb->size(0); ++i)
  //   for (unsigned int j = i; j < Cb->size(1); ++j)
  //     CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //         "testOperators9Ter: ", fabs((*Cb)(i, j) - 2 * a * (*A)(i, j)) < tol, true);

  std::cout << "-->  test operators9Ter ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators10() {
  std::cout << "--> Test: operator10." << std::endl;
  // double m = 2.2;
  // int i = 3;
  // auto tmp1 = std::make_shared<SimpleMatrix>(*T);
  // auto res = std::make_shared<SimpleMatrix>(3, 3, siconos::algebra::UblasType::TRIANGULAR);
  // *res = m * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() * m) < tol, true);
  // *res = i * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() * i) < tol, true);
  // // *res = *tmp1 / m;
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() -
  // // tmp1->getTriang() / m) < tol, true); *res = *tmp1 / i;
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() -
  // // tmp1->getTriang() / i) < tol, true);
  std::cout << "-->  test operators10 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators11() {
  std::cout << "--> Test: operator11." << std::endl;
  // double m = 2.2;
  // int i = 3;

  // auto tmp1 = std::make_shared<SimpleMatrix>(*S);
  // auto res = std::make_shared<SimpleMatrix>(3, 3, siconos::algebra::UblasType::SYMMETRIC);
  // *res = m * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getSym() - tmp1->getSym() * m) < tol, true);
  // *res = i * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getSym() - tmp1->getSym() * i) < tol, true);
  // // *res = *tmp1 / m;
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym() /
  // // m) < tol, true); *res = *tmp1 / i; CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ",
  // // norm_inf(res->getSym() - tmp1->getSym() / i) < tol, true);
  std::cout << "-->  test operator11 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators12() {
  std::cout << "--> Test: operator12." << std::endl;
  double m = 2.2;
  int i = 3;
  // auto tmp1 = std::make_shared<SimpleMatrix>(*SP);
  // auto res = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // *res = m * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() * m) < tol, true);
  // *res = i * *tmp1;
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() * i) < tol, true);
  // // *res = *tmp1 / m;
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() -
  // // tmp1->getSparse() / m) < tol, true); *res = *tmp1 / i;
  // // CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() -
  // // tmp1->getSparse() / i) < tol, true);
  std::cout << "-->  test operators12 ended with success." << std::endl;
}

void SimpleMatrixTest::testOperators13() {
  std::cout << "--> Test: operator13." << std::endl;
  //   double m = 2.2;
  //   int i = 3;
  //   auto tmp1 = std::make_shared<SimpleMatrix>(*Band);
  //   auto res = std::make_shared<SimpleMatrix>(*Band);//4,4,BANDED,1,1);
  //   *res = m * *tmp1;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()*m)<tol, true); *res = i ** tmp1;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()*i)<tol, true); *res = *tmp1 * m;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()*m)<tol, true); *res = *tmp1 * i;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()*i)<tol, true); *res = *tmp1 / m;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()/m)<tol, true); *res = *tmp1 / i;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded()-
  //   tmp1->getBanded()/i)<tol, true);
  std::cout << "-->  test operators13 ended with success." << std::endl;
}

void SimpleMatrixTest::testProd()  // y = A*x
{
  std::cout << "--> Test: ublas::prod. mat-vect" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->zero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  *y = siconos::algebra::prod(*A, *x);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  //  *y = siconos::algebra::prod(*A , *xB);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*y)(i) - sum)< tol, true);
  //  }

  // Block = Simple * Simple
  *yB = siconos::algebra::prod(*A, *x);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;

    for (unsigned int j = 0; j < A->size(1); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*yB)(i)-sum) < tol, true);
  }

  // Block = Simple * Block
  //  *yB = siconos::algebra::prod(*A ,*xB);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd: ", fabs((*yB)(i) - sum)< tol, true);
  //  }

  // Others or old stuff ...

  // auto tmp2 = std::make_shared<SimpleMatrix>(*T);
  // auto tmp3 = std::make_shared<SimpleMatrix>(*S);
  // auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp5 = std::make_shared<SimpleMatrix>(*Band2);
  // auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  // (*v)(0) = 1;
  // (*v)(1) = 2;
  // (*v)(2) = 3;
  // auto vv = std::make_shared<siconos::algebra::SiconosVector>(4);
  // (*vv)(0) = 1;
  // (*vv)(1) = 2;
  // (*vv)(2) = 3;
  // auto sv = std::make_shared<siconos::algebra::SparseVect>(3);
  // (*sv)(0) = 4;
  // (*sv)(1) = 5;
  // (*sv)(2) = 6;
  // auto sv2 = std::make_shared<siconos::algebra::SparseVect>(4);
  // (*sv2)(0) = 4;
  // (*sv2)(1) = 5;
  // (*sv2)(2) = 6;
  // auto w = std::make_shared<siconos::algebra::SiconosVector>(*sv);
  // auto ww = std::make_shared<siconos::algebra::SiconosVector>(*sv2);
  // auto res = std::make_shared<siconos::algebra::SiconosVector>(4);
  // auto res2 = std::make_shared<siconos::algebra::SiconosVector>(3);

  // // Triang * ...
  // *res2 = siconos::algebra::prod(*tmp2, *v);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *v->dense())) < tol, true);
  // *res2 = siconos::algebra::prod(*tmp2, *w);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *w->sparse())) < tol,
  //     true);
  // //   Sym * ...
  // *res2 = siconos::algebra::prod(*tmp3, *v);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *v->dense())) < tol, true);
  // *res2 = siconos::algebra::prod(*tmp3, *w);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *w->sparse())) < tol, true);
  // // Sparse * ...
  // *res = siconos::algebra::prod(*tmp4, *vv);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *vv->dense())) < tol, true);
  // *res = siconos::algebra::prod(*tmp4, *ww);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *ww->sparse())) < tol,
  //     true);
  // // Triang * ...
  // *res = siconos::algebra::prod(*tmp5, *v);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *v->dense())) < tol, true);
  // *res = siconos::algebra::prod(*tmp5, *w);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProd: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *w->sparse())) < tol, true);
  std::cout << "-->  test ublas::prod ended with success." << std::endl;
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
  yB->zero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  siconos::algebra::prod(*A, *x, *y);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  siconos::algebra::prod(*A, *xB, *y);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*y)(i)-sum) < tol, true);
  }

  // Block = Simple * Simple
  siconos::algebra::prod(*A, *x, *yB);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*yB)(i)-sum) < tol, true);
  }

  // Block = Simple * Block
  //  siconos::algebra::prod(*A ,*xB,*yB);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  //
  // Others or old stuff ...

  // auto tmp2 = std::make_shared<SimpleMatrix>(*T);
  // auto tmp3 = std::make_shared<SimpleMatrix>(*S);
  // auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp5 = std::make_shared<SimpleMatrix>(*Band2);
  // auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  // (*v)(0) = 1;
  // (*v)(1) = 2;
  // (*v)(2) = 3;
  // auto vv = std::make_shared<siconos::algebra::SiconosVector>(4);
  // (*vv)(0) = 1;
  // (*vv)(1) = 2;
  // (*vv)(2) = 3;
  // auto sv = std::make_shared<siconos::algebra::SparseVect>(3);
  // (*sv)(0) = 4;
  // (*sv)(1) = 5;
  // (*sv)(2) = 6;
  // auto sv2 = std::make_shared<siconos::algebra::SparseVect>(4);
  // (*sv2)(0) = 4;
  // (*sv2)(1) = 5;
  // (*sv2)(2) = 6;
  // auto w = std::make_shared<siconos::algebra::SiconosVector>(*sv);
  // auto ww = std::make_shared<siconos::algebra::SiconosVector>(*sv2);
  // auto res = std::make_shared<siconos::algebra::SiconosVector>(4);
  // auto res2 = std::make_shared<siconos::algebra::SiconosVector>(3);

  // // Triang * ...
  // siconos::algebra::prod(*tmp2, *v, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp2, *w, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *w->sparse())) < tol,
  //     true);
  // //   Sym * ...
  // siconos::algebra::prod(*tmp3, *v, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp3, *w, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *w->sparse())) < tol, true);
  // // Sparse * ...
  // siconos::algebra::prod(*tmp4, *vv, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *vv->dense())) < tol, true);
  // siconos::algebra::prod(*tmp4, *ww, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *ww->sparse())) < tol,
  //     true);
  // // Banded * ...
  // siconos::algebra::prod(*tmp5, *v, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp5, *w, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdBis: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *w->sparse())) < tol, true);
  std::cout << "-->  test ublas::prodBis ended with success." << std::endl;
}

void SimpleMatrixTest::testProdTer() {
  std::cout << "--> Test: ublas::prod. mat-vect (ter)" << std::endl;

  auto y = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto x = std::make_shared<siconos::algebra::SiconosVector>(size);
  x->setConstant(4.3);
  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(size - 2);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x2->setConstant(3.1);

  auto xB = std::make_shared<siconos::algebra::BlockVector>(x1, x2);
  auto yB = std::make_shared<siconos::algebra::BlockVector>(*xB);
  yB->zero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  // siconos::algebra::axpy_prod(*A, *x, *y, true);
  // double sum;
  // for (unsigned int i = 0; i < size; ++i)
  // {
  //   sum = 0;
  //   for (unsigned int j = 0; j < A->size(1); ++j)
  //     sum += (*A)(i, j) * (*x)(j);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum) < tol, true);
  // }

  // auto backUp= std::make_shared<siconos::algebra::SiconosVector>(*y);
  //  // Simple += Simple * Simple
  //  siconos::algebra::axpy_prod(*A, *x, *y, false);
  //  for (unsigned int i = 0; i < size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j = 0; j < A->size(1); ++j)
  //      sum += (*A)(i, j) * (*x)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum - (*backUp)(i)) < tol,
  //    true);
  //  }

  // Simple = Simple * Block
  //  siconos::algebra::axpy_prod(*A ,*xB,*y, true);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum)< tol, true);
  //  }
  //
  //*backUp = *y;
  // Simple += Simple * Block
  //  siconos::algebra::axpy_prod(*A ,*xB,*y, false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*y)(i) - sum - (*backUp)(i))< tol,
  //    true);
  //  }
  //
  //  // Block = Simple * Simple
  //  siconos::algebra::axpy_prod(*A ,*x,*yB, true);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*x)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  //
  //  // Block += Simple * Simple
  //  *backUp = *yB;
  //  siconos::algebra::axpy_prod(*A ,*x,*yB, false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*x)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum - (*backUp)(i))< tol,
  //    true);
  //  }
  //
  //  // Block = Simple * Block
  //  siconos::algebra::axpy_prod(*A ,*xB,*yB,true);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  //
  //  // Block += Simple * Block
  //  *backUp = *yB;
  //  siconos::algebra::axpy_prod(*A ,*xB,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += (*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdTer: ", fabs((*yB)(i) - sum - (*backUp)(i))< tol,
  //    true);
  //  }
  //  // Others or old stuff ...

  // auto tmp2 = std::make_shared<SimpleMatrix>(*T);
  // auto tmp3 = std::make_shared<SimpleMatrix>(*S);
  // auto tmp4 = std::make_shared<SimpleMatrix>(*SP);
  // auto tmp5 = std::make_shared<SimpleMatrix>(*Band2);
  // auto v = std::make_shared<siconos::algebra::SiconosVector>(3);
  // (*v)(0) = 1;
  // (*v)(1) = 2;
  // (*v)(2) = 3;
  // auto vv = std::make_shared<siconos::algebra::SiconosVector>(4);
  // (*vv)(0) = 1;
  // (*vv)(1) = 2;
  // (*vv)(2) = 3;
  // auto sv = std::make_shared<siconos::algebra::SparseVect>(3);
  // (*sv)(0) = 4;
  // (*sv)(1) = 5;
  // (*sv)(2) = 6;
  // auto sv2 = std::make_shared<siconos::algebra::SparseVect>(4);
  // (*sv2)(0) = 4;
  // (*sv2)(1) = 5;
  // (*sv2)(2) = 6;
  // auto w = std::make_shared<siconos::algebra::SiconosVector>(*sv);
  // auto ww = std::make_shared<siconos::algebra::SiconosVector>(*sv2);
  // auto res = std::make_shared<siconos::algebra::SiconosVector>(4);
  // auto res2 = std::make_shared<siconos::algebra::SiconosVector>(3);

  // // Triang * ...
  // siconos::algebra::prod(*tmp2, *v, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp2, *w, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp2->getTriang(), *w->sparse())) < tol,
  //     true);
  // //   Sym * ...
  // siconos::algebra::prod(*tmp3, *v, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp3, *w, *res2);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res2->dense() - ublas::prod(tmp3->getSym(), *w->sparse())) < tol, true);
  // // Sparse * ...
  // siconos::algebra::prod(*tmp4, *vv, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *vv->dense())) < tol, true);
  // siconos::algebra::prod(*tmp4, *ww, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp4->getSparse(), *ww->sparse())) < tol,
  //     true);
  // // Banded * ...
  // siconos::algebra::prod(*tmp5, *v, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *v->dense())) < tol, true);
  // siconos::algebra::prod(*tmp5, *w, *res);
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testProdTer: ",
  //     ublas::norm_2(*res->dense() - ublas::prod(tmp5->getBanded(), *w->sparse())) < tol, true);
  std::cout << "-->  test ublas::prodTer ended with success." << std::endl;
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
  yB->zero();

  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  y->zero();
  siconos::algebra::prod(*A, *x, *y, false);
  siconos::algebra::prod(*A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
  siconos::algebra::prod(*A, *xB, *y, false);
  siconos::algebra::prod(*A, *xB, *y, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * (*A)(i, j) * (*xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*y)(i)-sum) < tol, true);
  }

  // Block = Simple * Simple
  yB->zero();
  siconos::algebra::prod(*A, *x, *yB, false);
  siconos::algebra::prod(*A, *x, *yB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd4: ", fabs((*yB)(i)-sum) < tol, true);
  }

  // Block = Simple * Block
  yB->zero();
  //  siconos::algebra::prod(*A ,*xB,*yB,false);
  //  siconos::algebra::prod(*A ,*xB,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
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
  yB->zero();

  // Matrix - vector ublas::product
  double a = 3.0;
  // Simple = Simple * Simple
  y->zero();
  siconos::algebra::prod(a, *A, *x, *y, false);
  siconos::algebra::prod(a, *A, *x, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * a * (*A)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
  //  siconos::algebra::prod(a,*A ,*xB,*y,false);
  //  siconos::algebra::prod(a,*A ,*xB,*y,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += a*2*(*A)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*y)(i) - sum)< tol, true);
  //  }

  // Block = Simple * Simple
  yB->zero();
  //  siconos::algebra::prod(a,*A ,*x,*yB,false);
  //  siconos::algebra::prod(a,*A ,*x,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += a*2*(*A)(i,j)*(*x)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd5: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  //
  // Block = Simple * Block
  yB->zero();
  //  siconos::algebra::prod(a,*A ,*xB,*yB,false);
  //  siconos::algebra::prod(a,*A ,*xB,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
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
  yB->zero();

  auto tmp = std::make_shared<SimpleMatrix>(*A);
  tmp->transposeInPlace();
  // Matrix - vector ublas::product

  // Simple = Simple * Simple
  y->zero();
  siconos::algebra::prod(*x, *A, *y);
  siconos::algebra::prod(*x, *A, *y, false);
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*y)(i)-sum) < tol, true);
  }
  // Simple = Simple * Block
  y->zero();
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
  yB->zero();
  siconos::algebra::prod(*x, *A, *yB);
  siconos::algebra::prod(*x, *A, *yB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->size(1); ++j) sum += 2 * (*tmp)(i, j) * (*x)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i)-sum) < tol, true);
  }

  // Block = Simple * Block
  yB->zero();
  //  siconos::algebra::prod(*xB,*A ,*yB);
  //  siconos::algebra::prod(*xB,*A ,*yB,false);
  //  for (unsigned int i = 0; i< size; ++i)
  //  {
  //    sum = 0;
  //    for (unsigned int j=0; j< A->size(1); ++j)
  //      sum += 2*(*tmp)(i,j)*(*xB)(j);
  //    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd6: ", fabs((*yB)(i) - sum)< tol, true);
  //  }
  std::cout << "-->  test ublas::prod6 ended with success." << std::endl;
}

// void SimpleMatrixTest::testGemv()
// {
//   std::cout << "--> Test: gemv" <<std::endl;

//   auto y= std::make_shared<siconos::algebra::SiconosVector>(size, 1.0);
//   auto x= std::make_shared<siconos::algebra::SiconosVector>(size, 4.3);

//   auto backUp= std::make_shared<siconos::algebra::SiconosVector>(*y);

//   double a = 2.3;
//   double b = 1.5;
//   double sum;
//   gemv(a, *A, *x, b, *y);

//   for (unsigned int i = 0; i < size; ++i)
//   {
//     sum = b * (*backUp)(i);
//     for (unsigned int j = 0; j < A->size(1); ++j)
//       sum += a * (*A)(i, j) * (*x)(j) ;
//     CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv: ", fabs((*y)(i) - sum) < tol, true);
//   }

//   *y = *backUp;
//   gemvtranspose(a, *A, *x, b, *y);
//   for (unsigned int i = 0; i < size; ++i)
//   {
//     sum = b * (*backUp)(i);
//     for (unsigned int j = 0; j < A->size(0); ++j)
//       sum += a * (*A)(j, i) * (*x)(j);
//     CPPUNIT_ASSERT_EQUAL_MESSAGE("testgemv (trans): ", fabs((*y)(i) - sum) < tol, true);
//   }
//   std::cout << "-->  test gemv ended with success." <<std::endl;
// }

// void SimpleMatrixTest::testGemm()
// {
//   std::cout << "--> Test: gemm." <<std::endl;

//   double a = 2.3;
//   double b = 1.5;
//   *C = *A;
//   auto backUp = std::make_shared<SimpleMatrix>(*C);

//   gemm(a, *A, *B, b, *C);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm: ", norm_inf(*C->dense() - a *
//   ublas::prod(*A->dense(), *B->dense()) - b**backUp->dense()) < tol, true);

//   *C = *backUp;
//   gemmtranspose(a, *A, *B, b, *C);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testGemm (trans): ", norm_inf(*C->dense() - a *
//   siconos::algebra::prod(trans(*A->dense()), trans(*B->dense())) - b**backUp->dense()) <
//   tol, true); std::cout
//   << "-->  test gemm ended with success." <<std::endl;
// }

void SimpleMatrixTest::testFromAndFillCSC() {
//   std::cout << "Start SimpleMatrixTest::testFromAndFillCSC() " << std::endl;

//   auto Sparse4 = std::make_shared<SimpleMatrix>(*SP4);
//   Sparse4->updateNumericsMatrix();
//   auto NM = Sparse4->numericsMatrix();
//   NM_display(NM);
//   //  auto NM_1 = NM_create(4,4, NM_SPARSE);

//   auto Sparse1 = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
//   Sparse1->fromCSC(NM_csc(NM));
//   Sparse1->displayExpert();

//   auto NM_1 = NM_create(NM_SPARSE, 4, 4);
//   NM_1->matrix2->origin = NSM_CSC;
//   NM_csc_alloc(NM_1, Sparse4->nnz());
//   Sparse4->fillCSC(NM_csc(NM_1));
//   // NM_display(NM_1);  --> Note FP : fails when exiting the function ... To
//   // be investigating
//   // ...
//   NM_1 = NM_free(NM_1);
//   std::cout << "End SimpleMatrixTest::testFromAndFillCSC() " << std::endl;
// }

// void SimpleMatrixTest::testPLUFactorizationInPlace() {
//   std::cout << "--> Test: PLUFactorizationInPlace." << std::endl;

//   auto Dense = std::make_shared<SimpleMatrix>(*D);
//   Dense->display();
//   Dense->PLUFactorizationInPlace();
//   Dense->display();
//   // CPPUNIT_ASSERT_EQUAL_MESSAGE("testPLUFactorizationInPlace: ",  < tol, true);

//   auto Sparse = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
//   Sparse->eye();
//   Sparse->display();
//   Sparse->PLUFactorizationInPlace();
//   Sparse->display();

//   auto Sparse2 = std::make_shared<SimpleMatrix>(*SP4);
//   DEBUG_EXPR(Sparse2->display(););
//   Sparse2->PLUFactorizationInPlace();
//   DEBUG_EXPR(Sparse2->display(););

//   std::cout << "-->  test PLUFactorizationInPlace ended with success." << std::endl;
// }

// void SimpleMatrixTest::testFactorize() {
//   std::cout << "--> Test: Factorize (LU)." << std::endl;

//   auto Dense = std::make_shared<SimpleMatrix>(*D);
//   Dense->display();
//   Dense->Factorize();
//   Dense->display();
//   // CPPUNIT_ASSERT_EQUAL_MESSAGE("testFactorize: ",  < tol, true);

//   auto Sparse = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
//   Sparse->eye();
//   Sparse->display();
//   Sparse->Factorize();

//   Sparse = std::make_shared<SimpleMatrix>(*SP3);
//   Sparse->display();
//   Sparse->displayExpert();
//   Sparse->Factorize();
//   Sparse->display();

//   // Other types than SPARSE are not working since fillCSC is not implemented.
//   // std::cout << "--> Test: Factorize -- Triangle" <<std::endl;
//   // auto Triangle = std::make_shared<SimpleMatrix>(*T2);
//   // Triangle->display();
//   // Triangle->displayExpert();
//   // Triangle->Factorize();
//   // Triangle->display();

//   /* A last column full of zero caused memory corruption in cs_lusol
//      since ublas does not fill the last entry correctly*/
//   // Sparse = std::make_shared<SimpleMatrix>(*SP2);
//   // Sparse->display();
//   // Sparse->displayExpert();
//   // Sparse->PLUFactorizationInPlace();
//   // Sparse->display();

//   std::cout << "--> Test: Factorize (Cholesky)." << std::endl;

//   Dense = std::make_shared<SimpleMatrix>(*D);
//   /* conpute DD^T */
//   auto DenseT = std::make_shared<SimpleMatrix>(*Dense);
//   DenseT->trans();
//   auto DDT = std::make_shared<SimpleMatrix>(Dense->size(0), Dense->size(1));
//   siconos::algebra::prod(*Dense, *DenseT, *DDT);
//   DDT->display();
//   DDT->setIsSymmetric(true);
//   DDT->setIsPositiveDefinite(true);
//   DDT->Factorize();

//   DDT->display();
//   // CPPUNIT_ASSERT_EQUAL_MESSAGE("testFactorize: ",  < tol, true);

  std::cout << "-->  test Factorize ended with success." << std::endl;
}
void SimpleMatrixTest::testSolve() {
  std::cout << "\n--> Test: Solve. Dense. LU." << std::endl;

  // Test dense matrix
  auto Dense = std::make_shared<SimpleMatrix>(*D);
  auto b = std::make_shared<siconos::algebra::SiconosVector>(Dense->size(0));
  for (int i = 0; i < Dense->size(0); i++) {
    (*b)(i) = 1.0;
  }
  auto backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  auto D_backup = std::make_shared<SimpleMatrix>(*Dense);
  Dense->display();
  siconos::algebra::solveInPlace(*Dense, *b);
  Dense->display();
  b->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testSolve: ", (siconos::algebra::prod(*D_backup, *b) - *backup).norm() < tol, true);

  // // Test dense matrix and sparse rhs
  // Dense = std::make_shared<SimpleMatrix>(*D);
  // auto b_sparse = std::make_shared<siconos::algebra::SiconosVector>(Dense->size(0),
  // siconos::algebra::UblasType::SPARSE); for( int i =0; i <Dense->size(0); i++)
  // {
  //   (*b_sparse)(i)=1.0;
  // }
  // backup= std::make_shared<siconos::algebra::SiconosVector>(*b_sparse);
  // Dense->Solve(*b_sparse);
  // b_sparse->display();

  // CPPUNIT_ASSERT_EQUAL_MESSAGE("testSolve: ", (siconos::algebra::prod(*D_backup,*b_sparse) -
  // *backup).norm2() < tol, true);

  std::cout << "\n\n--> Test: Solve. Dense. Cholesky." << std::endl;
  Dense = std::make_shared<SimpleMatrix>(*D);
  b = std::make_shared<siconos::algebra::SiconosVector>(Dense->size(0));
  for (int i = 0; i < Dense->size(0); i++) {
    (*b)(i) = 1.0;
  }
  auto DenseT = std::make_shared<SimpleMatrix>(*Dense);
  DenseT->transposeInPlace();
  auto DDT = std::make_shared<SimpleMatrix>(Dense->size(0), Dense->size(1));
  siconos::algebra::prod(*Dense, *DenseT, *DDT);
  auto DDT_backup = std::make_shared<SimpleMatrix>(*DDT);
  // DDT->setIsSymmetric(true);
  // DDT->setIsPositiveDefinite(true);
  siconos::algebra::solveInPlace(*DDT, *b);
  b->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testSolve: ", (siconos::algebra::prod(*DDT_backup, *b) - *backup).norm() < tol, true);

  std::cout << "\n\n--> Test: Solve. Sparse. LU." << std::endl;

  // Test sparse matrix identity

  // auto Sparse = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // auto Sparse_backup =
  //     std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // Sparse_backup->eye();
  // Sparse->eye();
  // Sparse->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse_backup, *b) - *backup).norm2() < tol,
  //     true);

  // // test sparse matrix 3x3
  // Sparse = std::make_shared<SimpleMatrix>(*SP3);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // Sparse->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse, *b) - *backup).norm2() < tol, true);

  // // Solve again with another r.h.s.
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 2.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // Sparse->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse, *b) - *backup).norm2() < tol, true);

  // // test sparse matrix 4x4 SP4

  // Sparse = std::make_shared<SimpleMatrix>(*SP4);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // Sparse->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse, *b) - *backup).norm2() < tol, true);

  // std::cout << "\n\n--> Test: Solve. Sparse. LU.  Sparse rhs" << std::endl;

  // // test sparse matrix 3x3 sparse RhS. trivial solution (Id)
  // Sparse = std::make_shared<SimpleMatrix>(*SP3);
  // auto Sparse_rhs = std::make_shared<SimpleMatrix>(*SP3);
  // Sparse->Solve(*Sparse_rhs);
  // // std::cout << "Sparse_rhs :" << std::endl;
  // // Sparse_rhs->display();
  // // std::cout << "Sparse :" << std::endl;
  // // Sparse->display();
  // // std::cout << "A A^{-1}" << std::endl;
  // // (siconos::algebra::prod(*Sparse,*Sparse_rhs)).display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse, *Sparse_rhs) - *Sparse).normInf() < tol,
  //     true);

  // // test sparse matrix 3x3 sparse RhS. inverse
  // Sparse = std::make_shared<SimpleMatrix>(*SP3);
  // Sparse_rhs = std::make_shared<SimpleMatrix>(3, 3);
  // Sparse_rhs->eye();
  // Sparse->Solve(*Sparse_rhs);
  // auto Id = std::make_shared<SimpleMatrix>(3, 3);
  // Id->eye();

  // // Sparse_rhs->display();
  // // std::cout << "A A^{-1}" << std::endl;

  // // (siconos::algebra::prod(*Sparse,*Sparse_rhs)).display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse, *Sparse_rhs) - *Id).normInf() < tol,
  //     true);

  // std::cout << "\n\n--> Test: Solve. Sparse. Cholesky." << std::endl;

  // // Test sparse matrix identity
  // Sparse = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // Sparse_backup = std::make_shared<SimpleMatrix>(4, 4, siconos::algebra::UblasType::SPARSE);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // Sparse_backup->eye();
  // Sparse->eye();
  // Sparse->setIsSymmetric(true);
  // Sparse->setIsPositiveDefinite(true);
  // Sparse->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*Sparse_backup, *b) - *backup).norm2() < tol,
  //     true);

  // // test sparse matrix 3x3
  // Sparse = std::make_shared<SimpleMatrix>(*SP3);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // auto SparseT = std::make_shared<SimpleMatrix>(*SP3);
  // SparseT->trans();
  // auto SST = std::make_shared<SimpleMatrix>(Sparse->size(0), Sparse->size(1),
  //                                           siconos::algebra::UblasType::SPARSE);
  // siconos::algebra::prod(*Sparse, *SparseT, *SST);
  // SST->setIsSymmetric(true);
  // SST->setIsPositiveDefinite(true);
  // SST->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*SST, *b) - *backup).norm2() < tol, true);

  // // Solve again with another r.h.s.
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 2.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // SST->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*SST, *b) - *backup).norm2() < tol, true);

  // // test sparse matrix 4x4 SP4
  // Sparse = std::make_shared<SimpleMatrix>(*SP4);
  // b = std::make_shared<siconos::algebra::SiconosVector>(Sparse->size(0));
  // for (int i = 0; i < Sparse->size(0); i++) {
  //   (*b)(i) = 1.0;
  // }
  // backup = std::make_shared<siconos::algebra::SiconosVector>(*b);
  // SparseT = std::make_shared<SimpleMatrix>(*SP4);
  // SparseT->trans();
  // SST = std::make_shared<SimpleMatrix>(Sparse->size(0), Sparse->size(1),
  //                                      siconos::algebra::UblasType::SPARSE);
  // siconos::algebra::prod(*Sparse, *SparseT, *SST);
  // SST->setIsSymmetric(true);
  // SST->setIsPositiveDefinite(true);
  // SST->Solve(*b);
  // b->display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*SST, *b) - *backup).norm2() < tol, true);

  // std::cout << "\n\n--> Test: Solve. Sparse. Cholesky.  Sparse rhs" << std::endl;

  // // test sparse matrix 3x3 sparse RhS. trivial solution (Id)
  // Sparse = std::make_shared<SimpleMatrix>(*SP3);
  // SparseT = std::make_shared<SimpleMatrix>(*SP3);
  // SparseT->trans();
  // SST = std::make_shared<SimpleMatrix>(Sparse->size(0), Sparse->size(1),
  //                                      siconos::algebra::UblasType::SPARSE);
  // siconos::algebra::prod(*Sparse, *SparseT, *SST);
  // Sparse_rhs = std::make_shared<SimpleMatrix>(*SST);
  // // std::cout << "SST" << std::endl;
  // // SST->display();
  // SST->setIsSymmetric(true);
  // SST->setIsPositiveDefinite(true);
  // SST->Solve(*Sparse_rhs);
  // // Sparse_rhs->display();
  // // std::cout << "A A^{-1}" << std::endl;
  // // (siconos::algebra::prod(*SST,*Sparse_rhs)).display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*SST, *Sparse_rhs) - *SST).normInf() < tol, true);

  // // test sparse matrix 3x3 sparse RhS. inverse
  // SST = std::make_shared<SimpleMatrix>(Sparse->size(0), Sparse->size(1),
  //                                      siconos::algebra::UblasType::SPARSE);
  // siconos::algebra::prod(*Sparse, *SparseT, *SST);
  // Sparse_rhs = std::make_shared<SimpleMatrix>(3, 3, siconos::algebra::UblasType::SPARSE);
  // Sparse_rhs->eye();
  // // std::cout << "SST" << std::endl;
  // // SST->display();
  // SST->setIsSymmetric(true);
  // SST->setIsPositiveDefinite(true);
  // SST->Solve(*Sparse_rhs);
  // // Sparse_rhs->display();

  // // Sparse_rhs->display();
  // // std::cout << "A A^{-1}" << std::endl;

  // // (siconos::algebra::prod(*SST,*Sparse_rhs)).display();
  // CPPUNIT_ASSERT_EQUAL_MESSAGE(
  //     "testSolve: ", (siconos::algebra::prod(*SST, *Sparse_rhs) - *Id).normInf() < tol, true);

  std::cout << "-->  test Solve ended with success." << std::endl;
}

void SimpleMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SimpleMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
