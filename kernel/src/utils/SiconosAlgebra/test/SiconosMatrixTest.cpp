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
#include "SiconosMatrixTest.hpp"

#include "NumericsToolsNamespace.h"  // for NM_csc, NM_free ...
#include "SiconosAlgebraAddons.hpp"  // for prod
#include "SiconosMatrix.hpp"
#include "Tools.hpp"
#include "io.hpp"
// namespace ublas = boost::numeric::ublas;

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMatrixTest);

#define DEBUG_MESSAGES
#include "siconos_debug.h"

using SiconosMatrix = siconos::algebra::SiconosDenseMatrix;

void SiconosMatrixTest::setUp() {
  tol = 1e-9;

  fic1 = "mat1.dat";  // 2 X 2
  fic2 = "mat2.dat";  // 2 X 3
  SicM = std::make_shared<SiconosMatrix>(siconos::algebra::io::readDenseMatrix(fic1));

  // BlockMat
  size = 10;

  A = std::make_shared<SiconosMatrix>(siconos::algebra::io::readDenseMatrix("A.dat"));

  // Build a reference sparse matrix
  int n = 6;
  std::vector<siconos::algebra::Triplet> triplets;
  triplets.emplace_back(0, 0, 1.0);
  triplets.emplace_back(0, 1, 2.0);
  triplets.emplace_back(1, 0, 2.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(2, 2, 4.0);
  triplets.emplace_back(2, 3, 5.0);
  triplets.emplace_back(3, 2, 7.0);
  triplets.emplace_back(3, 3, 6.0);
  Asparse.resize(n, n);
  Asparse.setFromTriplets(triplets.begin(), triplets.end());
}

void SiconosMatrixTest::tearDown() {}

void SiconosMatrixTest::testNormInf() {
  std::cout << "--> Test: normInf." << std::endl;
  double n = siconos::algebra::normInf(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7, true);

  auto nS = siconos::algebra::normInf(Asparse);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf (sparse): ", nS == 13, true);

  auto vecnorm = siconos::algebra::normInfByColumn(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf : ", vecnorm(0) == 3., true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf : ", vecnorm(1) == 4., true);

  auto vecnorm_sp = siconos::algebra::normInfByColumn(Asparse);
  siconos::algebra::SiconosVector ref{6};
  ref << 2., 3., 7., 6., 0., 0.;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf (sparse): ", vecnorm_sp.isApprox(ref, 1e-16),
                               true);
  std::cout << "--> normInf test ended with success." << std::endl;
}

void SiconosMatrixTest::testSymm() {
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSym : ", siconos::algebra::isSymmetric(Asparse), false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSym : ", siconos::algebra::isSymmetric(*SicM), false);

  siconos::algebra::SiconosSparseMatrix Bsym = Asparse;

  Bsym.coeffRef(3, 2) = 5.;
  Bsym.makeCompressed();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSym : ", siconos::algebra::isSymmetric(Bsym), true);

  siconos::algebra::SiconosDenseMatrix B{3, 3};
  B << 1, 2, 0, 2, 4, 3, 0, 3, 8;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSym : ", siconos::algebra::isSymmetric(B), true);
}

void SiconosMatrixTest::testfillTriplet() {
  CSparseMatrix* triplet = cs_spalloc(4, 4, 10, 1, 1);

  siconos::algebra::SiconosDenseMatrix dense(2, 2);
  dense(0, 0) = 1.0;
  dense(0, 1) = 0.0;
  dense(1, 0) = 3.0;
  dense(1, 1) = 4.0;
  siconos::algebra::fillTriplet(dense, triplet, 0, 0, 1e-12);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", triplet->nz == 3, true);

  bool found_1_0 = false, found_3_1 = false, found_4_1_1 = false;
  for (CS_INT k = 0; k < triplet->nz; ++k) {
    CS_INT i = triplet->i[k];
    CS_INT j = triplet->p[k];
    double x = triplet->x[k];
    if (i == 0 && j == 0 && fabs(x - 1.0) < 1e-12) found_1_0 = true;
    if (i == 1 && j == 0 && fabs(x - 3.0) < 1e-12) found_3_1 = true;
    if (i == 1 && j == 1 && fabs(x - 4.0) < 1e-12) found_4_1_1 = true;
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_1_0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_3_1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_4_1_1, true);

  cs_spfree(triplet);

  std::cout << "✅ test_fillTriplet passed.\n";

  siconos::algebra::SiconosSparseMatrix sparseMat(3, 3);
  sparseMat.insert(0, 0) = 1.5;
  sparseMat.insert(2, 1) = -3.0;
  sparseMat.insert(1, 2) = 2.0;
  sparseMat.makeCompressed();

  CSparseMatrix* triplet2 = cs_spalloc(5, 5, 10, 1, 1);
  assert(triplet2);

  siconos::algebra::fillTriplet(sparseMat, triplet2, 0, 0);

  // Vérifie que le nombre d'éléments non nuls est correct
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", triplet2->nz == 3, true);

  // Vérifie que les valeurs ont bien été copiées
  bool found_0_0 = false, found_2_1 = false, found_1_2 = false;
  for (CS_INT k = 0; k < triplet2->nz; ++k) {
    CS_INT i = triplet2->i[k];
    CS_INT j = triplet2->p[k];
    double x = triplet2->x[k];
    if (i == 0 && j == 0 && fabs(x - 1.5) < 1e-12) found_0_0 = true;
    if (i == 2 && j == 1 && fabs(x + 3.0) < 1e-12) found_2_1 = true;
    if (i == 1 && j == 2 && fabs(x - 2.0) < 1e-12) found_1_2 = true;
  }

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_0_0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_2_1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", found_1_2, true);

  std::cout << "✅ test_fillTriplet_sparse passed.\n";
}

void SiconosMatrixTest::testSetBlock() {
  std::cout << "--> Test: testSetBlock." << std::endl;

  // Copy of a sub-block of a Simple into a Simple
  SiconosMatrix result{10, 8};
  result.setRandom();
  SiconosMatrix mIn{2, 3};
  mIn << 1, 2, 3, 4, 5, 6;
  siconos::algebra::setBlock(1, 2, mIn, result);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", result.block(1, 2, 2, 3).isApprox(mIn, tol),
                               true);

  siconos::algebra::SiconosSparseMatrix resultSparse = Asparse;  // Copy
  siconos::algebra::setBlock(1, 2, mIn, resultSparse);

  siconos::algebra::SiconosDenseMatrix blockDense = resultSparse.block(1, 2, 2, 3).toDense();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetBlock: ", blockDense.isApprox(mIn, tol), true);

  std::cout << "✅-->  setBlock test ended with success." << std::endl;
}

void SiconosMatrixTest::testProd() {
  std::cout << "--> Test: mat-block vect product)" << std::endl;

  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(4);
  x1->setConstant(2.3);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(6);
  x2->setConstant(3.1);
  siconos::algebra::BlockVector xB{};
  xB.insertPtr(x1);
  xB.insertPtr(x2);

  siconos::algebra::SiconosVector result{10};
  siconos::algebra::matrixBlockVector_prod(*A, xB, result);

  double sum;
  for (Eigen::Index i = 0; i < size; ++i) {
    sum = 0;
    for (Eigen::Index j = 0; j < A->cols(); ++j) sum += (*A)(i, j) * (xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((result)(i)-sum) < tol, true);
  }

  siconos::algebra::matrixBlockVector_prod(*A, xB, result, false);

  for (Eigen::Index i = 0; i < size; ++i) {
    sum = 0;
    for (Eigen::Index j = 0; j < A->cols(); ++j) sum += 2. * (*A)(i, j) * (xB)(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", fabs((result)(i)-sum) < tol, true);
  }
  std::cout << "✅-->  test Test: mat-block vect ended with success." << std::endl;
}

void SiconosMatrixTest::testProd2()  // y += trans(A)*x
{
  std::cout << "--> Test: tA.x prod to block " << std::endl;

  siconos::algebra::SiconosVector x{10};
  x.setConstant(2.3);

  auto x1 = std::make_shared<siconos::algebra::SiconosVector>(4);
  auto x2 = std::make_shared<siconos::algebra::SiconosVector>(6);
  siconos::algebra::BlockVector resultB{};
  resultB.insertPtr(x1);
  resultB.insertPtr(x2);
  resultB.setZero();

  siconos::algebra::transposeMatrixVector_prod_toBlock(x, *A, resultB);

  auto tmp = std::make_shared<SiconosMatrix>(*A);
  tmp->transposeInPlace();
  double sum;
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += (*tmp)(i, j) * x(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd2: ", fabs(resultB(i) - sum) < tol, true);
  }

  siconos::algebra::transposeMatrixVector_prod_toBlock(x, *A, resultB, false);
  for (unsigned int i = 0; i < size; ++i) {
    sum = 0;
    for (unsigned int j = 0; j < A->cols(); ++j) sum += 2. * (*tmp)(i, j) * x(j);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd2: ", fabs(resultB(i) - sum) < tol, true);
  }
  std::cout << "✅-->  test prod2 ended with success." << std::endl;
}

void SiconosMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
