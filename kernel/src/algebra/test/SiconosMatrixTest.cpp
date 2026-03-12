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

#include "SiconosAlgebraAddons.hpp"  // for prod
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "io.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMatrixTest);

#define DEBUG_MESSAGES

using SiconosMatrix = siconos::algebra::SiconosDenseMatrix;

void SiconosMatrixTest::setUp() {
  tol = 1e-9;

  fic1 = "mat1.dat";  // 2 X 2
  fic2 = "mat2.dat";  // 2 X 3
  SicM = std::make_shared<SiconosMatrix>(siconos::algebra::io::readDenseMatrix(fic1));

  // BlockMat
  size = 10;

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

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testfillTriplet : ", triplet2->nz == 3, true);

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

// How to build a SiconosSparseMatrix
void SiconosMatrixTest::testBuildSparse() {
  std::cout << "Start testBuildSparse \n";
  // 1. Random values.
  // Input = dimension and number of nonzeros values
  siconos::algebra::Index rows = 20;
  siconos::algebra::Index cols = 20;
  int nnz = 150;
  auto m = siconos::algebra::generateRandomSparseMatrix(rows, cols, nnz);
  std::cout << "nnz" << m.nonZeros() << "\n";
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build random - check size: ", m.rows() == rows, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build random - check size: ", m.cols() == cols, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build random - check nnz: ", m.nonZeros() == nnz, true);

  std::cout << "✅ testBuildSparse - Build random sparse matrix passed.\n";
  // 2. Explicitely set values (from triplets)
  typedef Eigen::Triplet<double> T;
  std::vector<siconos::algebra::Triplet> triplets;
  rows = 6;
  cols = 6;
  triplets.emplace_back(0, 0, 1.0);
  triplets.emplace_back(0, 1, 2.0);
  triplets.emplace_back(1, 0, 2.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(2, 2, 4.0);
  triplets.emplace_back(2, 3, 5.0);
  triplets.emplace_back(3, 2, 7.0);
  triplets.emplace_back(3, 3, 6.0);
  siconos::algebra::SiconosSparseMatrix m2(rows, cols);
  m2.setFromTriplets(triplets.begin(), triplets.end());
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build from triplet- check size: ", m2.rows() == rows, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build from triplet- check size: ", m2.cols() == cols, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build from triplet- check nnz: ", m2.nonZeros() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build from triplet- check coeff: ", m2.coeffRef(1, 1) == 3.,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build from triplet- check coeff: ", m2.coeffRef(3, 2) == 7.,
                               true);
  std::cout << "✅ testBuildSparse - Build sparse matrix from triplets passed.\n";

  siconos::algebra::SiconosSparseMatrix m3(rows, cols);
  m3.reserve(Eigen::VectorXi::Constant(cols, 6));
  m3.insert(0, 0) = 1.0;
  m3.insert(0, 1) = 2.0;
  m3.insert(1, 0) = 2.0;
  m3.insert(1, 1) = 3.0;
  m3.insert(2, 2) = 4.0;
  m3.insert(2, 3) = 5.0;
  m3.insert(3, 2) = 7.0;
  m3.insert(3, 3) = 6.0;
  m3.makeCompressed();  // optional
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build with insert- check size: ", m3.rows() == rows, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build with insert- check size: ", m3.cols() == cols, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build with insert- check nnz: ", m3.nonZeros() == 8, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build with insert- check coeff: ", m3.coeffRef(1, 1) == 3.,
                               true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("Build with insert- check coeff: ", m3.coeffRef(3, 2) == 7.,
                               true);
  std::cout << "✅ testBuildSparse - Build sparse matrix with inserts passed.\n";
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
  auto A = siconos::algebra::io::readDenseMatrix("A.dat");

  siconos::algebra::SiconosVector result{10};
  siconos::algebra::matrixBlockVector_prod(A, xB, result);

  auto xref = xB.toSiconosVector();
  auto ref = A * xref;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", result.isApprox(ref, tol), true);

  siconos::algebra::matrixBlockVector_prod(A, xB, result, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", result.isApprox(2. * ref, tol), true);

  // Call with a view
  siconos::algebra::MapType Aview = siconos::algebra::MapType(A.data(), 10, 10);
  siconos::algebra::matrixBlockVector_prod(Aview, xB, result);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProdBis: ", result.isApprox(ref, tol), true);

  std::cout << "✅ Test: mat(dense)-block vect ended with success." << std::endl;

  // Now sparse version
  std::vector<siconos::algebra::Triplet> triplets;
  triplets.emplace_back(0, 0, 1.0);
  triplets.emplace_back(0, 1, 2.0);
  triplets.emplace_back(1, 0, 2.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(2, 2, 4.0);
  triplets.emplace_back(2, 3, 5.0);
  triplets.emplace_back(3, 2, 7.0);
  triplets.emplace_back(3, 3, 6.0);
  siconos::algebra::SiconosSparseMatrix sparse{7, 6};
  sparse.setFromTriplets(triplets.begin(), triplets.end());

  auto x3 = std::make_shared<siconos::algebra::SiconosVector>(4);
  x3->setRandom();
  auto x4 = std::make_shared<siconos::algebra::SiconosVector>(2);
  x4->setRandom();
  siconos::algebra::BlockVector yB{};
  yB.insertPtr(x3);
  yB.insertPtr(x4);

  siconos::algebra::SiconosVector resultsparse{7};
  resultsparse.setZero();
  siconos::algebra::matrixBlockVector_prod(sparse, yB, resultsparse);
  auto ycopy = yB.toSiconosVector();
  siconos::algebra::SiconosDenseMatrix dense{sparse};
  auto ref2 = dense * ycopy;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test mat-block ", resultsparse.isApprox(ref2, tol), true);

  // +=
  siconos::algebra::matrixBlockVector_prod(sparse, yB, resultsparse, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test mat-block: ", resultsparse.isApprox(2. * ref2, tol),
                               true);
  std::cout << "✅ Test: mat(sparse)-block vect ended with success." << std::endl;
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
  auto A = siconos::algebra::io::readDenseMatrix("A.dat");

  siconos::algebra::transposeMatrixVector_prod_toBlock(x, A, resultB);

  auto ref = A.transpose() * x;
  auto res = resultB.toSiconosVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd2: ", res.isApprox(ref, tol), true);

  siconos::algebra::transposeMatrixVector_prod_toBlock(x, A, resultB, false);
  res = resultB.toSiconosVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd2: ", res.isApprox(2. * ref, tol), true);

  // with view
  siconos::algebra::MapType Aview = siconos::algebra::MapType(A.data(), 10, 10);
  siconos::algebra::transposeMatrixVector_prod_toBlock(x, Aview, resultB);
  res = resultB.toSiconosVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testProd2: ", res.isApprox(ref, tol), true);

  std::cout << "✅-->  test mat(dense)-vec prod (to block) ended with success." << std::endl;
  // Now, sparse
  auto sparse = siconos::algebra::generateRandomSparseMatrix(x.size(), 8, 12);
  siconos::algebra::BlockVector yB{2, 4};
  yB.setZero();

  siconos::algebra::transposeMatrixVector_prod_toBlock(x, sparse, yB);
  auto ycopy = yB.toSiconosVector();
  auto ref2 = sparse.transpose() * x;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test mat-block ", ycopy.isApprox(ref2, tol), true);

  // +=
  siconos::algebra::transposeMatrixVector_prod_toBlock(x, sparse, yB, false);
  ycopy = yB.toSiconosVector();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("test mat-block ", ycopy.isApprox(2. * ref2, tol), true);

  std::cout << "✅-->  test mat(sparse)-vec prod (to block) ended with success." << std::endl;
}

void SiconosMatrixTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of SiconosMatrix Tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
