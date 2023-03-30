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
#include "EigenProblemsTest.hpp"

// #include <boost/numeric/bindings/std/vector.hpp>
// #include <boost/numeric/bindings/ublas/matrix.hpp>
// #include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <iostream>
#include <limits>

#include "SiconosAlgebraTools.hpp"
#include "SiconosConfig.h"
#include "SiconosMatrixOp.hpp"
#include "SiconosMatrixVectorOp.hpp"
#include "SiconosVector.hpp"
#include "bindings_utils.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// Note FP : add tests for complex matrices and geev, if needed (?)

CPPUNIT_TEST_SUITE_REGISTRATION(EigenProblemsTest);

namespace ublas = boost::numeric::ublas;

using complex_matrix = ublas::matrix<std::complex<double>, ublas::column_major>;
using complex_vector = ublas::vector<std::complex<double>>;

void EigenProblemsTest::setUp() {
  size = 5;
  A = std::make_shared<siconos::algebra::SimpleMatrix>(size, size);
  // Initialize A with random values.
  A->randomize();
  Aref = std::make_shared<siconos::algebra::SimpleMatrix>(*A);

  Asym = std::make_shared<siconos::algebra::SimpleMatrix>(
      size, size, siconos::algebra::UblasType::SYMMETRIC);
  Asym->randomize();
}

void EigenProblemsTest::tearDown() {}

void EigenProblemsTest::testSyev() {
  std::cout << "--> Test: syev." << std::endl;

  // turn A into a symmetric matrix
  Asym->randomize();
  *A = *Asym;
  *Aref = *Asym;

  // Initialize EigenVectors with A
  auto EigenValues = std::make_shared<siconos::algebra::SiconosVector>(size);
  auto EigenVectors = std::make_shared<siconos::algebra::SimpleMatrix>(*A);
  //  *EigenVectors = *A;

  siconos::algebra::syev(*EigenValues, *EigenVectors);

  siconos::algebra::DenseVect error(size);
  error *= 0.0;

  for (unsigned int i = 0; i < size; ++i) {
    error.plus_assign(ublas::prod(*A->dense(), ublas::column(*EigenVectors->dense(), i)));
    error.minus_assign((*EigenValues->dense())(i)*ublas::column(*EigenVectors->dense(), i));
  }
  // Check ...

  auto nerr = ublas::norm_2(error);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testSyev 1: ", nerr < 2000. * std::numeric_limits<double>::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 2: ", (*A) == (*Aref), true);

  // Now compute only eigenvalues
  auto RefEigenValues = std::make_shared<siconos::algebra::SiconosVector>(*EigenValues);
  *EigenVectors = *A;
  *EigenValues *= 0.0;
  siconos::algebra::syev(*EigenValues, *EigenVectors, false);
  nerr = ((*EigenValues) - (*RefEigenValues)).norm2();

  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testSyev 3: ", nerr < 2000 * std::numeric_limits<double>::epsilon(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 4: ", (*A) == (*Aref), true);
  std::cout << "--> Syev test ended with success." << std::endl;
}

void EigenProblemsTest::testGeev1() {
  std::cout << "--> Test: geev1." << std::endl;
  // Compute only right eigenvectors.
  complex_matrix fake(1, 1), rightV(size, size);
  complex_vector eigenval(size);
  siconos::algebra::geev(*A, eigenval, fake, rightV);
  complex_vector error(size);
  for (unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for (unsigned int i = 0; i < size; ++i) {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i)));
    error.minus_assign(eigenval(i) * column(rightV, i));
  }

  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testGeev1 1: ", ublas::norm_2(error) < 10000 * std::numeric_limits<double>::epsilon(),
      true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 2: ", (*A) == (*Aref), true);
  // Now compute only eigenvalues
  complex_vector RefEigenValues(size);
  RefEigenValues = eigenval;
  eigenval *= 0.0;
  std::cout << ublas::norm_2(eigenval - RefEigenValues) << "\n";
  siconos::algebra::geev(*A, eigenval, fake, fake, false, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testGeev1 3: ",
      ublas::norm_2(eigenval - RefEigenValues) < 10 * std::numeric_limits<double>::epsilon(),
      true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 4: ", (*A) == (*Aref), true);

  std::cout << "--> geev1 test ended with success." << std::endl;
}

void EigenProblemsTest::testGeev2() {
  std::cout << "--> Test: geev2." << std::endl;
  // Compute only left eigenvectors.
  complex_matrix fake(1, 1), leftV(size, size);
  complex_vector eigenval(size);
  siconos::algebra::geev(*A, eigenval, leftV, fake, true, false);
  complex_vector error(size);
  for (unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for (unsigned int i = 0; i < size; ++i) {
    error.plus_assign(ublas::prod(conj(column(leftV, i)), *A->dense()));
    error.minus_assign(eigenval(i) * conj(column(leftV, i)));
  }
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testGeev2 1: ", ublas::norm_2(error) < 3000 * std::numeric_limits<double>::epsilon(),
      true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev2 2: ", (*A) == (*Aref), true);

  std::cout << "--> geev1 test ended with success." << std::endl;
}

void EigenProblemsTest::testGeev3() {
  std::cout << "--> Test: geev3." << std::endl;

  // Compute left and right eigenvectors.
  complex_matrix leftV(size, size), rightV(size, size);
  complex_vector eigenval(size);
  siconos::algebra::geev(*A, eigenval, leftV, rightV, true, true);
  complex_vector error(size);
  for (unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for (unsigned int i = 0; i < size; ++i) {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i)));
    error.minus_assign(eigenval(i) * column(rightV, i));
    error.plus_assign(ublas::prod(conj(column(leftV, i)), *A->dense()));
    error.minus_assign(eigenval(i) * conj(column(leftV, i)));
  }
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE(
      "testGeev3 1: ",
      ublas::norm_2(error) < size * 2000 * std::numeric_limits<double>::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev3 2: ", (*A) == (*Aref), true);

  std::cout << "--> geev3 test ended with success." << std::endl;
}

void EigenProblemsTest::End() {
  std::cout << "======================================" << std::endl;
  std::cout << " ===== End of EigenProblems tests ===== " << std::endl;
  std::cout << "======================================" << std::endl;
}
