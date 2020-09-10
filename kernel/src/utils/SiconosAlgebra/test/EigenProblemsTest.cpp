/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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
#include "SiconosConfig.h"

#include "EigenProblemsTest.hpp"
#include "SiconosAlgebra.hpp"
#include "SimpleMatrixFriends.hpp"
#include "bindings_utils.hpp"
#include <limits>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);


// Note FP : add tests for complex matrices and geev, if needed (?)

CPPUNIT_TEST_SUITE_REGISTRATION(EigenProblemsTest);

void EigenProblemsTest::setUp()
{
  size = 5;
  A.reset(new SimpleMatrix(size,size));
  // Initialize A with random values.
  A->randomize();
  Aref.reset(new SimpleMatrix(*A));
}

void EigenProblemsTest::tearDown()
{}

void EigenProblemsTest::testSyev()
{
  std::cout << "--> Test: syev." <<std::endl;

  // turn A to a symmetric matrix
  A->randomize_sym();
  *Aref = *A;

  // Initialize EigenVectors with A
  SP::SiconosVector EigenValues(new SiconosVector(size));
  SP::SimpleMatrix EigenVectors(new SimpleMatrix(*A));
//  *EigenVectors = *A;

  Siconos::eigenproblems::syev(*EigenValues, *EigenVectors);

  DenseVect error(size);
  error *= 0.0;

  for(unsigned int i = 0; i < size; ++i)
  {
    error.plus_assign(ublas::prod(*A->dense(), column(*EigenVectors->dense(), i)));
    error.minus_assign((*EigenValues->dense())(i)*column(*EigenVectors->dense(),i));
  }
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 2: ", (*A) == (*Aref), true);

  // Now compute only eigenvalues
  SP::SiconosVector RefEigenValues(new SiconosVector(*EigenValues));
  *EigenVectors = *A;
  *EigenValues *= 0.0;
  Siconos::eigenproblems::syev(*EigenValues, *EigenVectors, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 3: ", ((*EigenValues) - (*RefEigenValues)).norm2() < 10 * std::numeric_limits< double >::epsilon(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 4: ", (*A) == (*Aref), true);
  std::cout << "--> Syev test ended with success." <<std::endl;
}

void EigenProblemsTest::testGeev1()
{
  std::cout << "--> Test: geev1." <<std::endl;
  // Compute only right eigenvectors.
  complex_matrix fake(1,1), rightV(size,size);
  complex_vector eigenval(size);
  Siconos::eigenproblems::geev(*A, eigenval, fake, rightV);
  complex_vector error(size);
  for(unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for(unsigned int i = 0; i < size; ++i)
  {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i)));
    error.minus_assign(eigenval(i)*column(rightV,i));
  }

  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 2: ", (*A) == (*Aref), true);
  // Now compute only eigenvalues
  complex_vector RefEigenValues(size);
  RefEigenValues = eigenval;
  eigenval *= 0.0;
  Siconos::eigenproblems::geev(*A, eigenval, fake, fake, false, false);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 3: ", norm_2(eigenval - RefEigenValues) < 10 * std::numeric_limits< double >::epsilon(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 4: ", (*A) == (*Aref), true);

  std::cout << "--> geev1 test ended with success." <<std::endl;
}

void EigenProblemsTest::testGeev2()
{
  std::cout << "--> Test: geev2." <<std::endl;
  // Compute only left eigenvectors.
  complex_matrix fake(1,1), leftV(size,size);
  complex_vector eigenval(size);
  Siconos::eigenproblems::geev(*A, eigenval, leftV, fake, true, false);
  complex_vector error(size);
  for(unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for(unsigned int i = 0; i < size; ++i)
  {
    error.plus_assign(ublas::prod(conj(column(leftV, i)), *A->dense()));
    error.minus_assign(eigenval(i)*conj(column(leftV,i)));
  }
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev2 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev2 2: ", (*A) == (*Aref), true);

  std::cout << "--> geev1 test ended with success." <<std::endl;
}

void EigenProblemsTest::testGeev3()
{
  std::cout << "--> Test: geev3." <<std::endl;

  // Compute left and right eigenvectors.
  complex_matrix leftV(size,size), rightV(size,size);
  complex_vector eigenval(size);
  Siconos::eigenproblems::geev(*A, eigenval, leftV, rightV, true, true);
  complex_vector error(size);
  for(unsigned int i = 0; i < size; ++i) error(i) = 0.0;
  for(unsigned int i = 0; i < size; ++i)
  {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i)));
    error.minus_assign(eigenval(i)*column(rightV,i));
    error.plus_assign(ublas::prod(conj(column(leftV, i)), *A->dense()));
    error.minus_assign(eigenval(i)*conj(column(leftV,i)));
  }
  std::cout << norm_2(error) << std::endl;
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev3 1: ", norm_2(error) < size * 10 * std::numeric_limits< double >::epsilon(), true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev3 2: ", (*A) == (*Aref), true);

  std::cout << "--> geev3 test ended with success." <<std::endl;
}

void EigenProblemsTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of EigenProblems tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
