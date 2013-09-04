/* Siconos-Kernel, Copyright INRIA 2005-2013.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "KernelConfig.h"

#include "EigenProblemsTest.hpp"
#include "SiconosAlgebra.hpp"
#include "bindings_utils.hpp"
#include <limits>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "SiconosVector.hpp"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega) \
  if ((alpha) == (omega)) CPPUNIT_FAIL(message);

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
  // Initialize EigenVectors with A
  SP::SiconosVector EigenValues(new SiconosVector(size));
  SP::SimpleMatrix EigenVectors(new SimpleMatrix(*A));
//  *EigenVectors = *A;
  
  Siconos::eigenproblems::syev(*EigenValues, *EigenVectors);
  
  DenseVect error(size);
  error *= 0.0;

  for( unsigned int i = 0; i < size; ++i )
  {
    error.plus_assign(ublas::prod( *A->dense(), column(*EigenVectors->dense(), i) ));
    error.minus_assign((*EigenValues->dense())(i)*column(*EigenVectors->dense(),i));
  }
  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon() , true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 2: ", (*A) == (*Aref) , true);
  
  // Now compute only eigenvalues 
  SP::SiconosVector RefEigenValues(new SiconosVector(*EigenValues));
  *EigenVectors = *A;
  *EigenValues *= 0.0;
  Siconos::eigenproblems::syev(*EigenValues, *EigenVectors, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 3: ", ((*EigenValues) - (*RefEigenValues)).norm2() < 10 * std::numeric_limits< double >::epsilon(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSyev 4: ", (*A) == (*Aref) , true);
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
  error *= 0.0;
  for( unsigned int i = 0; i < size; ++i )
  {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i) ));
    error.minus_assign(eigenval(i)*column(rightV,i));
  }

  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon() , true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 2: ", (*A) == (*Aref) , true);
  // Now compute only eigenvalues 
  complex_vector RefEigenValues(size);
  RefEigenValues = eigenval;
  eigenval *= 0.0;
  Siconos::eigenproblems::geev(*A, eigenval, fake, fake, false, false);
  
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 3: ", norm_2(eigenval - RefEigenValues) < 10 * std::numeric_limits< double >::epsilon() , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev1 4: ", (*A) == (*Aref) , true);

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
  error *= 0.0;
  for( unsigned int i = 0; i < size; ++i )
  {
    error.plus_assign(ublas::prod(column(leftV, i), *A->dense() ));
    error.minus_assign(eigenval(i)*column(leftV,i));
  }

  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev2 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon() , true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev2 2: ", (*A) == (*Aref) , true);

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
  error *= 0.0;
  for( unsigned int i = 0; i < size; ++i )
  {
    error.plus_assign(ublas::prod(*A->dense(), column(rightV, i) ));
    error.minus_assign(eigenval(i)*column(rightV,i));
    error.plus_assign(ublas::prod(column(leftV, i), *A->dense() ));
    error.minus_assign(eigenval(i)*column(leftV,i));
  }

  // Check ...
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev3 1: ", norm_2(error) < 10 * std::numeric_limits< double >::epsilon() , true);
  // Check if A has not been modified
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGeev3 2: ", (*A) == (*Aref) , true);

  std::cout << "--> geev3 test ended with success." <<std::endl;
}

void EigenProblemsTest::End()
{
  std::cout << "======================================" <<std::endl;
  std::cout << " ===== End of EigenProblems tests ===== " <<std::endl;
  std::cout << "======================================" <<std::endl;
}
