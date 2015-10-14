/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
#include "NonSmoothSolverTest.hpp"


#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(NonSmoothSolverTest);


void NonSmoothSolverTest::setUp()
{
  iparam.resize(NB_PARAM);
  dparam.resize(NB_PARAM);

  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    iparam[i] = i + 1;
    dparam[i] = 2.0 * i + 1.2;
  }
}

void NonSmoothSolverTest::tearDown()
{}

// Default constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver0()
{
  std::cout << "===========================================" <<std::endl;
  std::cout << " ===== NonSmoothSolver tests start ...===== " <<std::endl;
  std::cout << "===========================================" <<std::endl;
  std::cout << "------- Default constructor test -------" <<std::endl;
  SP::NonSmoothSolver NSS(new NonSmoothSolver());
  std::vector<int> i0 = *NSS->intParameters();
  std::vector<double> d0 = *NSS->doubleParameters();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->name() == "undefined", true);
  SP::SolverOptions opt = NSS->numericsSolverOptions();
  int * i1 = opt->iparam;
  double * d1 = opt->dparam;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverE : ", opt->isSet == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverG : ", opt->nbParam == (int)NB_PARAM, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverH : ", !i1, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverI : ", !d1, false);
  std::cout << "------- Default constructor NonSmoothSolver ok -------" <<std::endl;
}

// Basic constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver1()
{
  std::cout << "------- Basic constructor test -------" <<std::endl;
  SP::NonSmoothSolver NSS(new NonSmoothSolver("NSGS", iparam, dparam));
  std::vector<int> i0 = *NSS->intParameters();
  std::vector<double> d0 = *NSS->doubleParameters();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->name() == "NSGS", true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverC : ", i0[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverD : ", d0[i] == dparam[i], true);
  }
  SP::SolverOptions opt = NSS->numericsSolverOptions();
  int * i1 = opt->iparam;
  double * d1 = opt->dparam;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverE : ", opt->isSet == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverF : ", strcmp(opt->solverName, "NSGS") == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverG : ", opt->nbParam == (int)NB_PARAM, true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverH : ", i1[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverI : ", d1[i] == dparam[i], true);
  }
  std::cout << "------- Basic Constructor NonSmoothSolver ok -------" <<std::endl;
}

// Copy constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver2()
{
  std::cout << "------- Copy constructor test -------" <<std::endl;
  SP::NonSmoothSolver NSS0(new NonSmoothSolver("NSGS", iparam, dparam));
  SP::NonSmoothSolver NSS(new NonSmoothSolver(*NSS0));

  std::vector<int> i0 = *NSS->intParameters();
  std::vector<double> d0 = *NSS->doubleParameters();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->name() == "NSGS", true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverC : ", i0[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverD : ", d0[i] == dparam[i], true);
  }
  SP::SolverOptions opt = NSS->numericsSolverOptions();
  int * i1 = opt->iparam;
  double * d1 = opt->dparam;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverE : ", opt->isSet == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverF : ", strcmp(opt->solverName, "NSGS") == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverG : ", opt->nbParam == (int)NB_PARAM, true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverH : ", i1[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverI : ", d1[i] == dparam[i], true);
  }
  std::cout << "------- Constructor copy NonSmoothSolver ok -------" <<std::endl;
}

void NonSmoothSolverTest::End()
{
  std::cout << "==========================================" <<std::endl;
  std::cout << " ===== End of NonSmoothSolver tests ===== " <<std::endl;
  std::cout << "==========================================" <<std::endl;
}
