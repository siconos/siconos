/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
#include "NonSmoothSolverTest.h"
using namespace std;

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
  cout << "===========================================" << endl;
  cout << " ===== NonSmoothSolver tests start ...===== " << endl;
  cout << "===========================================" << endl;
  cout << "------- Default constructor test -------" << endl;
  NonSmoothSolver * NSS = new NonSmoothSolver();
  std::vector<int> i0 = *NSS->getIntParametersPtr();
  std::vector<double> d0 = *NSS->getDoubleParametersPtr();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->getName() == "undefined", true);
  Solver_Options * opt = NSS->getNumericsSolverOptionsPtr();
  int * i1 = opt->iparam;
  double * d1 = opt->dparam;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverE : ", opt->isSet == 0, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverG : ", opt->nbParam == (int)NB_PARAM, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverH : ", i1 == NULL, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverI : ", d1 == NULL, true);

  delete NSS;
  cout << "------- Default constructor NonSmoothSolver ok -------" << endl;
}

// Basic constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver1()
{
  cout << "------- Basic constructor test -------" << endl;
  NonSmoothSolver * NSS = new NonSmoothSolver("NSGS", iparam, dparam);
  std::vector<int> i0 = *NSS->getIntParametersPtr();
  std::vector<double> d0 = *NSS->getDoubleParametersPtr();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->getName() == "NSGS", true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverC : ", i0[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverD : ", d0[i] == dparam[i], true);
  }
  Solver_Options * opt = NSS->getNumericsSolverOptionsPtr();
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
  delete NSS;
  cout << "------- Basic Constructor NonSmoothSolver ok -------" << endl;
}

// Copy constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver2()
{
  cout << "------- Copy constructor test -------" << endl;
  NonSmoothSolver * NSS0 = new NonSmoothSolver("NSGS", iparam, dparam);
  NonSmoothSolver * NSS = new NonSmoothSolver(*NSS0);

  std::vector<int> i0 = *NSS->getIntParametersPtr();
  std::vector<double> d0 = *NSS->getDoubleParametersPtr();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->getName() == "NSGS", true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverC : ", i0[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverD : ", d0[i] == dparam[i], true);
  }
  Solver_Options * opt = NSS->getNumericsSolverOptionsPtr();
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
  delete NSS;
  cout << "------- Constructor copy NonSmoothSolver ok -------" << endl;
}

// XML constructor
void NonSmoothSolverTest::testBuildNonSmoothSolver3()
{
  cout << "------- XML constructor test -------" << endl;

  // parse xml file:
  xmlDocPtr doc = xmlParseFile("NonSmoothSolverTest.xml");
  if (doc == NULL)
    XMLException::selfThrow("Document not parsed successfully");

  xmlNodePtr cur = xmlDocGetRootElement(doc);
  if (cur == NULL)
  {
    XMLException::selfThrow("empty document");
    xmlFreeDoc(doc);
  }

  // get rootNode
  if (xmlStrcmp(cur->name, (const xmlChar *) "SiconosModel"))
  {
    XMLException::selfThrow("document of the wrong type, root node !=SiconosModel");
    xmlFreeDoc(doc);
  }

  // look for NonSmoothSolver node
  xmlNodePtr node = SiconosDOMTreeTools::findNodeChild(cur, "NonSmoothSolver");
  NonSmoothSolverXML * solvXML = new NonSmoothSolverXML(node);
  NonSmoothSolver * NSS = new NonSmoothSolver(solvXML);
  std::vector<int> i0 = *NSS->getIntParametersPtr();
  std::vector<double> d0 = *NSS->getDoubleParametersPtr();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverA : ", NSS->isSolverSet(), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverB : ", NSS->getName() == "NSGS", true);
  for (unsigned int i = 0; i < NB_PARAM; ++i)
  {
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverC : ", i0[i] == iparam[i], true);
    CPPUNIT_ASSERT_EQUAL_MESSAGE("testBuildNonSmoothSolverD : ", d0[i] == dparam[i], true);
  }
  Solver_Options * opt = NSS->getNumericsSolverOptionsPtr();
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
  delete NSS;
  delete solvXML;
  cout << "------- XML Constructor NonSmoothSolver ok -------" << endl;
}

void NonSmoothSolverTest::End()
{
  cout << "==========================================" << endl;
  cout << " ===== End of NonSmoothSolver tests ===== " << endl;
  cout << "==========================================" << endl;
}
