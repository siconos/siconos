/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
#include "BlockMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(BlockMatrixTest);


void BlockMatrixTest::setUp()
{
  B = new SimpleMatrix(2, 2);
  C = new SimpleMatrix(2, 4);
  D = new SimpleMatrix(2, 1);
  E = new SimpleMatrix(3, 2);
  F = new SimpleMatrix(3, 4);
  G = new SimpleMatrix(3, 1);

  m.resize(6);
  m[0] = B ;
  m[1] = C ;
  m[2] = D ;
  m[3] = E ;
  m[4] = F ;
  m[5] = G ;
}

void BlockMatrixTest::tearDown()
{
  m.clear();
  delete B;
  delete C;
  delete D;
  delete E;
  delete F;
  delete G;
}

//______________________________________________________________________________

void BlockMatrixTest::testConstructor1()
{
  cout << "====================================" << endl;
  cout << "=== Block Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 1." << endl;
  BlockMatrix X(B, C, E, F);
  BlockMatrix * tmp = new BlockMatrix(X);


  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(tmp->getBlockPtr(0, 0)) == *B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(tmp->getBlockPtr(0, 1)) == *C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(tmp->getBlockPtr(1, 0)) == *E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *(tmp->getBlockPtr(1, 1)) == *F, true);
  delete tmp;
  cout << "--> Constructor 1 test ended with success." << endl;
}

void BlockMatrixTest::testConstructor2()
{
  cout << "--> Test: constructor 2." << endl;
  vector<LaGenMatDouble> mm(6);
  mm[0].resize(2, 2);
  mm[1].resize(2, 4);
  mm[2].resize(2, 1);
  mm[3].resize(3, 2);
  mm[4].resize(3, 4);
  mm[5].resize(3, 1);

  BlockMatrix X(m, 2, 3);
  vector<unsigned int> tmp;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", X.size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", X.size(1) == 7, true);
  tmp = X.getTabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp.size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp[1] == 5, true);
  tmp = X.getTabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp.size() == 3 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp[1] == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 ", tmp[2] == 7, true);

  cout << "--> Constructor 2 test ended with success." << endl;
}
void BlockMatrixTest::testConstructor3()
{
  cout << "--> Test: constructor 3." << endl;

  BlockMatrix X(m, 2, 3);
  vector<unsigned int> tmp;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.size(1) == 7, true);
  tmp = X.getTabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp.size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[1] == 5, true);
  tmp = X.getTabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp.size() == 3 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[1] == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[2] == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(1, 2) == G, true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

void BlockMatrixTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;

  BlockMatrix X(B, C, E, F);
  vector<unsigned int> tmp;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.size(0) == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.size(1) == 6, true);
  tmp = X.getTabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp.size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[1] == 5, true);
  tmp = X.getTabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp.size() == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[0] == 2 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", tmp[1] == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 ", X.getBlockPtr(1, 1) == F, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

void BlockMatrixTest::testGetLaGenMatDouble()
{
  cout << "--> Test: getLaGenMatDouble." << endl;
  BlockMatrix test(B, C, E, F);
  LaGenMatDouble mat = test.getLaGenMatDouble(0, 0);
  LaGenMatDouble mat1 = B->getLaGenMatDouble();
  double norm = Blas_Norm_Inf(mat - mat1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetLaGenMatDouble", norm < 1e-9 , true);
  cout << "--> getLaGenMatDouble test ended with success." << endl;
}

void BlockMatrixTest::testSetValue()
{
  cout << "--> Test: setValue." << endl;
  BlockMatrix test(B, C, E, F);
  LaGenMatDouble tmp(2, 2);
  tmp(0, 0) = 1;
  tmp(0, 1) = 2;
  tmp(1, 0) = 3;
  tmp(1, 1) = 4;
  test.setValue(tmp, 0, 0);
  LaGenMatDouble tmp2 = test.getLaGenMatDouble(0, 0);
  double norm = Blas_Norm_Inf(tmp - tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testSetValue", norm < 1e-9 , true);
  cout << "--> setValue test ended with success." << endl;
}

void BlockMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  BlockMatrix test(B, C, E, F);
  BlockMatrix test2(B, C, E, F);
  test2 *= 3;

  test = test2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", test.size(0) == test2.size(0), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", test.size(1) == test2.size(1), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", *(test.getBlockPtr(0, 0)) == 4 * *(test2.getBlockPtr(0, 0)), true);

  cout << "-->  test assignment ended with success." << endl;
}


void BlockMatrixTest::testOperators()
{
  cout << "--> Test: operators." << endl;
  //   BlockMatrix test(B,C,E,F);
  //   BlockMatrix test2(B,C,E,F);

  //   test+=test2;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators a:", *(test.getBlockPtr(0,0)) == 2**(test2.getBlockPtr(0,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators b:", *(test.getBlockPtr(0,1)) == 2**(test2.getBlockPtr(0,1)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators c:", *(test.getBlockPtr(1,0)) == 2**(test2.getBlockPtr(1,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators d:", *(test.getBlockPtr(1,1)) == 2**(test2.getBlockPtr(1,1)),true);
  //   cout << "--> Test: operators." << endl;

  //   test-=test2;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators e:", *(test.getBlockPtr(0,0)) == *(test2.getBlockPtr(0,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators f:", *(test.getBlockPtr(0,1)) == *(test2.getBlockPtr(0,1)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators g:", *(test.getBlockPtr(1,0)) == *(test2.getBlockPtr(1,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators h:", *(test.getBlockPtr(1,1)) == *(test2.getBlockPtr(1,1)),true);
  //   cout << "--> Test: operators." << endl;

  //   test*=3;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators e:", *(test.getBlockPtr(0,0)) == 3**(test2.getBlockPtr(0,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators f:", *(test.getBlockPtr(0,1)) == 3**(test2.getBlockPtr(0,1)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators g:", *(test.getBlockPtr(1,0)) == 3**(test2.getBlockPtr(1,0)),true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators h:", *(test.getBlockPtr(1,1)) == 3**(test2.getBlockPtr(1,1)),true);

  cout << "-->  test operators ended with success." << endl;
}


void BlockMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of BlockMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
