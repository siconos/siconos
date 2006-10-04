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
#include "MyBlockMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(MyBlockMatrixTest);


void MyBlockMatrixTest::setUp()
{
  B = new MySimpleMatrix(2, 2);
  C = new MySimpleMatrix(2, 4);
  D = new MySimpleMatrix(2, 1);
  E = new MySimpleMatrix(3, 2);
  F = new MySimpleMatrix(3, 4);
  G = new MySimpleMatrix(3, 1);

  m.resize(6);
  m[0] = B ;
  m[1] = C ;
  m[2] = D ;
  m[3] = E ;
  m[4] = F ;
  m[5] = G ;

  tRow.resize(2);
  tRow[0] = 2;
  tRow[1] = 5;
  tCol.resize(3);
  tCol[0] = 2;
  tCol[1] = 6;
  tCol[2] = 7;

  mapRef = new mapped(2, 1);
  (*mapRef)(0, 0) = B;
  (*mapRef)(1, 0) = E;

}

void MyBlockMatrixTest::tearDown()
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

void MyBlockMatrixTest::testConstructor0() // constructor with a vector of MySiconosMatrix*
{
  cout << "====================================" << endl;
  cout << "=== Block Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 0." << endl;
  MyBlockMatrix * test = new MyBlockMatrix(m, 2, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size1() == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size2() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(0, 2) == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(1, 1) == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getBlockPtr(1, 2) == G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getTabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getTabCol() == tCol, true);
  delete test;
  cout << "--> Constructor 0 test ended with success." << endl;
}

void MyBlockMatrixTest::testConstructor1() // Copy constructor, from a MyBlockMatrix
{
  cout << "--> Test: constructor 1." << endl;
  MyBlockMatrix * ref = new MyBlockMatrix(m, 2, 3);
  MyBlockMatrix * test = new MyBlockMatrix(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size1() == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->size2() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(0, 2) == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(1, 1) == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getBlockPtr(1, 2) == G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getTabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test->getTabCol() == tCol, true);
  cout << "--> Constructor (copy) test ended with success." << endl;
}

void MyBlockMatrixTest::testConstructor2() // Copy constructor, from a MySiconosMatrix
{
  cout << "--> Test: constructor 2." << endl;
  MySiconosMatrix * ref = new MyBlockMatrix(m, 2, 3);
  MyBlockMatrix * test = new MyBlockMatrix(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size1() == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->size2() == 7, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(0, 2) == D, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(1, 1) == F, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getBlockPtr(1, 2) == G, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getTabRow() == tRow, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test->getTabCol() == tCol, true);
  cout << "--> Constructor (copy) test ended with success." << endl;
}

void MyBlockMatrixTest::testConstructor3() // Constructor from a mapped
{
  cout << "--> Test: constructor 3." << endl;
  MyBlockMatrix * test = new MyBlockMatrix(*mapRef);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size1() == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->size2() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", test->getBlockPtr(1, 0) == E, true);
  std::vector<unsigned int> tmp = test->getTabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[1] == 3, true);
  tmp.clear();
  tmp = test->getTabCol();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp.size() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[0] == 2, true);
  delete test;
  cout << "--> Constructor 3 test ended with success." << endl;
}

void MyBlockMatrixTest::testConstructor4() // Constructor from 4 MySiconosMatrix*
{
  cout << "--> Test: constructor 4." << endl;
  MyBlockMatrix * test = new MyBlockMatrix(B, C, E, F);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->isBlock() == true, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size1() == 5 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->size2() == 6, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getBlockPtr(0, 0) == B, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getBlockPtr(0, 1) == C, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getBlockPtr(1, 0) == E, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getBlockPtr(1, 1) == F, true);
  std::vector<unsigned int> tmp = test->getTabRow();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[1] == 5, true);
  tmp.clear();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp.size() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[0] == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ", tmp[1] == 6, true);
  delete test;
  cout << "--> Constructor 4 test ended with success." << endl;
}

// void MyBlockMatrixTest::testResize()
// {
//   cout << "--> Test: resize." << endl;
//   MySiconosMatrix * test = new MyBlockMatrix(m,2,3);
//   test->resize(3,4);
//   (*test)(2,0) = B;
//   (*test)(2,1) = C;
//   (*test)(2,2) = D;
//   (*test)(2,3) = C;
//   (*test)(0,3) = C;
//   (*test)(1,3) = F;

//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->isBlock() == true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size1() == 7, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->size2() == 11, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(0,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(0,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(0,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(0,3) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(1,0) == E, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(1,1) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(1,2) == G, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(1,3) == F, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(2,0) == B, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(2,1) == C, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(2,2) == D, true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ",test->getBlockPtr(2,3) == C, true);

//   delete test;
//   cout << "--> resize test ended with success." << endl;
// }

void MyBlockMatrixTest::testNormInf()
{
  cout << "--> Test: normInf." << endl;
  MySiconosMatrix * test = new MyBlockMatrix(m, 2, 3);
  double n = 12;
  (*test)(4, 3) = n;
  (*test)(2, 1) = n - 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", test->normInf() == n , true);
  cout << "--> normInf test ended with success." << endl;
}

void MyBlockMatrixTest::testZero()
{
  cout << "--> Test: zero." << endl;
  MySiconosMatrix* A = new MySimpleMatrix(2, 2);
  A->eye();
  MySiconosMatrix* H = new MySimpleMatrix(2, 4);
  H->eye();
  MySiconosMatrix* I = new MySimpleMatrix(5, 2);
  I->eye();
  MySiconosMatrix* J = new MySimpleMatrix(5, 4);
  J->eye();

  std::vector<MySiconosMatrix*> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  MySiconosMatrix * test = new MyBlockMatrix(v, 2, 2);
  test->zero();
  unsigned int n1 = test->size1();
  unsigned int n2 = test->size2();
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*test)(i, j) == 0, true);
  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*A)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*J)(i, j) == 0, true);
  }

  delete test;
  delete A;
  delete H;
  delete I;
  delete J;
  cout << "--> zero test ended with success." << endl;
}

void MyBlockMatrixTest::testEye()
{
  cout << "--> Test: eye." << endl;
  MySiconosMatrix* A = new MySimpleMatrix(2, 2);
  MySiconosMatrix* H = new MySimpleMatrix(2, 4);
  MySiconosMatrix* I = new MySimpleMatrix(5, 2);
  MySiconosMatrix* J = new MySimpleMatrix(5, 4);

  std::vector<MySiconosMatrix*> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  MySiconosMatrix * test = new MyBlockMatrix(v, 2, 2);
  test->eye();
  unsigned int n1 = test->size1();
  unsigned int n2 = test->size2();
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*test)(i, j) == 0, true);

  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*A)(i, j) == 0, true);
    }
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*H)(i, j) == 0, true);
  }
  for (unsigned int i = 0; i < 5; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*I)(i, j) == 0, true);
    for (unsigned int j = 0; j < 4; ++j)
    {
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 1, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*J)(i, j) == 0, true);
    }
  }

  delete test;
  delete A;
  delete H;
  delete I;
  delete J;
  cout << "--> eye test ended with success." << endl;
}
// Add tests with getDense ...

void MyBlockMatrixTest::testGetSetRowCol()
{
  cout << "--> Test: get, set Row and Col." << endl;
  MySimpleVector * tmp = new MySimpleVector(6);
  MySimpleVector * tmp1 = new MySimpleVector(6);
  (*tmp1)(0) = 1;
  (*tmp1)(2) = 2;
  MySiconosMatrix* A = new MySimpleMatrix(2, 2);
  MySiconosMatrix* H = new MySimpleMatrix(2, 4);
  MySiconosMatrix* I = new MySimpleMatrix(5, 2);
  MySiconosMatrix* J = new MySimpleMatrix(5, 4);
  std::vector<MySiconosMatrix*> v(4);
  v[0] = A ;
  v[1] = H ;
  v[2] = I ;
  v[3] = J ;
  MySiconosMatrix * test = new MyBlockMatrix(v, 2, 2);

  test->setRow(1, *tmp1);
  test->getRow(1, *tmp);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp == *tmp1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(1, 0) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*H)(1, 0) == 2, true);

  MySimpleVector * tmp2 = new MySimpleVector(7);
  MySimpleVector * tmp3 = new MySimpleVector(7);
  (*tmp3)(0) = 1;
  (*tmp3)(2) = 2;
  test->setCol(1, *tmp3);
  test->getCol(1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp2 == *tmp3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*A)(0, 1) == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", (*I)(0, 1) == 2, true);


  delete test;
  delete A;
  delete H;
  delete I;
  delete J;
  delete tmp;
  delete tmp1;
  delete tmp2;
  delete tmp3;
  cout << "--> get, set Row and Col tests ended with success." << endl;
}

void MyBlockMatrixTest::testGetBlock()
{
  cout << "--> Test: getBlock." << endl;
  MySiconosMatrix * test = new MyBlockMatrix(m, 2, 3);
  MySiconosMatrix * block = new MySimpleMatrix(3, 1);
  (*G)(0, 0) = 1;
  (*G)(1, 0) = 12;
  (*G)(2, 0) = 4;

  test->getBlock(2, 2, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(0, 0) == 1 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(1, 0) == 12 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(2, 0) == 4 , true);
  test->zero();
  delete block;
  delete test;
  cout << "--> getBlock test ended with success." << endl;
}

void MyBlockMatrixTest::testMatrixCopy()
{
  cout << "--> Test: matrixCopy." << endl;
  MySiconosMatrix * test = new MyBlockMatrix(m, 2, 3);
  MySiconosMatrix * block = new MySimpleMatrix(3, 1);
  (*block)(0, 0) = 1;
  (*block)(1, 0) = 2;
  (*block)(2, 0) = 3;
  test->matrixCopy(*block, 2, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *G == *block , true);

  delete block;
  delete test;
  cout << "-->  matrixCopy test ended with success." << endl;
}

void MyBlockMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  MySiconosMatrix * test = new MyBlockMatrix(m, 2, 3);
  MySiconosMatrix * ref = new MySimpleMatrix(5, 7);
  for (unsigned int i = 0; i < 5 ; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      (*ref)(i, j) = i + j;

  *test = *ref;
  for (unsigned int i = 0; i < 2; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*B)(i, j) == i + j , true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*C)(i, j) == i + j + 2 , true);
    for (unsigned int j = 0; j < 1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*D)(i, j) == i + j + 6 , true);
  }
  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j < 2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*E)(i, j) == 2 + i + j , true);
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*F)(i, j) == i + j + 4 , true);
    for (unsigned int j = 0; j < 1; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*G)(i, j) == i + j + 8 , true);
  }

  ref = new MyBlockMatrix(m, 2, 3);
  *test = *ref;
  for (unsigned int i = 0; i < 5; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*test)(i, j) == 0, true);
  cout << "-->  test assignment ended with success." << endl;
}

void MyBlockMatrixTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;
  MySiconosMatrix * tmp = new MyBlockMatrix(m, 2, 3);
  MySiconosMatrix * ref = new MySimpleMatrix(5, 7);
  for (unsigned int i = 0; i < 5 ; ++i)
    for (unsigned int j = 0; j < 7; ++j)
      (*ref)(i, j) = i + j;
  tmp->zero();
  // operators with a MySimpleMatrix
  *tmp += *ref;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*ref)(i, j) , true);

  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*ref)(i, j)) , true);

  int mult = 2;
  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 4.4 * ((*ref)(i, j)) , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*ref)(i, j)) , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*ref)(i, j) , true);

  *tmp *= 3;
  *tmp -= *ref;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*ref)(i, j) , true);

  // operators with a MyBlockMatrix

  MySiconosMatrix * ref2 = new MyBlockMatrix(m, 3, 2);
  *ref2 = *ref; // ref2 is initialized with ref

  tmp->zero();
  *tmp += *ref2;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*ref2)(i, j) , true);

  mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*ref2)(i, j)) , true);

  mult = 2;
  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 4.4 * ((*ref2)(i, j)) , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*ref2)(i, j)) , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*ref2)(i, j) , true);

  *tmp *= 3;
  *tmp -= *ref2;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*ref2)(i, j) , true);
  delete ref2;
  delete ref;
  delete tmp;
  cout << "-->  test operators1 ended with success." << endl;
}

void MyBlockMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of MyBlockMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
