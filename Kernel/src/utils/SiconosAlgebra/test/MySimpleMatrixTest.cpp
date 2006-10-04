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
#include "MySimpleMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(MySimpleMatrixTest);


void MySimpleMatrixTest::setUp()
{

  fic1 = "mat1.dat";
  fic2 = "mat2.dat";
  SicM = new MySimpleMatrix(fic1, 1);
  SimM = new MySimpleMatrix(fic2, 1);
  std::vector<double> v3(2, 0);
  std::vector<double> v4(2, 0);
  std::vector<double> v5(3, 0);
  v4[0] = 6;
  v4[1] = 9;
  v5[0] = 8;
  v5[1] = 9;
  v5[2] = 10;

  vect1 = new MySimpleVector(v3, 2);
  vect2 = new MySimpleVector(v4, 2); // vect2 != vect1, but vect2 == SimM second column
  vect3 = new MySimpleVector(v5, 3); // vect3 != vect1, but vect3 == SimM second row

  // Dense
  D = new DenseMat(2, 2);
  for (unsigned i = 0; i < D->size1(); ++ i)
    for (unsigned j = 0; j < D->size2(); ++ j)
      (*D)(i, j) = 3 * i + j;
  // Triang
  T = new TriangMat(3, 3);
  for (unsigned i = 0; i < T->size1(); ++ i)
    for (unsigned j = i; j < T->size2(); ++ j)
      (*T)(i, j) = 3 * i + j;
  // Sym
  S = new SymMat(3, 3);
  for (unsigned i = 0; i < S->size1(); ++ i)
    for (unsigned j = i; j < S->size2(); ++ j)
      (*S)(i, j) = 3 * i + j;
  // Sparse
  SP = new SparseMat(4, 4);
  for (unsigned i = 0; i < SP->size1(); ++ i)
    for (unsigned j = 0; j < SP->size2(); ++ j)
      (*SP)(i, j) = 3 * i + j;

  // Banded
  B = new BandedMat(4, 4, 1, 1);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      (*B)(i, j) = 3 * i + j;
}

void MySimpleMatrixTest::tearDown()
{
  delete vect1;
  delete vect2;
  delete SicM;
  delete SimM;
}

//______________________________________________________________________________

void MySimpleMatrixTest::testConstructor0() // constructor with TYP and dim
{
  cout << "====================================" << endl;
  cout << "=== Simple Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 0." << endl;
  MySimpleMatrix * test = new MySimpleMatrix(2, 3);

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->size2() == 3, true);
  delete test;
  cout << "--> Constructor 0 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor1() // Copy constructor, from a MySimpleMatrix
{
  cout << "--> Test: constructor 1." << endl;
  MySimpleMatrix test(*SimM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test == *SimM, true);
  cout << "--> Constructor (copy) test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor2() // Copy constructor, from a MySiconosMatrix
{
  cout << "--> Test: constructor 2." << endl;
  MySimpleMatrix *  test = new MySimpleMatrix(*SicM);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", *test == *SicM, true);
  cout << "--> Constructor (copy) test ended with success." << endl;
}

// void MySimpleMatrixTest::testConstructor3()
// {
//   cout << "--> Test: constructor 3." << endl;
//   MySimpleMatrix * test = new MySimpleMatrix(SPARSE);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor3 : ",test->getNum() == 4, true);
//   delete test;
//   cout << "--> Constructor 3 test ended with success." << endl;
// }

void MySimpleMatrixTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(*D);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", test->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : ", norm_inf(test->getDense() - *D) == 0, true);
  delete test;
  cout << "--> Constructor 4 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(*T);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", test->getNum() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", norm_inf(test->getTriang() - *T) == 0, true);
  delete test;
  cout << "--> Constructor 5 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor6()
{
  cout << "--> Test: constructor 6." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(*S);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", test->getNum() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor6 : ", norm_inf(test->getSym() - *S) == 0, true);
  delete test;
  cout << "--> Constructor 6 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor7()
{
  cout << "--> Test: constructor 7." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(*SP);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", test->getNum() == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor7 : ", norm_inf(test->getSparse() - *SP) == 0, true);
  delete test;
  cout << "--> Constructor 7 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor8()
{
  cout << "--> Test: constructor 8." << endl;
  cout << "--> Constructor 8 test ended with success." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(*B);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor8 : ", test->getNum() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor8 : ", norm_inf(test->getBanded() - *B) == 0, true);
  delete test;
}

void MySimpleMatrixTest::testConstructor9()
{
  cout << "--> Test: constructor 9." << endl;
  const std::vector<double> v1(4, 1);
  const std::vector<double> v(6, 1);
  //DENSE
  MySiconosMatrix *test1 = new MySimpleMatrix(v1, 2, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->size2() == 2, true);
  //TRIANGULAR
  MySiconosMatrix *test2 = new MySimpleMatrix(v1, 2, 2, TRIANGULAR);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->getNum() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->size2() == 2, true);
  //SYMMETRIC
  MySiconosMatrix *test3 = new MySimpleMatrix(v1, 2, 2, SYMMETRIC);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->getNum() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->size2() == 2, true);
  //BANDED
  MySiconosMatrix *test4 = new MySimpleMatrix(v, 2, 2, BANDED, 1, 1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test4->getNum() == 5, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test4->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test4->size2() == 2, true);
  delete test4;
  delete test3;
  delete test2;
  delete test1;
  cout << "--> Constructor 9 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor10()
{
  cout << "--> Test: constructor 10." << endl;
  MySiconosMatrix * test = new MySimpleMatrix(fic1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor10 : ", *test == *SicM, true);
  delete test;
  cout << "--> Constructor 10 test ended with success." << endl;
}

// Add tests with getDense ...

void MySimpleMatrixTest::testGetSetRowCol()
{
  cout << "--> Test: get, set Row and Col." << endl;
  MySimpleVector * tmp = new MySimpleVector(3);
  SimM->getRow(1, *tmp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp == *vect3, true);

  MySimpleVector * tmp2 = new MySimpleVector(2);
  SimM->getCol(1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp2 == *vect2, true);

  MySiconosMatrix * M = new MySimpleMatrix(*SimM);
  M->setRow(0, *vect3);
  tmp->zero();
  M->getRow(0, *tmp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testGetSetRowCol: ", *tmp == *vect3, true);
  M->setCol(0, *vect1);
  M->getCol(0, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp2 == *vect1, true);
  delete M;
  delete tmp2;
  delete tmp;
  cout << "--> get, set Row and Col tests ended with success." << endl;
}

void MySimpleMatrixTest::testZero()
{
  cout << "--> Test: zero." << endl;
  MySiconosMatrix * tmp = new MySimpleMatrix(*SimM);
  tmp->zero();
  unsigned int n1 = tmp->size1();
  unsigned int n2 = tmp->size2();
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testZero : ", (*tmp)(i, j) == 0, true);
  cout << "--> zero test ended with success." << endl;
}

void MySimpleMatrixTest::testEye()
{
  cout << "--> Test: eye." << endl;
  MySiconosMatrix * tmp = new MySimpleMatrix(*SimM);
  tmp->eye();
  unsigned int n1 = tmp->size1();
  unsigned int n2 = tmp->size2();
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      if (i != j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 0, true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testEye : ", (*tmp)(i, j) == 1, true);
  cout << "--> eye test ended with success." << endl;
}

void MySimpleMatrixTest::testResize()
{
  cout << "--> Test: resize." << endl;
  MySiconosMatrix * tmp = new MySimpleMatrix(*SicM);
  tmp->resize(3, 4);
  unsigned int n1 = SicM->size1();
  unsigned int n2 = SicM->size2();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size1() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size2() == 4, true);

  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i, j) == (*SicM)(i, j), true);
  for (unsigned int i = n1; i < 3; ++i)
    for (unsigned int j = 0; j < 4; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i, j) == 0, true);
  for (unsigned int j = n2; j < 4; ++j)
    for (unsigned int i = 0; i < 3; ++i)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i, j) == 0, true)
      ;
  // Check the effect of bool = false (ie preserve == false in boost resize)
  //   tmp->resize(6,8, false);
  //   tmp->display();
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size1() == 6, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size2() == 8, true);
  //   for(unsigned int i = 0; i<6; ++i)
  //     for(unsigned int j=0;j<8;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  //   // Reduction ...
  //   tmp->resize(1,2, false);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size1() == 1, true);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", tmp->size2() == 2, true);
  //   for(unsigned int i = 0; i<1; ++i)
  //     for(unsigned int j=0;j<2;++j)
  //       CPPUNIT_ASSERT_EQUAL_MESSAGE("testResize : ", (*tmp)(i,j) == 0 , true);
  delete tmp;
  cout << "--> resize test ended with success." << endl;
}

void MySimpleMatrixTest::testNormInf()
{
  cout << "--> Test: normInf." << endl;
  double n = SicM->normInf();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testNormInf: ", n == 7 , true);
  cout << "--> normInf test ended with success." << endl;
}

void MySimpleMatrixTest::testGetBlock()
{
  cout << "--> Test: getBlock." << endl;
  MySiconosMatrix * full = new MySimpleMatrix(*SimM);
  full->resize(4, 5);
  MySiconosMatrix * block = new MySimpleMatrix(2, 2);
  full->getBlock(1, 1, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(0, 0) == 9 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(0, 1) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(1, 0) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(1, 1) == 0 , true);
  delete block;
  delete full;
  cout << "--> getBlock test ended with success." << endl;
}

void MySimpleMatrixTest::testMatrixCopy()
{
  cout << "--> Test: blockMatrixCopy." << endl;
  // 1 -- Copy into a dense matrix of a ...
  // Dense
  MySiconosMatrix * full = new MySimpleMatrix(*SimM);
  full->resize(4, 5);
  full->matrixCopy(*SicM, 2, 3);
  MySiconosMatrix * block = new MySimpleMatrix(2, 2);
  full->getBlock(2, 3, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *block == *SicM , true);
  delete block;
  // Triang
  MySiconosMatrix * Tmp = new MySimpleMatrix(*T);
  full->matrixCopy(*Tmp, 1, 2);
  block = new MySimpleMatrix(3, 3);
  full->getBlock(1, 2, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *block == *Tmp , true);
  delete block;
  delete Tmp;
  // Sym
  Tmp = new MySimpleMatrix(*S);
  full->matrixCopy(*Tmp, 1, 2);
  block = new MySimpleMatrix(3, 3);
  full->getBlock(1, 2, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *block == *Tmp , true);
  delete block;
  delete Tmp;
  // Sparse
  Tmp = new MySimpleMatrix(*SP);
  full->matrixCopy(*Tmp, 0, 1);
  block = new MySimpleMatrix(4, 4);
  full->getBlock(0, 1, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *block == *Tmp , true);
  delete block;
  delete Tmp;
  // Banded
  Tmp = new MySimpleMatrix(*B);
  full->matrixCopy(*Tmp, 0, 1);
  block = new MySimpleMatrix(4, 4);
  full->getBlock(0, 1, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testMatrixCopy: ", *block == *Tmp , true);
  delete block;
  delete Tmp;
  delete full;
  cout << "-->  matrixCopy test ended with success." << endl;
}

void MySimpleMatrixTest::testTrans()
{
  cout << "--> Test: trans." << endl;
  // Dense
  MySiconosMatrix * ref = new MySimpleMatrix(*D);
  MySiconosMatrix * tRef = new MySimpleMatrix(*ref);
  *tRef = trans(*ref);
  for (unsigned int i = 0; i < ref->size1(); ++i)
    for (unsigned int j = 0 ; j < ref->size2(); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i) , true);

  delete tRef;
  delete ref;
  // Sym
  ref = new MySimpleMatrix(*S);
  tRef = new MySimpleMatrix(*ref);
  *tRef = trans(*ref);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef) == (*ref) , true);
  delete tRef;
  delete ref;
  // Sparse
  ref = new MySimpleMatrix(*SP);
  tRef = new MySimpleMatrix(*ref);
  *tRef = trans(*ref);
  for (unsigned int i = 0; i < ref->size1(); ++i)
    for (unsigned int j = 0 ; j < ref->size2(); ++j)
      if (i == j)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i, j) == (*ref)(j, i) , true);

  delete tRef;
  delete ref;
  // Banded
  //   ref = new MySimpleMatrix(*B);
  //   tRef = new MySimpleMatrix(*ref);
  //   *tRef = trans(*ref);
  //   for(unsigned int i = 0; i<ref->size1(); ++i)
  //     for(unsigned int j = 0 ; j< ref->size2(); ++j)
  //       if(i==j)
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(i,j) , true);
  //       else
  //  CPPUNIT_ASSERT_EQUAL_MESSAGE("testTrans: ", (*tRef)(i,j) == (*ref)(j,i) , true);
  //   delete tRef;
  //   delete ref;
  cout << "-->  test trans ended with success." << endl;
}

void MySimpleMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  MySiconosMatrix * ref = new MySimpleMatrix(*D);
  MySiconosMatrix * tRef = new MySimpleMatrix(*SicM);
  // Dense = any type:
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new MySimpleMatrix(*T);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new MySimpleMatrix(*S);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new MySimpleMatrix(*SP);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;

  ref = new MySimpleMatrix(*B);
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Triang = Triang
  ref = new MySimpleMatrix(*T);
  tRef = new MySimpleMatrix(*T);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Sym = Sym
  ref = new MySimpleMatrix(*S);
  tRef = new MySimpleMatrix(*S);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Sparse = Sparse
  ref = new MySimpleMatrix(*SP);
  tRef = new MySimpleMatrix(*SP);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  // Banded = Banded
  ref = new MySimpleMatrix(*B);
  tRef = new MySimpleMatrix(*B);
  tRef->zero();
  *tRef = *ref;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment: ", (*tRef) == (*ref) , true);
  delete ref;
  delete tRef;
  cout << "-->  test assignment ended with success." << endl;
}

void MySimpleMatrixTest::testOperators1()
{
  cout << "--> Test: operators1." << endl;
  //+=, -=, *=, /= for a dense matrix
  MySiconosMatrix * tmp = new MySimpleMatrix(*D);
  *tmp += *SicM;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SicM)(i, j) + (*D)(i, j) , true);

  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*SicM)(i, j) + (*D)(i, j)) , true);

  int mult = 2;
  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 4.4 * ((*SicM)(i, j) + (*D)(i, j)) , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.2 * ((*SicM)(i, j) + (*D)(i, j)) , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == ((*SicM)(i, j) + (*D)(i, j)) , true);

  *tmp -= *SicM;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getDense() - *D) == 0 , true);
  delete tmp;
  cout << "-->  test operators1 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators2()
{
  cout << "--> Test: operators2." << endl;
  double tol = 1e-12;
  // +=, -=, *=, /= triangular
  MySiconosMatrix * tmp = new MySimpleMatrix(*T);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*T);
  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*T)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*T)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*T)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*T)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getTriang() - *T) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 2 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators2 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators3()
{
  cout << "--> Test: operators3." << endl;
  double tol = 1e-12;
  // +=, -=, *=, /= triangular
  MySiconosMatrix * tmp = new MySimpleMatrix(*S);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*S);
  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*S)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*S)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*S)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*S)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSym() - *S) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 3 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators3 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators4()
{
  cout << "--> Test: operators4." << endl;
  double tol = 1e-12;
  // +=, -=, *=, /= triangular
  MySiconosMatrix * tmp = new MySimpleMatrix(*SP);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*SP);
  MySiconosMatrix * tmp3 = new MySimpleMatrix(*T);
  tmp3->resize(4, 4);
  MySiconosMatrix * tmp4 = new MySimpleMatrix(*B);
  MySiconosMatrix * tmp5 = new MySimpleMatrix(*S);
  tmp5->resize(4, 4);
  *tmp += *tmp2;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*SP)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol , true);

  *tmp *= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*SP)(i, j)) < tol , true);

  *tmp /= mult;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*SP)(i, j)) < tol , true);

  *tmp /= mult0;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0 ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*SP)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getSparse() - *SP) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 4 , true);

  // += -= a triangular
  *tmp += *tmp3;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
  {
    for (unsigned int j = 0; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) , true);
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) + (*tmp3)(i, j) , true);
  }

  *tmp -= *tmp3;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) , true);

  // += -= a banded
  *tmp += *tmp4;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*B)(i, j) + (*SP)(i, j) , true);

  *tmp -= *tmp4;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) , true);

  // += -= a sym
  *tmp += *tmp5;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
  {
    for (unsigned int j = 0; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) + (*tmp5)(j, i) , true);
    for (unsigned int j = i ; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) + (*tmp5)(i, j) , true);
  }

  *tmp -= *tmp5;
  for (unsigned int i = 0; i < tmp->size1(); ++i)
    for (unsigned int j = 0; j < tmp->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == (*SP)(i, j) , true);


  delete tmp;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  cout << "-->  test operators4 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators5()
{
  cout << "--> Test: operators5." << endl;
  double tol = 1e-12;
  // +=, -=, *=, /= triangular
  MySiconosMatrix * tmp = new MySimpleMatrix(*B);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*B);
  *tmp += *tmp2;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2.0 * (*B)(i, j) , true);

  int mult = 2;
  double mult0 = 2.2;
  *tmp *= mult0;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*B)(i, j)) < tol , true);

  *tmp *= mult;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult * mult0 * (*B)(i, j)) < tol , true);

  *tmp /= mult;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", ((*tmp)(i, j) - 2.0 * mult0 * (*B)(i, j)) < tol , true);

  *tmp /= mult0;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*tmp)(i, j) == 2 * (*B)(i, j) , true);

  *tmp -= *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(tmp->getBanded() - *B) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", tmp->getNum() == 5 , true);

  delete tmp;
  delete tmp2;
  cout << "-->  test operators5 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators6()
{
  cout << "--> Test: operator6." << endl;

  // Dense +,-,* Dense
  MySiconosMatrix * tmp = new MySimpleMatrix(*D);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*SicM);
  MySiconosMatrix * res = new MySimpleMatrix(2, 2);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*SicM)(i, j) + (*D)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*D)(i, j) - (*SicM)(i, j)), true);

  *res = *tmp * *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(0, 0) == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(0, 1) == 4, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(1, 0) == 15, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(1, 1) == 22, true);
  delete tmp;
  delete tmp2;
  delete res;

  // Triang +,-,* Triang
  tmp = new MySimpleMatrix(*T);
  tmp2 = new MySimpleMatrix(*T);
  res = new MySimpleMatrix(3, 3, TRIANGULAR);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*T)(i, j) + (*T)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == 0, true);

  *res = *tmp * *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - prod(*T, *T)) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  // Sym +,-,* Sym
  tmp = new MySimpleMatrix(*S);
  tmp2 = new MySimpleMatrix(*S);
  res = new MySimpleMatrix(3, 3, SYMMETRIC);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*S)(i, j) + (*S)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == 0, true);

  *res = *tmp * *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - prod(*S, *S)) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  // Sparse +,-,* Sparse
  tmp = new MySimpleMatrix(*SP);
  tmp2 = new MySimpleMatrix(*SP);
  res = new MySimpleMatrix(4, 4, SPARSE);
  *res = *tmp + *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*SP)(i, j) + (*SP)(i, j)), true);

  *res = *tmp - *tmp2;
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == 0, true);

  *res = *tmp * *tmp2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - prod(*SP, *SP)) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  // Banded +,- Banded
  tmp = new MySimpleMatrix(*B);
  tmp2 = new MySimpleMatrix(*B);
  res = new MySimpleMatrix(4, 4, BANDED);
  *res = *tmp + *tmp2;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == ((*B)(i, j) + (*B)(i, j)), true);

  *res = *tmp - *tmp2;
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == 0, true);

  delete tmp;
  delete tmp2;
  delete res;

  cout << "-->  test operators6 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators7()
{
  cout << "--> Test: operator7." << endl;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*D);
  tmp1->resize(4, 4);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*T);
  tmp2->resize(4, 4);
  MySiconosMatrix * tmp3 = new MySimpleMatrix(*S);
  tmp3->resize(4, 4);
  MySiconosMatrix * tmp4 = new MySimpleMatrix(*SP);
  MySiconosMatrix * tmp5 = new MySimpleMatrix(*B);

  MySiconosMatrix * res = new MySimpleMatrix(4, 4);

  // dense + ...
  // ... triang
  *res = add(*tmp1, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j), true);
  }
  // ... Sym
  *res = add(*tmp1, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = add(*tmp1, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp4)(i, j), true);
  // ... Banded
  *res = add(*tmp1, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*B)(i, j) , true);

  // dense - ...
  // ... triangular
  *res = sub(*tmp1, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) - (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j), true);
  }
  // ... Sym
  *res = sub(*tmp1, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) - (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) - (*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = sub(*tmp1, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) - (*tmp4)(i, j), true);
  // ... Banded
  *res = sub(*tmp1, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) - (*B)(i, j) , true);

  // triang + ...
  // ... dense
  *res = add(*tmp2, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j), true);
  }
  // ... Sym
  *res = add(*tmp2, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) + (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = add(*tmp2, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j), true);
  }

  // ... Banded
  *res = add(*tmp2, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) + (*B)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*B)(i, j) , true);

  // triang - ...
  // ... dense
  *res = sub(*tmp2, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) - (*tmp1)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp1)(i, j), true);
  }
  // ... Sym
  *res = sub(*tmp2, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) - (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = sub(*tmp2, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) - (*tmp4)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp4)(i, j), true);
  }

  // ... Banded
  *res = sub(*tmp2, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) - (*B)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*B)(i, j) , true);

  // sym + ...
  // ... dense
  *res = add(*tmp3, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp3)(j, i), true);
  }
  // ... triang
  *res = add(*tmp3, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) + (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = add(*tmp3, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp3)(j, i), true);
  }

  // ... Banded
  *res = add(*tmp3, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) + (*B)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i) + (*B)(i, j) , true);

  // sym - ...
  // ... dense
  *res = sub(*tmp3, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) - (*tmp1)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i) - (*tmp1)(i, j), true);
  }
  // ... triang
  *res = sub(*tmp3, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) - (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i), true);
  }
  // ... Sparse
  *res = sub(*tmp3, * tmp4);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) - (*tmp4)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i) - (*tmp4)(i, j), true);
  }

  // ... Banded
  *res = sub(*tmp3, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) - (*B)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i) - (*B)(i, j) , true);

  // sparse + ...
  // ... dense
  *res = add(*tmp4, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp4)(i, j), true);
  // ... triang
  *res = add(*tmp4, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j), true);
  }
  // ... Sym
  *res = add(*tmp4, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp3)(j, i), true);
  }
  // ... Banded
  *res = add(*tmp4, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*B)(i, j) , true);

  // sparse - ...
  // ... dense
  *res = sub(*tmp4, * tmp1);
  for (unsigned int i = 0; i < res->size1(); ++i)
    for (unsigned int j = 0 ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) - (*tmp1)(i, j), true);
  // ... triangular
  *res = sub(*tmp4, * tmp2);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) - (*tmp2)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j), true);
  }
  // ... Sym
  *res = sub(*tmp4, * tmp3);
  for (unsigned int i = 0; i < res->size1(); ++i)
  {
    for (unsigned int j = i ; j < res->size2(); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) - (*tmp3)(i, j), true);
    for (unsigned int j = 0 ; j < i; ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) - (*tmp3)(j, i), true);
  }
  // ... Banded
  *res = sub(*tmp4, * tmp5);
  for (signed i = 0; i < signed(B->size1()); ++ i)
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) - (*B)(i, j) , true);

  // Banded + ...
  // ... dense
  *res = add(*tmp5, * tmp1);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j) + (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp1)(i, j), true);
  }
  // ... triang
  *res = add(*tmp5, * tmp2);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j) + (*tmp5)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp2)(i, j), true);
  }
  // ...sym
  *res = add(*tmp5, * tmp3);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j), true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j) + (*tmp5)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i) + (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(i, j), true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp3)(j, i), true);
  }
  //... sparse
  *res = add(*tmp5, * tmp4);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j) + (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp4)(i, j), true);
  }
  // Banded - ...
  // ... dense

  *res = sub(*tmp5, * tmp1);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp1)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) - (*tmp1)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp1)(i, j), true);
  }

  // ... triang
  *res = sub(*tmp5, * tmp2);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp2)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) - (*tmp2)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp2)(i, j), true);
  }

  // ...sym
  *res = sub(*tmp5, * tmp3);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp3)(i, j), true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp3)(j, i), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      if (j >= i)
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) - (*tmp3)(i, j) , true);
      else
        CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == (*tmp5)(i, j) - (*tmp3)(j, i) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      if (j >= i) CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp3)(i, j), true);
      else  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp3)(j, i), true);
  }

  //... sparse
  *res = sub(*tmp5, * tmp4);
  for (signed i = 0; i < signed(B->size1()); ++ i)
  {
    for (signed j = 0; j < std::max(i - 1, 0); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp4)(i, j), true);
    for (signed j = std::max(i - 1, 0); j < std::min(i + 2, signed(B->size2())); ++ j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp4)(i, j) + (*tmp5)(i, j) , true);
    for (signed j = std::min(i + 2, signed(B->size2())); j < signed(B->size2()); ++j)
      CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", (*res)(i, j) == -(*tmp4)(i, j), true);
  }

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete res;

  cout << "-->  test operators7 ended with success." << endl;
}
void MySimpleMatrixTest::testOperators8()
{
  cout << "--> Test: operator8." << endl;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*D);
  tmp1->resize(4, 4);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*T);
  tmp2->resize(4, 4);
  MySiconosMatrix * tmp3 = new MySimpleMatrix(*S);
  tmp3->resize(4, 4);
  MySiconosMatrix * tmp4 = new MySimpleMatrix(*SP);
  MySiconosMatrix * tmp5 = new MySimpleMatrix(*B);

  MySiconosMatrix * res = new MySimpleMatrix(4, 4);

  // Dense * ...
  // triang
  *res = prod(*tmp1, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp2->getTriang())) == 0, true);
  // Sym
  *res = prod(*tmp1, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp3->getSym())) == 0, true);
  // Sparse
  *res = prod(*tmp1, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp4->getSparse())) == 0, true);
  // Banded
  *res = prod(*tmp1, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp5->getBanded())) == 0, true);
  // triang * ...
  // dense
  *res = prod(*tmp2, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp1->getDense())) == 0, true);
  // Sym
  *res = prod(*tmp2, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp3->getSym())) == 0, true);
  // Sparse
  *res = prod(*tmp2, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp4->getSparse())) == 0, true);
  // Banded
  *res = prod(*tmp2, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp5->getBanded())) == 0, true);
  // sym * ...
  // dense
  *res = prod(*tmp3, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp1->getDense())) == 0, true);
  // triang
  *res = prod(*tmp3, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp2->getTriang())) == 0, true);
  // Sparse
  *res = prod(*tmp3, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp4->getSparse())) == 0, true);
  // Banded
  *res = prod(*tmp3, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), tmp5->getBanded())) == 0, true);
  // Sparse * ...
  // dense
  *res = prod(*tmp4, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp1->getDense())) == 0, true);
  // triang
  *res = prod(*tmp4, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp2->getTriang())) == 0, true);
  // Sym
  *res = prod(*tmp4, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp3->getSym())) == 0, true);
  // Banded
  *res = prod(*tmp4, *tmp5);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp5->getBanded())) == 0, true);
  // Banded * ...
  // dense
  *res = prod(*tmp5, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp1->getDense())) == 0, true);
  // triang
  *res = prod(*tmp5, *tmp2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp2->getTriang())) == 0, true);
  // Sparse
  *res = prod(*tmp5, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp4->getSparse())) == 0, true);
  // Sym
  *res = prod(*tmp5, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp3->getSym())) == 0, true);

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  delete res;
  cout << "-->  test operators8 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators9()
{
  cout << "--> Test: operator9." << endl;
  double m = 2.2;
  int i = 3;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*D);
  MySiconosMatrix * res = new MySimpleMatrix(2, 2);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense()*m) == 0, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense()*i) == 0, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense()*m) == 0, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense()*i) == 0, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense() / m) == 0, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - tmp1->getDense() / i) == 0, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators9 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators10()
{
  cout << "--> Test: operator10." << endl;
  double m = 2.2;
  int i = 3;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*T);
  MySiconosMatrix * res = new MySimpleMatrix(3, 3, TRIANGULAR);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*m) == 0, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*i) == 0, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*m) == 0, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang()*i) == 0, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() / m) == 0, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getTriang() - tmp1->getTriang() / i) == 0, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators10 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators11()
{
  cout << "--> Test: operator11." << endl;
  double m = 2.2;
  int i = 3;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*S);
  MySiconosMatrix * res = new MySimpleMatrix(3, 3, SYMMETRIC);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*m) == 0, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*i) == 0, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*m) == 0, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym()*i) == 0, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym() / m) == 0, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSym() - tmp1->getSym() / i) == 0, true);
  delete tmp1;
  delete res;
  cout << "-->  test operator11 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators12()
{
  cout << "--> Test: operator12." << endl;
  double m = 2.2;
  int i = 3;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*SP);
  MySiconosMatrix * res = new MySimpleMatrix(4, 4, SPARSE);
  *res = m ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*m) == 0, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*i) == 0, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*m) == 0, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse()*i) == 0, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() / m) == 0, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getSparse() - tmp1->getSparse() / i) == 0, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators12 ended with success." << endl;
}

void MySimpleMatrixTest::testOperators13()
{
  cout << "--> Test: operator13." << endl;
  double m = 2.2;
  int i = 3;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*B);
  MySiconosMatrix * res = new MySimpleMatrix(4, 4, BANDED);
  *res = m * *tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded()*m) == 0, true);
  *res = i ** tmp1;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded()*i) == 0, true);
  *res = *tmp1 * m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded()*m) == 0, true);
  *res = *tmp1 * i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded()*i) == 0, true);
  *res = *tmp1 / m;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded() / m) == 0, true);
  *res = *tmp1 / i;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getBanded() - tmp1->getBanded() / i) == 0, true);
  delete tmp1;
  delete res;
  cout << "-->  test operators13 ended with success." << endl;
}


void MySimpleMatrixTest::testMultTranspose()
{
  cout << "--> Test: multTranspose." << endl;
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*D);
  tmp1->resize(4, 4);
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*T);
  tmp2->resize(4, 4);
  MySiconosMatrix * tmp3 = new MySimpleMatrix(*S);
  tmp3->resize(4, 4);
  MySiconosMatrix * tmp4 = new MySimpleMatrix(*SP);
  MySiconosMatrix * tmp5 = new MySimpleMatrix(*B);

  MySiconosMatrix * res = new MySimpleMatrix(4, 4);
  //Dense * ...
  *res = multTranspose(*tmp1, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), trans(tmp1->getDense()))) == 0, true);
  *res = multTranspose(*tmp1, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), tmp3->getSym())) == 0, true);
  *res = multTranspose(*tmp1, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp1->getDense(), trans(tmp4->getSparse()))) == 0, true);
  // Triang * ...
  *res = multTranspose(*tmp2, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), trans(tmp1->getDense()))) == 0, true);
  *res = multTranspose(*tmp2, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), tmp3->getSym())) == 0, true);
  *res = multTranspose(*tmp2, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp2->getTriang(), trans(tmp4->getSparse()))) == 0, true);
  // Sym * ...
  *res = multTranspose(*tmp3, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), trans(tmp1->getDense()))) == 0, true);
  *res = multTranspose(*tmp3, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp3 * *tmp3, true);
  *res = multTranspose(*tmp3, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp3->getSym(), trans(tmp4->getSparse()))) == 0, true);
  // Sparse* ...
  *res = multTranspose(*tmp4, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), trans(tmp1->getDense()))) == 0, true);
  *res = multTranspose(*tmp4, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), tmp3->getSym())) == 0, true);
  *res = multTranspose(*tmp4, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp4->getSparse(), trans(tmp4->getSparse()))) == 0, true);
  // Banded * ...
  *res = multTranspose(*tmp5, *tmp1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), trans(tmp1->getDense()))) == 0, true);
  *res = multTranspose(*tmp5, *tmp3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), tmp3->getSym())) == 0, true);
  *res = multTranspose(*tmp5, *tmp4);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", norm_inf(res->getDense() - prod(tmp5->getBanded(), trans(tmp4->getSparse()))) == 0, true);

  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  cout << "-->  test multTranspose ended with success." << endl;
}
void MySimpleMatrixTest::testPow()
{
  cout << "--> Test: pow." << endl;
  // Dense
  MySiconosMatrix * tmp1 = new MySimpleMatrix(*D);
  MySiconosMatrix * res = new MySimpleMatrix(2, 2);
  *res = pow(*tmp1, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp1**tmp1**tmp1, true);
  delete res;
  // Triang
  MySiconosMatrix * tmp2 = new MySimpleMatrix(*T);
  res = new MySimpleMatrix(3, 3, TRIANGULAR);
  *res = pow(*tmp2, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp2**tmp2**tmp2, true);
  delete res;
  // Sym
  MySiconosMatrix * tmp3 = new MySimpleMatrix(*S);
  res = new MySimpleMatrix(3, 3, SYMMETRIC);
  *res = pow(*tmp3, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp3**tmp3**tmp3, true);
  delete res;
  // Sparse
  MySiconosMatrix * tmp4 = new MySimpleMatrix(*SP);
  res = new MySimpleMatrix(4, 4, SPARSE);
  *res = pow(*tmp4, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == *tmp4**tmp4**tmp4, true);
  delete res;
  // Banded
  MySiconosMatrix * tmp5 = new MySimpleMatrix(*B);
  res = new MySimpleMatrix(4, 4);
  *res = pow(*tmp5, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators: ", *res == prod(*tmp5, *tmp5) == 0, true);
  delete res;
  delete tmp1;
  delete tmp2;
  delete tmp3;
  delete tmp4;
  delete tmp5;
  cout << "-->  test pow ended with success." << endl;
}

void MySimpleMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of MySimpleMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
