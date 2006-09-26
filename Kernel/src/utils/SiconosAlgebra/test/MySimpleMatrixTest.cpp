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
  SimM->display();
  SicM->display();
  std::vector<double> v3(2, 0);
  std::vector<double> v4(2, 0);
  std::vector<double> v5(3, 0);
  v4[0] = 6;
  v4[1] = 9;
  v5[0] = 8;
  v5[1] = 9;
  v5[2] = 10;

  vect1 = new MySimpleVector(DENSE, v3, 2);
  vect2 = new MySimpleVector(DENSE, v4, 2); // vect2 != vect1, but vect2 == SimM second column
  vect3 = new MySimpleVector(DENSE, v5, 3); // vect3 != vect1, but vect3 == SimM second row

  //   unsigned int i,j;
  //   unsigned int Arow = 10;
  //   unsigned int Acol = 10;
  //   unsigned int Brow = 10;
  //   unsigned int Bcol = 5;
  //   unsigned int Vsize = 10;

  //   A = new MySimpleMatrix(Arow, Acol);
  //   B = new MySimpleMatrix(Brow, Bcol);
  //   C = new MySimpleMatrix(5, 5);


  //   /* init A, C */
  //   srand( (unsigned)time( NULL ) );
  //   for (i=0; i < Arow; i++)
  //     for (unsigned int j=0; j < Acol; j++)
  //       {
  //  (*A)(i,j) = rand()%10+20;
  //       }

  //   /* init B */
  //   for (i=0; i < Brow; i++)
  //     for (j=0; j < Bcol; j++)
  //       (*B)(i,j) = rand()%100-50;
  //   for (i=0; i < 5; i++)
  //     for (j=0; j < 5; j++)
  //       (*C)(i,j) = rand()%100-50;


  //   /* init SV */
  //   vector<double> vtmp(Vsize);
  //   for (i=0; i < Vsize; i++)
  //     vtmp.at(i) = rand()%10+20;

  //   SV = new SimpleVector(vtmp);
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
  MySimpleMatrix * test = new MySimpleMatrix(DENSE, 2, 3);

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
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", *test == *SicM, true);
  cout << "--> Constructor (copy) test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor3()
{
  cout << "--> Test: constructor 3." << endl;
  MySimpleMatrix * test = new MySimpleMatrix(SPARSE);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor0 : ", test->getNum() == 4, true);
  delete test;
  cout << "--> Constructor 3 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  cout << "--> Constructor 4 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  cout << "--> Constructor 5 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor6()
{
  cout << "--> Test: constructor 6." << endl;
  cout << "--> Constructor 6 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor7()
{
  cout << "--> Test: constructor 7." << endl;
  cout << "--> Constructor 7 test ended with success." << endl;
}
void MySimpleMatrixTest::testConstructor8()
{
  cout << "--> Test: constructor 8." << endl;
  cout << "--> Constructor 8 test ended with success." << endl;
}

void MySimpleMatrixTest::testConstructor9()
{
  cout << "--> Test: constructor 9." << endl;
  const std::vector<double> v1(4, 1);
  const std::vector<double> v(6, 1);
  //DENSE
  MySiconosMatrix *test1 = new MySimpleMatrix(DENSE, v1, 2, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->getNum() == 1, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test1->size2() == 2, true);
  //TRIANGULAR
  MySiconosMatrix *test2 = new MySimpleMatrix(TRIANGULAR, v1, 2, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->getNum() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test2->size2() == 2, true);
  //SYMMETRIC
  MySiconosMatrix *test3 = new MySimpleMatrix(SYMMETRIC, v1, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->getNum() == 3, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->size1() == 2, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor9 : ", test3->size2() == 2, true);
  //BANDED
  MySiconosMatrix *test4 = new MySimpleMatrix(BANDED, v, 2, 2, 1, 1);
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
  MySimpleVector * tmp = new MySimpleVector(DENSE, 3);
  SimM->getRow(1, *tmp);
  SimM->display();
  tmp->display();
  vect3->display();
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetSetRowCol : ", *tmp == *vect3, true);

  MySimpleVector * tmp2 = new MySimpleVector(DENSE, 2);
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
  MySiconosMatrix * block = new MySimpleMatrix(DENSE, 2, 2);
  full->getBlock(1, 1, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(0, 0) == 9 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(0, 1) == 10 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(1, 0) == 0 , true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetBlock: ", (*block)(1, 1) == 0 , true);
  delete block;
  delete full;
  cout << "--> getBlock test ended with success." << endl;
}

void MySimpleMatrixTest::testBlockMatrixCopy()
{
  cout << "--> Test: blockMatrixCopy." << endl;
  // Dense
  MySiconosMatrix * full = new MySimpleMatrix(*SimM);
  full->resize(4, 5);
  full->display();
  full->blockMatrixCopy(*SicM, 2, 3);
  full->display();
  MySiconosMatrix * block = new MySimpleMatrix(DENSE, 2, 2);
  full->getBlock(2, 3, *block);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy: ", *block == *SicM , true);
  delete block;
  delete full;
  // Triang
  // Sym
  // Sparse
  // Banded
  cout << "-->  blockMatrixCopy test ended with success." << endl;
}

void MySimpleMatrixTest::testOperators()
{
  cout << "--> Test: operators." << endl;
  //   MySimpleMatrix test(*A);
  //   MySimpleMatrix test2(10,10);
  //   test+=*A;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators: B.blockMatrixCopy( A, x, y) ", test == (2**A),true);
  //   test-=*A;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators", test == (*A),true);
  //   test*=3;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators", test == (3**A),true);
  //   MySimpleMatrix X(2,3);
  //   X(0,0) = 1;  X(0,1) = 2;  X(0,2) = 3;
  //   X(1,0) = 0;  X(1,1) = 1;  X(1,2) = 1;
  //   MySimpleMatrix Y(3,3);
  //   Y(0,0) = 2;  Y(0,1) = 3;  Y(0,2) = 1;
  //   Y(1,0) = 0;  Y(1,1) = 1;  Y(1,2) = 0;
  //   Y(2,0) = 1;  Y(2,1) = 0;  Y(2,2) = 2;
  //   MySimpleMatrix Z(2,3);
  //   Z(0,0) = 5;  Z(0,1) = 5;  Z(0,2) = 7;
  //   Z(1,0) = 1;  Z(1,1) = 1;  Z(1,2) = 2;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", X*Y == Z,true);
  //   test2 = 2 * test;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 6**A,true);
  //   test2 = test * 2;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 6**A,true);
  //   test2 = test / 2;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 1.5**A,true);
  //   test2 = test + *A;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 4**A,true);
  //   test2 = test - *A;
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 2**A,true);
  //   test2 = pow(*A,3);
  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == ((*A)*(*A)*(*A)),true);
  cout << "-->  test operators ended with success." << endl;
}

void MySimpleMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  cout << "-->  test assignment ended with success." << endl;
}

void MySimpleMatrixTest::testMultTranspose()
{
  cout << "--> Test: multTranspose." << endl;
  //   MySimpleMatrix X(2,3);
  //   X(0,0) = 1;  X(0,1) = 2;  X(0,2) = 3;
  //   X(1,0) = 0;  X(1,1) = 1;  X(1,2) = 1;
  //   MySimpleMatrix Y(3,3);
  //   Y(0,0) = 2;  Y(0,1) = 0;  Y(0,2) = 1;
  //   Y(1,0) = 3;  Y(1,1) = 1;  Y(1,2) = 0;
  //   Y(2,0) = 1;  Y(2,1) = 0;  Y(2,2) = 2;
  //   MySimpleMatrix Z(2,3);
  //   Z(0,0) = 5;  Z(0,1) = 5;  Z(0,2) = 7;
  //   Z(1,0) = 1;  Z(1,1) = 1;  Z(1,2) = 2;

  //   CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", X.multTranspose(Y)==Z, true);
  cout << "-->  test multTranspose ended with success." << endl;
}

// void MySimpleMatrixTest::testBlockMatrixCopy1()
// {
//   cout << "--> Test: BlockMatrixCopy1." << endl;
//   MySimpleMatrix Aa(1,2), Bb(3,3), Cc(3,3);
//   Aa(0,0) = 1.0; Aa(0,1) = 2.0;
//   Bb(0,0) = 0.0; Bb(0,1) = 0.0; Bb(0,2) = 0.0;
//   Bb(1,0) = 0.0; Bb(1,1) = 0.0; Bb(1,2) = 0.0;
//   Bb(2,0) = 0.0; Bb(2,1) = 0.0; Bb(2,2) = 0.0;
//   Cc(0,0) = 0.0; Cc(0,1) = 0.0; Cc(0,2) = 0.0;
//   Cc(1,0) = 1.0; Cc(1,1) = 2.0; Cc(1,2) = 0.0;
//   Cc(2,0) = 0.0; Cc(2,1) = 0.0; Cc(2,2) = 0.0;

//   Bb.blockMatrixCopy( Aa, 1, 0);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy1 : B.blockMatrixCopy( A, x, y) ", Bb == Cc, true);
//   cout << "-->  test BlockMatrixCopy1 ended with success." << endl;
// }

// void MySimpleMatrixTest::testBlockMatrixCopy2()
// {
//   cout << "--> Test: BlockMatrixCopy2." << endl;
//   MySimpleMatrix Aa(1,1), Bb(1,1), Cc(1,1);
//   Aa(0,0) = 1.0;
//   Bb(0,0) = 0.0;
//   Cc(0,0) = 1.0;

//   Bb.blockMatrixCopy( Aa, 0, 0);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy2 : B.blockMatrixCopy( A, x, y) ", Bb == Cc,true);
//   cout << "-->  test BlockMatrixCopy2 ended with success." << endl;
// }

// void MySimpleMatrixTest::testBlockMatrixCopy3()
// {
//   cout << "--> Test: BlockMatrixCopy3." << endl;
//   MySiconosMatrix * Aa = new MySimpleMatrix(2,2);
//   MySiconosMatrix * Bb = new MySimpleMatrix(4,4);
//   MySiconosMatrix * Cc = new MySimpleMatrix(4,4);

//   (*Aa)(0,0) = 1.0; (*Aa)(0,1) = 2.0;
//   (*Aa)(1,0) = 3.0; (*Aa)(1,1) = 4.0;
//   (*Bb)(0,0) = 0.0; (*Bb)(0,1) = 0.0; (*Bb)(0,2) = 0.0; (*Bb)(0,3) = 0.0;
//   (*Bb)(1,0) = 0.0; (*Bb)(1,1) = 0.0; (*Bb)(1,2) = 0.0; (*Bb)(1,3) = 0.0;
//   (*Bb)(2,0) = 0.0; (*Bb)(2,1) = 0.0; (*Bb)(2,2) = 0.0; (*Bb)(2,3) = 0.0;
//   (*Bb)(3,0) = 0.0; (*Bb)(3,1) = 0.0; (*Bb)(3,2) = 0.0; (*Bb)(3,3) = 0.0;
//   (*Cc)(0,0) = 1.0; (*Cc)(0,1) = 2.0; (*Cc)(0,2) = 0.0; (*Cc)(0,3) = 0.0;
//   (*Cc)(1,0) = 3.0; (*Cc)(1,1) = 4.0; (*Cc)(1,2) = 0.0; (*Cc)(1,3) = 0.0;
//   (*Cc)(2,0) = 0.0; (*Cc)(2,1) = 0.0; (*Cc)(2,2) = 0.0; (*Cc)(2,3) = 0.0;
//   (*Cc)(3,0) = 0.0; (*Cc)(3,1) = 0.0; (*Cc)(3,2) = 0.0; (*Cc)(3,3) = 0.0;

//   (*Bb).blockMatrixCopy( (*Aa), 0, 0);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy3 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc),true);
//   delete Aa;
//   delete Bb;
//   delete Cc;
//   cout << "-->  test BlockMatrixCopy3 ended with success." << endl;
// }

// void MySimpleMatrixTest::testBlockMatrixCopy4()
// {
//   cout << "--> Test: BlockMatrixCopy4." << endl;
//   MySiconosMatrix *Aa = new MySimpleMatrix(2,2), *Bb = new MySimpleMatrix(4,4), *Cc = new MySimpleMatrix(4,4);
//   (*Aa)(0,0) = 1.0; (*Aa)(0,1) = 2.0;
//   (*Aa)(1,0) = 3.0; (*Aa)(1,1) = 4.0;
//   (*Bb)(0,0) = 0.0; (*Bb)(0,1) = 0.0; (*Bb)(0,2) = 0.0; (*Bb)(0,3) = 0.0;
//   (*Bb)(1,0) = 0.0; (*Bb)(1,1) = 0.0; (*Bb)(1,2) = 0.0; (*Bb)(1,3) = 0.0;
//   (*Bb)(2,0) = 0.0; (*Bb)(2,1) = 0.0; (*Bb)(2,2) = 0.0; (*Bb)(2,3) = 0.0;
//   (*Bb)(3,0) = 0.0; (*Bb)(3,1) = 0.0; (*Bb)(3,2) = 0.0; (*Bb)(3,3) = 0.0;
//   (*Cc)(0,0) = 0.0; (*Cc)(0,1) = 0.0; (*Cc)(0,2) = 0.0; (*Cc)(0,3) = 0.0;
//   (*Cc)(1,0) = 0.0; (*Cc)(1,1) = 0.0; (*Cc)(1,2) = 0.0; (*Cc)(1,3) = 0.0;
//   (*Cc)(2,0) = 0.0; (*Cc)(2,1) = 0.0; (*Cc)(2,2) = 1.0; (*Cc)(2,3) = 2.0;
//   (*Cc)(3,0) = 0.0; (*Cc)(3,1) = 0.0; (*Cc)(3,2) = 3.0; (*Cc)(3,3) = 4.0;

//   (*Bb).blockMatrixCopy( (*Aa), 2, 2);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc),true);
//   delete Aa;
//   delete Bb;
//   delete Cc;
//   cout << "-->  test BlockMatrixCopy4 ended with success." << endl;
// }

// void MySimpleMatrixTest::testBlockMatrixCopyException1()
// {
//   MySiconosMatrix *Aa = new MySimpleMatrix(2,2), *Bb = new MySimpleMatrix(4,4);
//   (*Bb).blockMatrixCopy( (*Aa), 3, 3);
//   delete Aa;
//   delete Bb;
// }

// void MySimpleMatrixTest::testBlockMatrixCopyException2()
// {
//   MySiconosMatrix *Aa = new MySimpleMatrix(6,6), *Bb = new MySimpleMatrix(4,4);
//   (*Bb).blockMatrixCopy( (*Aa), 0, 0);
//   delete Aa;
//   delete Bb;
// }


// void MySimpleMatrixTest::testLinearSolve()
// {
//   cout << "--> Test: linearSolve." << endl;
//   MySimpleMatrix X(10,5);
//   A->linearSolve((*B), X);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve : A==A", (*A)==(*A), true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve : B==B", (*B)==(*B), true);
//   CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve 3: ", ((*A)*X)==(*B), true);
//   cout << "--> linearSolve test ended with success." << endl;
// }

// void MySimpleMatrixTest::testSizeException()
// {
//   *A+(*B);
// }

// void MySimpleMatrixTest::testConstructorException()
// {
//   LaVectorDouble LVD(10);
//   for (int i=0; i < LVD.size(); i++)
//     LVD(i) = i;
//   MySimpleMatrix X(LVD,4,2);
// }

void MySimpleMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of MySimpleMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
