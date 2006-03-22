/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#include "SimpleMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SimpleMatrixTest);


void SimpleMatrixTest::setUp()
{
  unsigned int i, j;
  unsigned int Arow = 10;
  unsigned int Acol = 10;
  unsigned int Brow = 10;
  unsigned int Bcol = 5;
  unsigned int Vsize = 10;

  A = new SimpleMatrix(Arow, Acol);
  B = new SimpleMatrix(Brow, Bcol);
  C = new SimpleMatrix(5, 5);


  /* init A, C */
  srand((unsigned)time(NULL));
  for (i = 0; i < Arow; i++)
    for (unsigned int j = 0; j < Acol; j++)
    {
      (*A)(i, j) = rand() % 10 + 20;
    }

  /* init B */
  for (i = 0; i < Brow; i++)
    for (j = 0; j < Bcol; j++)
      (*B)(i, j) = rand() % 100 - 50;
  for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++)
      (*C)(i, j) = rand() % 100 - 50;


  /* init SV */
  vector<double> vtmp(Vsize);
  for (i = 0; i < Vsize; i++)
    vtmp.at(i) = rand() % 10 + 20;

  SV = new SimpleVector(vtmp);
}

void SimpleMatrixTest::tearDown()
{
  delete A;
  delete B;
  delete C;
  delete SV;
}

//______________________________________________________________________________

void SimpleMatrixTest::testConstructor1()
{
  cout << "====================================" << endl;
  cout << "=== Simple Matrix tests start ...=== " << endl;
  cout << "====================================" << endl;
  cout << "--> Test: constructor 1." << endl;
  SimpleMatrix test(*A);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : ", test == *A, true);
  cout << "--> Constructor 1 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor2()
{
  cout << "--> Test: constructor 2." << endl;
  unsigned int row = 120;
  unsigned int col = 500;
  SimpleMatrix X(row, col);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : X.size(0) == row", X.size(0) == row, true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : X.size(1) == col", X.size(1) == col, true);
  cout << "--> Constructor 2 test ended with success." << endl;
}
void SimpleMatrixTest::testConstructor3()
{
  cout << "--> Test: constructor 3." << endl;
  LaGenMatDouble tmp;
  tmp = A->getLaGenMatDouble();
  SimpleMatrix test(tmp);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor2 : ", test == *A, true);
  cout << "--> Constructor 3 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor4()
{
  cout << "--> Test: constructor 4." << endl;
  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;

  SimpleMatrix X(LVD, 5, 2);

  SimpleMatrix Xtrue(5, 2);
  for (unsigned int i = 0; i < 5; i++)
    for (unsigned int j = 0; j < 2; j++)
      Xtrue(i, j) = 2 * i + j;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor4 : X == Xtrue", X == Xtrue, true);
  cout << "--> Constructor 4 test ended with success." << endl;
}

void SimpleMatrixTest::testConstructor5()
{
  cout << "--> Test: constructor 5." << endl;
  A->write("tmp.dat", "ascii");
  SimpleMatrix test("tmp.dat", 1);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", test == *A, true);
  B->write("tmp2.dat", "binary");
  SimpleMatrix test2("tmp2.dat", 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor5 : ", test2 == *B, true);
  cout << "--> Constructor 5 test ended with success." << endl;
}

void SimpleMatrixTest::testGetRow()
{
  cout << "--> Test: GetRow." << endl;
  SimpleMatrix *X = new SimpleMatrix(*A);
  SimpleVector *V2 = new SimpleVector(10);
  X->getRow(2, *V2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetRow : X == A ", *X == (*A), true);
  X->getRow(3, *SV);
  SimpleVector * V3 = new SimpleVector(10);
  X->getRow(3, *V3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetRow : V2 != V3 ", *V2 == *V3, false);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testGetRow : V3 == SV ", *V3 == *SV, true);
  delete V2;
  delete V3;
  delete X;
  cout << "--> GetRow test ended with success." << endl;
}

void SimpleMatrixTest::testReadWriteBinary()
{
  cout << "--> Test: readWriteBinary." << endl;
  SimpleMatrix X(*B);
  (*A).write("totoMat.bin", "binary");
  X.read("totoMat.bin", "binary");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWriteBinary : A == X ", (*A) == X, true);
  cout << "--> readWriteBinary test ended with success." << endl;
}

void SimpleMatrixTest::testReadWriteAscii()
{
  cout << "--> Test: readWriteAscii." << endl;
  SimpleMatrix X(*B);
  (*A).write("totoMat.ascii", "ascii");
  X.read("totoMat.ascii", "ascii");
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testReadWriteAscii : A == X ", (*A) == X, true);
  cout << "-->  test readWriteAscii ended with success." << endl;
}

void SimpleMatrixTest::testAssignment()
{
  cout << "--> Test: assignment." << endl;
  SimpleMatrix X(10, 10);
  X = (*A);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", X == *A, true);
  cout << "-->  test assignment ended with success." << endl;
}

void SimpleMatrixTest::testMultTranspose()
{
  cout << "--> Test: multTranspose." << endl;
  SimpleMatrix X(2, 3);
  X(0, 0) = 1;
  X(0, 1) = 2;
  X(0, 2) = 3;
  X(1, 0) = 0;
  X(1, 1) = 1;
  X(1, 2) = 1;
  SimpleMatrix Y(3, 3);
  Y(0, 0) = 2;
  Y(0, 1) = 0;
  Y(0, 2) = 1;
  Y(1, 0) = 3;
  Y(1, 1) = 1;
  Y(1, 2) = 0;
  Y(2, 0) = 1;
  Y(2, 1) = 0;
  Y(2, 2) = 2;
  SimpleMatrix Z(2, 3);
  Z(0, 0) = 5;
  Z(0, 1) = 5;
  Z(0, 2) = 7;
  Z(1, 0) = 1;
  Z(1, 1) = 1;
  Z(1, 2) = 2;

  CPPUNIT_ASSERT_EQUAL_MESSAGE("testAssignment :", X.multTranspose(Y) == Z, true);
  cout << "-->  test multTranspose ended with success." << endl;
}

void SimpleMatrixTest::testBlockMatrixCopy1()
{
  cout << "--> Test: BlockMatrixCopy1." << endl;
  SimpleMatrix Aa(1, 2), Bb(3, 3), Cc(3, 3);
  Aa(0, 0) = 1.0;
  Aa(0, 1) = 2.0;
  Bb(0, 0) = 0.0;
  Bb(0, 1) = 0.0;
  Bb(0, 2) = 0.0;
  Bb(1, 0) = 0.0;
  Bb(1, 1) = 0.0;
  Bb(1, 2) = 0.0;
  Bb(2, 0) = 0.0;
  Bb(2, 1) = 0.0;
  Bb(2, 2) = 0.0;
  Cc(0, 0) = 0.0;
  Cc(0, 1) = 0.0;
  Cc(0, 2) = 0.0;
  Cc(1, 0) = 1.0;
  Cc(1, 1) = 2.0;
  Cc(1, 2) = 0.0;
  Cc(2, 0) = 0.0;
  Cc(2, 1) = 0.0;
  Cc(2, 2) = 0.0;

  Bb.blockMatrixCopy(Aa, 1, 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy1 : B.blockMatrixCopy( A, x, y) ", Bb == Cc, true);
  cout << "-->  test BlockMatrixCopy1 ended with success." << endl;
}

void SimpleMatrixTest::testBlockMatrixCopy2()
{
  cout << "--> Test: BlockMatrixCopy2." << endl;
  SimpleMatrix Aa(1, 1), Bb(1, 1), Cc(1, 1);
  Aa(0, 0) = 1.0;
  Bb(0, 0) = 0.0;
  Cc(0, 0) = 1.0;

  Bb.blockMatrixCopy(Aa, 0, 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy2 : B.blockMatrixCopy( A, x, y) ", Bb == Cc, true);
  cout << "-->  test BlockMatrixCopy2 ended with success." << endl;
}

void SimpleMatrixTest::testBlockMatrixCopy3()
{
  cout << "--> Test: BlockMatrixCopy3." << endl;
  SiconosMatrix * Aa = new SimpleMatrix(2, 2);
  SiconosMatrix * Bb = new SimpleMatrix(4, 4);
  SiconosMatrix * Cc = new SimpleMatrix(4, 4);

  (*Aa)(0, 0) = 1.0;
  (*Aa)(0, 1) = 2.0;
  (*Aa)(1, 0) = 3.0;
  (*Aa)(1, 1) = 4.0;
  (*Bb)(0, 0) = 0.0;
  (*Bb)(0, 1) = 0.0;
  (*Bb)(0, 2) = 0.0;
  (*Bb)(0, 3) = 0.0;
  (*Bb)(1, 0) = 0.0;
  (*Bb)(1, 1) = 0.0;
  (*Bb)(1, 2) = 0.0;
  (*Bb)(1, 3) = 0.0;
  (*Bb)(2, 0) = 0.0;
  (*Bb)(2, 1) = 0.0;
  (*Bb)(2, 2) = 0.0;
  (*Bb)(2, 3) = 0.0;
  (*Bb)(3, 0) = 0.0;
  (*Bb)(3, 1) = 0.0;
  (*Bb)(3, 2) = 0.0;
  (*Bb)(3, 3) = 0.0;
  (*Cc)(0, 0) = 1.0;
  (*Cc)(0, 1) = 2.0;
  (*Cc)(0, 2) = 0.0;
  (*Cc)(0, 3) = 0.0;
  (*Cc)(1, 0) = 3.0;
  (*Cc)(1, 1) = 4.0;
  (*Cc)(1, 2) = 0.0;
  (*Cc)(1, 3) = 0.0;
  (*Cc)(2, 0) = 0.0;
  (*Cc)(2, 1) = 0.0;
  (*Cc)(2, 2) = 0.0;
  (*Cc)(2, 3) = 0.0;
  (*Cc)(3, 0) = 0.0;
  (*Cc)(3, 1) = 0.0;
  (*Cc)(3, 2) = 0.0;
  (*Cc)(3, 3) = 0.0;

  (*Bb).blockMatrixCopy((*Aa), 0, 0);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy3 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc), true);
  delete Aa;
  delete Bb;
  delete Cc;
  cout << "-->  test BlockMatrixCopy3 ended with success." << endl;
}

void SimpleMatrixTest::testBlockMatrixCopy4()
{
  cout << "--> Test: BlockMatrixCopy4." << endl;
  SiconosMatrix *Aa = new SimpleMatrix(2, 2), *Bb = new SimpleMatrix(4, 4), *Cc = new SimpleMatrix(4, 4);
  (*Aa)(0, 0) = 1.0;
  (*Aa)(0, 1) = 2.0;
  (*Aa)(1, 0) = 3.0;
  (*Aa)(1, 1) = 4.0;
  (*Bb)(0, 0) = 0.0;
  (*Bb)(0, 1) = 0.0;
  (*Bb)(0, 2) = 0.0;
  (*Bb)(0, 3) = 0.0;
  (*Bb)(1, 0) = 0.0;
  (*Bb)(1, 1) = 0.0;
  (*Bb)(1, 2) = 0.0;
  (*Bb)(1, 3) = 0.0;
  (*Bb)(2, 0) = 0.0;
  (*Bb)(2, 1) = 0.0;
  (*Bb)(2, 2) = 0.0;
  (*Bb)(2, 3) = 0.0;
  (*Bb)(3, 0) = 0.0;
  (*Bb)(3, 1) = 0.0;
  (*Bb)(3, 2) = 0.0;
  (*Bb)(3, 3) = 0.0;
  (*Cc)(0, 0) = 0.0;
  (*Cc)(0, 1) = 0.0;
  (*Cc)(0, 2) = 0.0;
  (*Cc)(0, 3) = 0.0;
  (*Cc)(1, 0) = 0.0;
  (*Cc)(1, 1) = 0.0;
  (*Cc)(1, 2) = 0.0;
  (*Cc)(1, 3) = 0.0;
  (*Cc)(2, 0) = 0.0;
  (*Cc)(2, 1) = 0.0;
  (*Cc)(2, 2) = 1.0;
  (*Cc)(2, 3) = 2.0;
  (*Cc)(3, 0) = 0.0;
  (*Cc)(3, 1) = 0.0;
  (*Cc)(3, 2) = 3.0;
  (*Cc)(3, 3) = 4.0;

  (*Bb).blockMatrixCopy((*Aa), 2, 2);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc), true);
  delete Aa;
  delete Bb;
  delete Cc;
  cout << "-->  test BlockMatrixCopy4 ended with success." << endl;
}

void SimpleMatrixTest::testBlockMatrixCopyException1()
{
  SiconosMatrix *Aa = new SimpleMatrix(2, 2), *Bb = new SimpleMatrix(4, 4);
  (*Bb).blockMatrixCopy((*Aa), 3, 3);
  delete Aa;
  delete Bb;
}

void SimpleMatrixTest::testBlockMatrixCopyException2()
{
  SiconosMatrix *Aa = new SimpleMatrix(6, 6), *Bb = new SimpleMatrix(4, 4);
  (*Bb).blockMatrixCopy((*Aa), 0, 0);
  delete Aa;
  delete Bb;
}
void SimpleMatrixTest::testOperators()
{
  cout << "--> Test: operators." << endl;
  SimpleMatrix test(*A);
  SimpleMatrix test2(10, 10);
  test += *A;
  CPPUNIT_ASSERT_EQUAL_MESSAGE(" testOperators: B.blockMatrixCopy( A, x, y) ", test == (2 * *A), true);
  test -= *A;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators", test == (*A), true);
  test *= 3;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testOperators", test == (3 * *A), true);
  SimpleMatrix X(2, 3);
  X(0, 0) = 1;
  X(0, 1) = 2;
  X(0, 2) = 3;
  X(1, 0) = 0;
  X(1, 1) = 1;
  X(1, 2) = 1;
  SimpleMatrix Y(3, 3);
  Y(0, 0) = 2;
  Y(0, 1) = 3;
  Y(0, 2) = 1;
  Y(1, 0) = 0;
  Y(1, 1) = 1;
  Y(1, 2) = 0;
  Y(2, 0) = 1;
  Y(2, 1) = 0;
  Y(2, 2) = 2;
  SimpleMatrix Z(2, 3);
  Z(0, 0) = 5;
  Z(0, 1) = 5;
  Z(0, 2) = 7;
  Z(1, 0) = 1;
  Z(1, 1) = 1;
  Z(1, 2) = 2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", X * Y == Z, true);
  test2 = 2 * test;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 6 * *A, true);
  test2 = test * 2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 6 * *A, true);
  test2 = test / 2;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 1.5 * *A, true);
  test2 = test + *A;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 4 * *A, true);
  test2 = test - *A;
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == 2 * *A, true);
  test2 = pow(*A, 3);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", test2 == ((*A) * (*A) * (*A)), true);
  cout << "-->  test operators ended with success." << endl;
}

void SimpleMatrixTest::testLinearSolve()
{
  cout << "--> Test: linearSolve." << endl;
  SimpleMatrix X(10, 5);
  A->linearSolve((*B), X);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve : A==A", (*A) == (*A), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve : B==B", (*B) == (*B), true);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testLinearSolve 3: ", ((*A)*X) == (*B), true);
  cout << "--> linearSolve test ended with success." << endl;
}

void SimpleMatrixTest::testSizeException()
{
  *A + (*B);
}

void SimpleMatrixTest::testConstructorException()
{
  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;
  SimpleMatrix X(LVD, 4, 2);
}

void SimpleMatrixTest::End()
{
  cout << "======================================" << endl;
  cout << " ===== End of SimpleMatrix Tests ===== " << endl;
  cout << "======================================" << endl;
}
