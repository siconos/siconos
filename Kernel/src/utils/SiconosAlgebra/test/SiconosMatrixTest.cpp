/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
#include "SiconosMatrixTest.h"

#define CPPUNIT_ASSERT_NOT_EQUAL(message, alpha, omega)      \
            if ((alpha) == (omega)) CPPUNIT_FAIL(message);

// on place cette classe de test dans le registry
CPPUNIT_TEST_SUITE_REGISTRATION(SiconosMatrixTest);


void SiconosMatrixTest::setUp()
{
  int i, j;
  int Arow = 10;
  int Acol = 10;
  int Brow = 10;
  int Bcol = 5;
  int Vsize = 10;

  A = new SiconosMatrix(Arow, Acol);
  B = new SiconosMatrix(Brow, Bcol);
  /*C==A*/
  C = new SiconosMatrix(Arow, Acol);


  /* init A, C */
  srand((unsigned)time(NULL));
  for (i = 0; i < Arow; i++)
    for (int j = 0; j < Acol; j++)
    {
      (*A)(i, j) = rand() % 10 + 20;
      (*C)(i, j) = (*A)(i, j);
    }

  /* init B */
  for (i = 0; i < Brow; i++)
    for (j = 0; j < Bcol; j++)
      (*B)(i, j) = rand() % 100 - 50;


  /* init SV */
  vector<double> vtmp(Vsize);
  for (i = 0; i < Vsize; i++)
    vtmp.at(i) = rand() % 10 + 20;

  SV = new SimpleVector(vtmp);
}

void SiconosMatrixTest::tearDown()
{
  delete A;
  delete B;
  delete C;
  delete SV;
}

//______________________________________________________________________________

void SiconosMatrixTest::testConstructor1()
{
  unsigned int row = 120;
  unsigned int col = 500;
  SiconosMatrix X(row, col);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : X.size(0) == row", X.size(0), row);
  CPPUNIT_ASSERT_EQUAL_MESSAGE("testConstructor1 : X.size(1) == col", X.size(1), col);
  cout << " SiconosMatrixTest >>> testConstructor1 .................................................. OK\n ";
}

void SiconosMatrixTest::testConstructor2()
{

  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;

  SiconosMatrix X(LVD, 5, 2);

  SiconosMatrix Xtrue(5, 2);
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 2; j++)
      Xtrue(i, j) = 2 * i + j;

  CPPUNIT_ASSERT_MESSAGE("testConstructor2 : X == Xtrue", X == Xtrue);
  cout << "SiconosMatrixTest >>> testConstructor2 .................................................. OK\n ";
}

void SiconosMatrixTest::copyConstructor()
{
  SiconosMatrix *X = new SiconosMatrix(*A);
  cout <<  "DSKDSLDK" << endl;
  CPPUNIT_ASSERT_MESSAGE("copyConstructor2 : X == A", *X == (*A));
  cout << "SiconosMatrixTest >>> copyConstructor .................................................. OK\n ";
  delete X;
}

void SiconosMatrixTest::testEquality()
{
  CPPUNIT_ASSERT_MESSAGE("testEquality : A==A", (*A) == (*A));
  CPPUNIT_ASSERT_MESSAGE("testEquality : B==B", (*B) == (*B));
  CPPUNIT_ASSERT_MESSAGE("testEquality : C==C", (*C) == (*C));

  CPPUNIT_ASSERT_MESSAGE("testEquality : A!=B", (*A) != (*B));
  CPPUNIT_ASSERT_MESSAGE("testEquality : B!=C", (*B) != (*C));
  CPPUNIT_ASSERT_MESSAGE("testEquality : A==C", (*A) == (*C));
  (*C)(5, 5) = (*C)(5, 5) * 2;
  CPPUNIT_ASSERT_MESSAGE("testEquality : A!=C", (*A) != (*C));
  cout << "SiconosMatrixTest >>> testEquality ...................................................... OK\n ";
}


void SiconosMatrixTest::testLinearSolve()
{
  SiconosMatrix X(10, 10);
  SiconosMatrix Xbefore = X;
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : Xbefore==X", Xbefore == X);
  X = (*A).linearSolve((*B));
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : A==A", (*A) == (*A));
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : B==B", (*B) == (*B));
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : X!=A", X != (*A));
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : X!=B", X != (*B));
  CPPUNIT_ASSERT_MESSAGE("testLinearSolve : Xbefore!=X", Xbefore != X);
  cout << "SiconosMatrixTest >>> testLinearSolve ................................................... OK\n ";
}

void SiconosMatrixTest::testReadWriteBinary()
{
  SiconosMatrix X = (*B);
  (*A).write("totoMat.bin", "binary");
  X.read("totoMat.bin", "binary");

  CPPUNIT_ASSERT_MESSAGE("testReadWriteBinary : A == X ", (*A) == X);
  SiconosMatrix Y("totoMat.bin", false);
  CPPUNIT_ASSERT_MESSAGE("testReadWriteBinary : A == Y ", *A == Y);
  cout << "SiconosMatrixTest >>> testReadWriteBinary ............................................... OK\n ";
}

void SiconosMatrixTest::testReadWriteAscii()
{
  SiconosMatrix X = (*B);
  (*A).write("totoMat.ascii", "ascii");
  X.read("totoMat.ascii", "ascii");

  CPPUNIT_ASSERT_MESSAGE("testReadWriteAscii : A == X ", (*A) == X);
  SiconosMatrix Y("totoMat.ascii", true);
  CPPUNIT_ASSERT_MESSAGE("testReadWriteAscii : A == Y ", (*A) == Y);
  cout << "SiconosMatrixTest >>> testReadWriteAscii ................................................ OK\n ";
}

void SiconosMatrixTest::testAffectation()
{
  SiconosMatrix X = (*A);
  SiconosMatrix Y = *A;
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A == X ", (*A) == X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : &A!=&X ", A != &X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A == Y ", (*A) == Y);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : &A!=&Y ", A != &Y);
  (*A)(1, 1) = (*A)(1, 1) * 2;
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A != X ", (*A) != X);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : A != Y ", (*A) != Y);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : X == Y ", X == Y);

  cout << "SiconosMatrixTest >>> testAffectation ................................................... OK\n ";
}

void SiconosMatrixTest::testOperator()
{
  cout << endl;
  (*B) = (*A) + (*A);
  (*B) = (*A) - (*A);
  (*B) = (*A) * 2;
  (*B) = 2 * (*A);
  (*B) = (*A) / 2;
  (*B) = (*A) ^ 0;
  (*B) = (*A) ^ 5;
  (*B) = (*A).multTranspose((*A));
  (*B) = (*A).multTranspose((*B));

  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A == A*2 ", (*A) + (*A) == (*A) * 2.0);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A == 2*A ", (*A) + (*A) == 2 * (*A));
  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A+A+A+A == A*5 ", (*A) + (*A) + (*A) + (*A) + (*A) == (*A) * 5.0);
  CPPUNIT_ASSERT_MESSAGE("testOperator : A+A+A+A+A == 5.0*A ", (*A) + (*A) + (*A) + (*A) + (*A) == 5.0 * (*A));
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^1 == A ", ((*A) ^ 1) == (*A));
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^2 == A*A ", ((*A) ^ 2) == (*A) * (*A));
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^3 == A*A*A ", ((*A) ^ 3) == (*A) * (*A) * (*A));
  CPPUNIT_ASSERT_MESSAGE("testOperator : A^10 == A*A*A*A*A*A*A*A*A*A ", ((*A) ^ 10) == (*A) * (*A) * (*A) * (*A) * (*A) * (*A) * (*A) * (*A) * (*A) * (*A));

  cout << "SiconosMatrixTest >>> testOperator ...................................................... OK\n ";
}

void SiconosMatrixTest::testGetRow()
{
  SiconosMatrix *X = new SiconosMatrix(*A);
  SimpleVector *V2 = new SimpleVector(10);
  X->getRow(2, *V2);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : X == A ", *X == (*A));
  X->getRow(3, *SV);
  SimpleVector * V3 = new SimpleVector(10);
  X->getRow(3, *V3);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : V2 != V3 ", *V2 != *V3);
  CPPUNIT_ASSERT_MESSAGE("testAffectation : V3 == SV ", *V3 == *SV);
  delete V2;
  delete V3;
  delete X;
  cout << "SiconosMatrixTest >>> testGetRow ........................................................ OK\n ";
}

void SiconosMatrixTest::testSizeException()
{
  *A + (*B);
}

void SiconosMatrixTest::testConstructorException()
{
  LaVectorDouble LVD(10);
  for (int i = 0; i < LVD.size(); i++)
    LVD(i) = i;
  SiconosMatrix X(LVD, 4, 2);
}

void SiconosMatrixTest::testBlockMatrixCopy1()
{
  SiconosMatrix Aa(1, 2), Bb(3, 3), Cc(3, 3);
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
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy1 : B.blockMatrixCopy( A, x, y) ", Bb == Cc);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy1 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy2()
{
  SiconosMatrix Aa(1, 1), Bb(1, 1), Cc(1, 1);
  Aa(0, 0) = 1.0;
  Bb(0, 0) = 0.0;
  Cc(0, 0) = 1.0;

  Bb.blockMatrixCopy(Aa, 0, 0);
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy2 : B.blockMatrixCopy( A, x, y) ", Bb == Cc);
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy2 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy3()
{
  SiconosMatrix * Aa = new SiconosMatrix(2, 2);
  SiconosMatrix * Bb = new SiconosMatrix(4, 4);
  SiconosMatrix * Cc = new SiconosMatrix(4, 4);

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
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy3 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc));
  delete Aa;
  delete Bb;
  delete Cc;
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy3 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopy4()
{
  SiconosMatrix *Aa = new SiconosMatrix(2, 2), *Bb = new SiconosMatrix(4, 4), *Cc = new SiconosMatrix(4, 4);
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
  CPPUNIT_ASSERT_MESSAGE("testBlockMatrixCopy4 : B.blockMatrixCopy( A, x, y) ", (*Bb) == (*Cc));
  delete Aa;
  delete Bb;
  delete Cc;
  cout << "SiconosMatrixTest >>> testBlockMatrixCopy4 ........................................................ OK\n ";
}

void SiconosMatrixTest::testBlockMatrixCopyException1()
{
  SiconosMatrix *Aa = new SiconosMatrix(2, 2), *Bb = new SiconosMatrix(4, 4);
  (*Bb).blockMatrixCopy((*Aa), 3, 3);
  delete Aa;
  delete Bb;
}

void SiconosMatrixTest::testBlockMatrixCopyException2()
{
  SiconosMatrix *Aa = new SiconosMatrix(6, 6), *Bb = new SiconosMatrix(4, 4);
  (*Bb).blockMatrixCopy((*Aa), 0, 0);
  delete Aa;
  delete Bb;
}
